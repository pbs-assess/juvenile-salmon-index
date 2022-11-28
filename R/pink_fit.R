## Impact of fitting pink even/odd separately?


library(tidyverse)
library(sdmTMB)
library(ggplot2)


is.even <- function(x) x %% 2 == 0


# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    week = lubridate::week(date),
    vessel = as.factor(vessel),
    # redefine bins for all species except chinook
    target_depth_bin = case_when(
      species %in% c("SOCKEYE", "PINK") & 
        target_depth_bin %in% c("30", "45", "60") ~ ">30",
      species %in% c("COHO", "CHUM") & 
        target_depth_bin %in% c("45", "60") ~ ">45",
      species == "CHINOOK" & target_depth_bin == "60" ~ ">60",
      TRUE ~ as.character(target_depth_bin)
    ),
    even_year = ifelse(is.even(year), TRUE, FALSE),
    species = case_when(
      species == "PINK" & even_year == TRUE ~ "PINK_EVEN",
      species == "PINK" & even_year == FALSE ~ "PINK_ODD",
      TRUE ~ species
    )
  ) %>% 
  droplevels() 


# keep early years in dataset for now but exclude 2022 and consider dropping 
# years with v. limited summer data
dat_in <- dat %>% 
  filter(!year == "2022",
         !bath_depth_mean_m < 0)


coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -137, ymin = 47, xmax = -121.25, ymax = 57) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# make meshes (differ among groups because off pooling of depth strata)
dat_tbl_temp <- dat_in %>% 
  filter(species %in% c("PINK_EVEN", "PINK_ODD")) %>% 
  group_by(species) %>% 
  group_nest() %>% 
  mutate(
    spde = purrr::map(data, make_mesh,  c("utm_x_1000", "utm_y_1000"), 
                      cutoff = 15, type = "kmeans"),
    bspde = purrr::map(spde, add_barrier_mesh,
                       coast_utm, range_fraction = 0.1,
                       # scaling = 1000 since UTMs were rescaled above
                       proj_scaling = 1000)
  )



pink_st_odd <- sdmTMB(
  n_juv ~ 1 +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    day_night +
    survey_f,
  offset = dat_tbl_temp$data[[2]]$effort,
  data = dat_tbl_temp$data[[2]],
  mesh = dat_tbl_temp$bspde[[2]],
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
  ),
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    nlminb_loops = 2,
    newton_loops = 1
  )
)

pink_st_even <- sdmTMB(
  n_juv ~ 1 +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    day_night +
    survey_f,
  offset = dat_tbl_temp$data[[1]]$effort,
  data = dat_tbl_temp$data[[1]],
  mesh = dat_tbl_temp$bspde[[1]],
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
  ),
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    nlminb_loops = 2,
    newton_loops = 1
  )
)


pink_st_even2 <- sdmTMB(
  n_juv ~ 1 +
    # dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    # day_night +
    survey_f,
  offset = dat_tbl_temp$data[[1]]$effort,
  data = dat_tbl_temp$data[[1]],
  mesh = dat_tbl_temp$bspde[[1]],
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = FALSE,
  # priors = sdmTMBpriors(
  #   matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
  # ),
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    nlminb_loops = 2,
    newton_loops = 1
  )
)


### CHECK SIMULATIONS ----------------------------------------------------------

dat_tbl <- readRDS(here::here("data", "fits", "st_mod_all_sp.rds"))

sp_dat <- dat_tbl %>% 
  select(species, data) %>% 
  unnest(cols = data) 

ggplot(sp_dat) +
  geom_boxplot(aes(x = as.factor(month), y = log(n_juv))) +
  facet_wrap(~species)

sim_dat <- dat_tbl$sims[[4]] %>%
  as.data.frame() %>% 
  rename_with(., ~ paste("n_juv", .x, sep = "_") %>% tolower()) 

pink_dat <- sp_dat %>% 
  filter(species == "PINK") %>% 
  cbind(., sim_dat[, 1:5]) %>%
  rename(n_juv_obs = n_juv) %>% 
  pivot_longer(
    cols = starts_with("n_juv"), 
    names_to = "n_juv",
    values_to = "count"
  )
  
ggplot(pink_dat) +
  geom_boxplot(aes(x = as.factor(month), y = log(count), fill = n_juv)) 


### CHECK PREDS ----------------------------------------------------------------

pink_st_list <- list(dat_tbl$st_mod[[4]],
                     pink_st_odd,
                     pink_st_even)
names(pink_st_list) <- c("full", "odd", "even")

# new df 
week_dat <- data.frame(
  week = seq(0, 52, length.out = 100),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year = 2012L
) 
week_dat_odd <- week_dat %>% mutate(year = 2011L)

pred_full <-  predict(dat_tbl$st_mod[[4]], 
                      newdata = week_dat, 
                      se_fit = T, 
                      re_form = NA) %>%
  mutate(model = "full")
pred_even <-  predict(pink_st_even, 
                      newdata = week_dat, se_fit = T, re_form = NA) %>%
  mutate(model = "even")
pred_odd <-  predict(pink_st_odd, 
                      newdata = week_dat_odd, se_fit = T, re_form = NA) %>%
  mutate(model = "odd")

pink_preds <- list(pred_full, pred_even, pred_odd) %>% 
  bind_rows() %>% 
  mutate(
    up = est + 1.96 * est_se,
    lo = est - 1.96 * est_se,
    exp_est = exp(est),
    exp_up = exp(up),
    exp_lo = exp(lo)
  ) %>% 
  glimpse()

ggplot(
  pink_preds,
  aes(week, exp_est, ymin = exp_lo, ymax = exp_up,
             fill = model)
) +
  ggsidekick::theme_sleek() +
  ylab("Abundance Index") +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) 


### ESTIMATES FROM NON-SPATIAL -------------------------------------------------

pink_dat <- dat_tbl$data[[4]]

mgcv_f <- mgcv::gam(
  n_juv ~ 1 +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    day_night +
    survey_f,
  data = pink_dat,
  family = mgcv::nb()
)
sdm_f <- sdmTMB(
  n_juv ~ 1 +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    day_night +
    survey_f,
  data = pink_dat,
  family = sdmTMB::nbinom2(),
  spatial = "off"
)


### REFIT W/ FEs ---------------------------------------------------------------

pink_dat <- dat_in %>% 
  filter(grepl("PINK", species))
pink_spde <- make_mesh(pink_dat,  
                       c("utm_x_1000", "utm_y_1000"),
                       type = "kmeans",
                       n_knots = 200)

pink_spde2 <- make_mesh(pink_dat,  
                       c("utm_x_1000", "utm_y_1000"),
                       type = "cutoff",
                       cutoff = 30)

pink_st_new <- sdmTMB(
  n_juv ~ 1 +
    as.factor(year) +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    day_night +
    survey_f,
  offset = pink_dat$effort,
  data = pink_dat,
  mesh = pink_spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = TRUE,
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1
  )
)

pink_st_new2 <- sdmTMB(
  n_juv ~ 1 +
    as.factor(year) +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth +
    day_night +
    survey_f,
  offset = pink_dat$effort,
  data = pink_dat,
  mesh = pink_spde2,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = TRUE,
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1
  )
)
