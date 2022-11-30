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
    even_year = ifelse(is.even(year), TRUE, FALSE)
  ) %>% 
  droplevels() 


# keep early years in dataset for now but exclude 2022 and consider dropping 
# years with v. limited summer data
dat_in <- dat %>% 
  filter(!year < 1998,
         !year == "2022",
         !bath_depth_mean_m < 0)
dat_pink <- dat_in %>% filter(species == "PINK")
dat_pink_odd <-  dat_in %>% filter(species == "PINK" & even_year == FALSE)
dat_pink_even <-  dat_in %>% filter(species == "PINK" & even_year == TRUE)

# make meshes
dat_list <- list(dat_pink, dat_pink_odd, dat_pink_even)
mesh_list <- purrr::map(
  dat_list, 
  ~ {
    inla_mesh_raw <- INLA::inla.mesh.2d(
      loc = cbind(.x$utm_x_1000, .x$utm_y_1000),
      max.edge = c(1, 5) * 500,
      cutoff = 20,
      offset = c(20, 200)
    ) 
    make_mesh(.x, 
              c("utm_x_1000", "utm_y_1000"),
              mesh = inla_mesh_raw) 
  }
)

fit_list <- purrr::map2(
  dat_list, mesh_list, 
  ~ {
    sdmTMB(
      n_juv ~ 1 +
        as.factor(year) +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = .x$effort,
      data = .x,
      mesh = .y,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "iid",
      time = "year",
      anisotropy = TRUE,
      share_range = FALSE,
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  }
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

pink_st_list <- fit_list[c(1, 3)] #ignore odd due to convergence issues
names(pink_st_list) <- c("full", "even")

# new df 
week_dat <- data.frame(
  week = seq(0, 52, length.out = 100),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year = 2012L
) 

pred_full <-  predict(pink_st_list[[1]], 
                      newdata = week_dat, 
                      se_fit = T, 
                      re_form = NA) %>%
  mutate(model = "full")
pred_even <-  predict(pink_st_list[[2]], 
                      newdata = week_dat, se_fit = T, re_form = NA) %>%
  mutate(model = "even")

pink_preds <- list(pred_full, pred_even#, pred_odd
                   ) %>% 
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
  aes(week, exp_est, #ymin = exp_lo, ymax = exp_up,
             fill = model)
) +
  ggsidekick::theme_sleek() +
  ylab("Abundance Index") +
  geom_line() +
  # geom_ribbon(alpha = 0.3) +
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
