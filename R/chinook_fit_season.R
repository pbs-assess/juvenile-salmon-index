### Chinook only fit
## Experiment w/ spatiotemporal models allowing for spatially varying seasonal 
## effects
## Apr 6, 2023

# alternative branch
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "zeta-intercept")

library(tidyverse)
library(sdmTMB)


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) 


# make key so that missing year-season combinations can be indexed
year_season_key <- expand.grid(
  year = unique(dat_in$year) %>% sort(),
  season_f = unique(dat_in$season_f) %>% sort()
) %>% 
  arrange(year, season_f) %>% 
  mutate(
    # make consecutive index for model fitting
    ys_index = seq(1, nrow(.), by = 1),
    year_season_f = paste(year, season_f, sep = "_") %>% as.factor()
  )


# downscale data and predictive grid
dat <- dat_in %>% 
  mutate(
    # min_date = as.numeric(as.Date(min(date))) - 1,
    # cont_date = as.numeric(as.Date(date)) - min_date,
    # adjust year for spring surveys since they are catching previous cohort
    # adj_year = ifelse(season_f == "sp", year - 1, year),
    year_f = as.factor(year),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    scale_season = scale(as.numeric(season_f))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night)
  ) %>% 
  left_join(., year_season_key %>% select(ys_index, year_season_f)) %>% 
  filter(species == "chinook") %>% 
  droplevels()



# ggplot(dat %>% filter(season_f == "wi")) +
#   geom_point(aes(x = utm_x, y = utm_y, size = n_juv), alpha = 0.4) +
#   scale_size(trans = "sqrt") +
#   facet_wrap(~year)
  

## mesh
dat_coords <- dat %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(1, 5) * 500,
  cutoff = 20,
  offset = c(20, 200)
) 
spde <- make_mesh(
  dat,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 


# fit season as spatially varying effect
test1 <-  sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 1 + season_f,
  time = "year",
  spatiotemporal = "ar1",
  anisotropy = TRUE,  
  share_range = TRUE, 
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)
test6a <- sdmTMB(
  n_juv ~ 0 + year_f + season_f +# target_depth + day_night + survey_f +
    #dist_to_coast_km + 
    (1 | year_season_f),
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  # time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1
    ),
  silent = FALSE
)
test6 <- update(test6a, time = "year", spatiotemporal = "ar1")
test8a <- sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f + 
    dist_to_coast_km,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time_varying = ~ 1 + season_f,
  time_varying_type = "rw0",
  time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    start = list(
      ln_tau_V = matrix(log(c(0.1, 0.1, 0.1)), nrow = 3, ncol = 1)
    ),
    map = list(
      ln_tau_V = factor(c(NA, NA, NA))
    )
  ),
  silent = FALSE
)
test8 <- update(test8a, spatiotemporal = "ar1")

test9a <- sdmTMB(
  n_juv ~ 0 + season_f + target_depth + day_night + survey_f + 
    dist_to_coast_km,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f + year_f,
  time_varying = ~ 1,
  time_varying_type = "rw0",
  time = "ys_index",
  spatiotemporal = "off",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    map = list(
      # 1 per season, fix all years to same value
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
      )
    )
  ),
  silent = FALSE
)

missing_surveys <- year_season_key$ys_index[!(year_season_key$ys_index %in% 
                                                dat$ys_index)]
test9 <- update(test9a, extra_time = missing_surveys)

# purrr::map(list(test8a, test8b), tidy, "ran_pars")
# 
# 
# dat_tbl <- tibble(
#   model = c("season_ri", "season_fe"),
#   fits = list(test6, test8)
# )
# saveRDS(dat_tbl, here::here("data", "fits", "chinook_sp_varying.rds"))
# 
# dat_tbl <- readRDS(here::here("data", "fits", "chinook_sp_varying.rds"))


sims <- simulate(test9a, nsim = 50)
sdmTMBextra::dharma_residuals(sims, test9a)


dat_tbl$resids <- purrr::map(dat_tbl$fits, function(x) {
  dat$resids <- residuals(x)
  return(dat)
}) 

sum(dat$n_juv == 0) / length(dat$n_juv)
sum(dat_tbl$sims[[1]] == 0)/length(dat_tbl$sims[[1]])
sum(dat_tbl$sims[[2]] == 0)/length(dat_tbl$sims[[2]])


# par estimates 
purrr::map(dat_tbl$fits, tidy, "ran_pars")
tt <- purrr::map(dat_tbl$fits, 
           ~ tidy(.x, "fixed") %>% print(n = 30))


## SPATIAL GRID  ---------------------------------------------------------------

# shape file for coastline
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -129, ymin = 48.25, xmax = -124, ymax = 51.2) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# predictive grid
grid_list <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  purrr::map(
    .,
    ~ {.x %>% 
        mutate(utm_x_1000 = X / 1000,
               utm_y_1000 = Y / 1000,
               dist_to_coast_km = shore_dist / 1000)}
  )

summer_grid <- grid_list$ipes_grid
fall_grid <- grid_list$wcvi_grid

pred_grid_list <- rbind(
  year_season_key %>% mutate(survey_f = "hss"),
  year_season_key %>% mutate(survey_f = "ipes")
) %>% 
#   expand.grid(
#   year = unique(dat$year),
#   survey_f = unique(dat$survey_f),
#   season_f = unique(dat$season_f)
# ) %>%
  mutate(
    year_f = as.factor(year),
    id = row_number(),
    survey_f = survey_f %>% as.factor()
  ) %>% 
  split(., .$id)


## SPATIAL PREDICTIONS ---------------------------------------------------------

spatial_grid <- pred_grid_list %>%
  purrr::map(., function (x) {
    fall_grid %>% 
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        season_f = x$season_f,
        ys_index = x$ys_index,
        target_depth = 0,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  filter(
    !survey_f == "ipes"
  ) %>% 
  select(-c(depth, slope, shore_dist)) 

spatial_preds <- predict(test9, newdata = spatial_grid, 
                         se_fit = FALSE, re_form = NULL)

 
plot_map_raster <- function(dat, column = est) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c()
}

plot_map_raster(spatial_preds %>% filter(year == "2012"), est) +
  facet_wrap(~season_f)

spatial_preds %>% 
  filter(year == "2014") %>% 
  pivot_longer(cols = starts_with("zeta_s_season")) %>% 
  plot_map_raster(., value) +
  facet_wrap(~name) +
  scale_fill_gradient2()

spatial_preds %>% 
  filter(year == "2014") %>% 
  pivot_longer(cols = starts_with("zeta_s_year")) %>% 
  plot_map_raster(., value) +
  facet_wrap(~name) +
  scale_fill_gradient2()

year_seq <- c("1999", "2004", "2009", "2014", "2019")
plot_map_raster(spatial_preds %>% 
                  filter(year %in% year_seq), 
                est) +
  facet_grid(year~season_f)


## INDICES ---------------------------------------------------------------------

# add unique years and seasons
index_grid <- pred_grid_list %>%
  purrr::map(., function (x) {
    dum_grid <- if (x$season_f == "su") summer_grid else fall_grid
    
    dum_grid %>% 
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        season_f = x$season_f,
        ys_index = x$ys_index,
        target_depth = 0,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  mutate(
    fake_survey = case_when(
      year > 2016 & survey_f == "hss" & season_f == "su" ~ "1",
      year < 2017 & survey_f == "ipes" & season_f == "su" ~ "1",
      TRUE ~ "0")
  ) %>% 
  select(-c(depth, slope, shore_dist))

index_grid_hss <- index_grid %>% filter(survey_f == "hss")


# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m)
sp_scalar <- 1000^2 * 13

preds <- predict(test9, 
                 newdata = index_grid_hss, 
                 se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
index_full <- get_index(preds, area = sp_scalar, bias_correct = TRUE)


# define real years
summer_years <- dat %>% 
  filter(season_f == "su", 
         # remove 2021 since only partial survey
         !year_f == "2021") %>%
  pull(year) %>% 
  unique()
fall_years <- dat %>% 
  filter(season_f == "wi",
         # remove 2020 since most of survey was outside of grid
         !year_f == "2020") %>%
  pull(year) %>% 
  unique()

indices <- index_full %>%
  left_join(., year_season_key, by = "ys_index") %>% 
  filter(!season_f == "sp") %>% 
  mutate(
    season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
    survey = case_when(
      season == "fall" & year %in% fall_years ~ "sampled",
      season == "summer" & year %in% summer_years ~ "sampled",
      TRUE ~ "no survey"
    )
  )
saveRDS(indices, here::here("data", "fits", "ck_season_indices.rds"))

ggplot(indices, aes(year, log_est)) +
  geom_point(aes(fill = survey),
             shape = 21,
             position = position_dodge(width=0.3)) +
  # geom_pointrange(aes(ymin = lwr, ymax = upr, fill = model),
  #                 shape = 21, position = position_dodge(width=0.3)) +
  labs(x = "Year", y = "Log Abundance") +
  ggsidekick::theme_sleek() +
  facet_wrap(~season)

# check among season correlations, should be weaker for poly
fall_dat <- indices %>% filter(season == "fall") %>% pull(log_est)
summ_dat <- indices %>% filter(season == "summer") %>% pull(log_est)
cor(fall_dat, summ_dat)


## compare non-st indices ------------------------------------------------------

test8a <- sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f + 
    dist_to_coast_km,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time_varying = ~ 1 + season_f,
  time_varying_type = "rw0",
  time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    start = list(
      ln_tau_V = matrix(log(c(0.1, 0.1, 0.1)), nrow = 3, ncol = 1)
    ),
    map = list(
      ln_tau_V = factor(c(NA, NA, NA))
    )
  ),
  silent = FALSE
)
test8b <- sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f + 
    dist_to_coast_km,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time_varying = ~ 1 + season_f,
  time_varying_type = "rw0",
  time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    start = list(
      ln_tau_V = matrix(log(c(0.3, 0.3, 0.3)), nrow = 3, ncol = 1)
    ),
    map = list(
      ln_tau_V = factor(c(NA, NA, NA))
    )
  ),
  silent = FALSE
)
purrr::map(list(test8a, test8b), tidy, "ran_pars")

preds_tvA <- predict(test8a, 
                    newdata = exp_grid_hss %>% 
                      filter(season_f == "su"), 
                    se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
preds_tvA_fa <- predict(test8a, 
                       newdata = exp_grid_hss %>% 
                         filter(season_f == "wi"), 
                       se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
preds_tvB <- predict(test8b, 
                     newdata = exp_grid_hss %>% 
                       filter(season_f == "su"), 
                     se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
preds_tvB_fa <- predict(test8b, 
                        newdata = exp_grid_hss %>% 
                          filter(season_f == "wi"), 
                        se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)

index_list <- furrr::future_map(
  list(preds_tvA,
       preds_tvA_fa,
       preds_tvB,
       preds_tvB_fa),
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
indices <- list(
  index_list[[1]] %>% 
    mutate(model = "A",
           season = "summer"),
  index_list[[3]] %>% 
    mutate(model = "B",
           season = "summer"),
  index_list[[2]] %>% 
    mutate(model = "A",
           season = "fall"),
  index_list[[4]] %>% 
    mutate(model = "B",
           season = "fall")
) %>% 
  bind_rows() %>% 
  mutate(
    survey = case_when(
      season == "fall" & year %in% fall_years ~ "sampled",
      season == "summer" & year %in% summer_years ~ "sampled",
      TRUE ~ "no survey"
    ),
    lwr_log = log_est - (1.96 * se),
    upr_log = log_est + (1.96 * se)
  )

ggplot(indices, aes(year, log_est)) +
  # geom_point(aes(colour = model, shape = survey),
  #            # shape = 21,
  #            position = position_dodge(width=0.3)) +
  geom_pointrange(aes(ymin = lwr_log, ymax = upr_log, fill = model,
                      shape = survey),
                  position = position_dodge(width=0.3)) +
  labs(x = "Year", y = "Log Abundance") +
  ggsidekick::theme_sleek() +
  facet_wrap(model~season)

cor(index_list[[1]]$log_est, index_list[[2]]$log_est)
cor(index_list[[3]]$log_est, index_list[[4]]$log_est)


## model comp ------------------------------------------------------------------


ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


cv_ri <- sdmTMB_cv(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f +
    (1 | season_year_f),
  offset = "effort",
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  # time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1
  ),
  silent = FALSE,
  k_folds = 5,
  use_initial_fit = TRUE
)

cv_tv <- sdmTMB_cv(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f,
  offset = "effort",
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time_varying = ~ 1 + season_f,
  time_varying_type = "rw0",
  time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,
  share_range = TRUE,
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1,
    start = list(
      ln_tau_V = matrix(log(c(0.1, 0.1, 0.1)), nrow = 3, ncol = 1)
    ),
    map = list(
      ln_tau_V = factor(c(NA, NA, NA))
    )
  ),
  silent = FALSE,
  k_folds = 5,
  use_initial_fit = TRUE
)

#> cv_ri$elpd
# [1] -0.4294949
# > cv_tv$elpd
# [1] -0.4294955