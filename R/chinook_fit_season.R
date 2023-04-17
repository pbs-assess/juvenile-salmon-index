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


# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    # min_date = as.numeric(as.Date(min(date))) - 1,
    # cont_date = as.numeric(as.Date(date)) - min_date,
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    scale_season = scale(as.numeric(season_f))[ , 1],
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night)
  ) %>% 
  filter(species == "chinook") %>% 
  droplevels()

ggplot(dat %>% filter(season_f == "wi")) +
  geom_point(aes(x = utm_x, y = utm_y, size = n_juv), alpha = 0.4) +
  scale_size(trans = "sqrt") +
  facet_wrap(~year)
  

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
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f +
    (1 | season_year_f),
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
test8 <- sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time_varying = ~ 1 + season_f,
  time_varying_type = "rw0",
  time = "year",
  spatiotemporal = "ar1",
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

# test9 <- sdmTMB(
#   n_juv ~ 0 + year_f + season_f + s(cont_date) +
#     target_depth + day_night,
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "off",
#   spatial_varying = ~ 0 + season_f,
#   # time_varying = ~ 0 + season_f,
#   # time_varying_type = "ar1",
#   # time = "year",
#   # spatiotemporal = "ar1",
#   anisotropy = FALSE,
#   share_range = TRUE,
#   priors = sdmTMBpriors(
#     phi = halfnormal(0, 10),
#     matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
#     matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
#   ),
#   control = sdmTMBcontrol(
#     newton_loops = 1
#   ),
#   silent = FALSE
# )


dat_tbl <- tibble(
  model = c("season_ri", "season_fe"),
  fits = list(test6, test8)
)
saveRDS(dat_tbl, here::here("data", "fits", "chinook_sp_varying.rds"))

dat_tbl <- readRDS(here::here("data", "fits", "chinook_sp_varying.rds"))


dat_tbl$sims <- purrr::map(dat_tbl$fits, simulate, nsim = 50)
qq_list <- purrr::map2(dat_tbl$sims, dat_tbl$fits, sdmTMBextra::dharma_residuals)


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


## check spatial predictions ---------------------------------------------------

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
               utm_y_1000 = Y / 1000)}
  )
summer_grid <- grid_list$ipes_grid
fall_grid <- grid_list$wcvi_grid

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(dat$year),
  survey_f = unique(dat$survey_f),
  season_f = unique(dat$season_f)
) %>%
  mutate(
    year_f = as.factor(year),
    id = row_number(),
    # scale_week = ifelse(season_f == "wi", 1, -0.25),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    # remove fall ipes surveys (doesn't meet definition)
    !(season_f %in% c("wi") & survey_f == "ipes"),
    !season_f == "sp"#,
    # season_year_f %in% unique(dat$season_year_f)
  ) %>% 
  split(., .$id) %>%
  purrr::map(., function (x) {
    dum_grid <- if (x$season_f == "su") summer_grid else fall_grid
    
    dum_grid %>% 
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        season_f = x$season_f,
        season_year_f = x$season_year_f,
        # scale_week = x$scale_week,
        # scale_week2 = scale_week^2,
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
# remove season-year factor levels not present in original dataframe
levels(exp_grid$season_year_f)[!levels(exp_grid$season_year_f) %in%
                                 levels(dat$season_year_f)] <- NA

# constrain to high seas surveys
exp_grid_hss <- exp_grid %>% 
  filter(survey_f == "hss")

spatial_preds <- purrr::map(
  dat_tbl$fits, 
  ~ predict(.x, newdata = exp_grid, se_fit = FALSE, re_form = NULL)
)
sp_preds_all <- rbind(
  spatial_preds[[1]] %>% mutate(model = "ri"),
  spatial_preds[[2]] %>% mutate(model = "fe")
)
 
plot_map_raster <- function(dat, column = est) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c()
}
plot_map_raster(sp_preds_all %>% filter(year == "2012"), est) +
  facet_grid(model~season_f)

sp_preds_all %>% 
  filter(year == "2014",
         model == "fe") %>% 
  pivot_longer(cols = starts_with("zeta_s")) %>% 
  plot_map_raster(., value) +
  facet_wrap(~name)

plot_map_raster(sp_preds_all %>% 
                  filter(year %in% c("1999", "2004", "2009"), season_f == "wi"), 
                est) +
  facet_grid(year~model)


## indices from plausible models -----------------------------------------------

# check how predictions incorporate RIs
exp_grid_no_ri <- exp_grid_hss %>% 
  filter(season_f == "su", 
         is.na(season_year_f))
pred1 <- predict(dat_tbl$fits[[1]], newdata = exp_grid_no_ri[1:100, ],
                 se_fit = FALSE,
                 re_form_iid = NULL)
pred2 <- predict(dat_tbl$fits[[1]], newdata = exp_grid_no_ri[1:100, ],
                 se_fit = FALSE,
                 re_form_iid = NA)
pred3 <- predict(dat_tbl$fits[[1]], 
                 newdata = exp_grid_no_ri[1:100, ] %>% 
                   mutate(season_year_f = "su_2019"),
                 se_fit = FALSE,
                 re_form_iid = NULL)


# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m)
sp_scalar <- 1000^2 * 13

preds_ri <- predict(dat_tbl$fits[[1]], 
                 newdata = exp_grid_hss %>% 
                   filter(season_f == "su"), 
                 se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
preds_tv <- predict(dat_tbl$fits[[2]], 
                 newdata = exp_grid_hss %>% 
                   filter(season_f == "su"), 
                 se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
preds_ri_fa <- predict(dat_tbl$fits[[1]], 
                    newdata = exp_grid_hss %>% 
                      filter(season_f == "wi"), 
                    se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
preds_tv_fa <- predict(dat_tbl$fits[[2]], 
                    newdata = exp_grid_hss %>% 
                      filter(season_f == "wi"), 
                    se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)

index_list <- furrr::future_map(
  list(preds_ri, preds_tv, preds_ri_fa, preds_tv_fa),
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)

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

indices <- list(
  index_list[[1]] %>% 
    mutate(model = "ri",
           season = "summer"),
  index_list[[2]] %>% 
    mutate(model = "tv",
           season = "summer"),
  index_list[[3]] %>% 
    mutate(model = "ri",
           season = "fall"),
  index_list[[4]] %>% 
    mutate(model = "tv",
           season = "fall")
) %>% 
  bind_rows() %>% 
  mutate(
    survey = case_when(
      season == "fall" & year %in% fall_years ~ "sampled",
      season == "summer" & year %in% summer_years ~ "sampled",
      TRUE ~ "no survey"
    )
  )
saveRDS(indices, here::here("data", "fits", "ck_season_indices.rds"))

ggplot(indices, aes(year, log_est)) +
  geom_point(aes(colour = model, shape = survey),
             # shape = 21,
             position = position_dodge(width=0.3)) +
  # geom_pointrange(aes(ymin = lwr, ymax = upr, fill = model),
  #                 shape = 21, position = position_dodge(width=0.3)) +
  labs(x = "Year", y = "Log Abundance") +
  ggsidekick::theme_sleek() +
  facet_wrap(~season)

# check among season correlations, should be weaker for poly
cor(index_list[[1]]$log_est, index_list[[3]]$log_est)
cor(index_list[[2]]$log_est, index_list[[4]]$log_est)


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