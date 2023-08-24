## Fit models w/ mvrfrw across different seasons instead of RW RIs

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


# downscale data and predictive grid
dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  filter(species == "chum")


# make key so that missing year-season combinations can be indexed 
# (only necessary for original RIs)
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
    year_f = as.factor(year),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night),
    # observation level intercept used for RI
    obs_f = 1:nrow(.) %>% as.factor()
  ) %>% 
  left_join(., year_season_key %>% select(ys_index, year_season_f)) %>% 
  droplevels()


## mesh shared among species
dat_coords <- dat %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(2, 10) * 500,
  cutoff = 30,
  offset = c(10, 50)
)  
spde <- make_mesh(
  dat,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 

# specify missing season-year levels
missing_surveys <- year_season_key$ys_index[!(year_season_key$ys_index %in% 
                                                dat$ys_index)]


# original model
fit_rw <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth,
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
  extra_time = missing_surveys,
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

# year variability entirely driven by upsilon
fit_mvrfrw <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  mvrw_category = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat$season_f)))
      )
    )
  ),
  silent = FALSE
)

# year RW intercepts (convergence issues)
# fit_mvrfrw2 <- sdmTMB(
#   n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
#     scale_depth,
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "off",
#   spatial_varying = ~ 0 + season_f,
#   time_varying =  ~ 1,
#   time_varying_type = "rw",
#   time = "year",
#   spatiotemporal = "rw",
#   anisotropy = FALSE,
#   mvrw_category = "season_f",
#   control = sdmTMBcontrol(
#     map = list(
#       ln_tau_Z = factor(
#         rep(1, times = length(unique(dat$season_f)))
#       )
#     )
#   ),
#   silent = FALSE
# )

# year FEs 
fit_mvrfrw3 <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth + year_f,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  # time_varying =  ~ 1,
  # time_varying_type = "rw",
  time = "year",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  mvrw_category = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat$season_f)))
      )
    )
  ),
  silent = FALSE
)

# add anisotropy
fit_mvrfrw4 <- update(fit_mvrfrw3, anisotropy = TRUE)

# save mod fits
saveRDS(list(fit_rw, fit_mvrfrw, fit_mvrfrw3), 
        here::here("data", "fits", "chum_fit_mvrfrw_list.rds"))
fit_list <- readRDS(here::here("data", "fits", "chum_fit_mvrfrw_list.rds"))


# check pars
pars <- fit_mvrfrw2$tmb_obj$env$parList()

pars$upsilon_stc
pars$mvrw_u
pars$mvrw_u

matplot(t(pars$mvrw_u), type = "p", pch = 21)


# compare indices
sp_scalar <- 1 * (13 / 1000)

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
  mutate(
    year_f = as.factor(year),
    id = row_number(),
    survey_f = survey_f %>% as.factor()
  ) %>% 
  split(., .$id)
#add unique years and seasons
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
        day_night = "DAY",
        scale_depth = (target_depth - mean(dat$target_depth)) /
          sd(dat$target_depth),
        scale_dist = (dist_to_coast_km - mean(dat$dist_to_coast_km)) /
          sd(dat$dist_to_coast_km)
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


# predictions for index (set to HSS reference values)
index_grid_hss <- index_grid %>% filter(survey_f == "hss")

ind_preds <- purrr::map(
  fit_list[2:3],
  ~ {
    predict(.x,
            newdata = index_grid_hss,
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)

index_list <- purrr::map(
  ind_preds,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)

# original index 
ind_orig <- readRDS(here::here("data", "fits", "season_index_list.rds"))
ind_rw <- ind_orig[[2]] %>% 
  left_join(., year_season_key, by = "ys_index") %>% 
  mutate(model = "original_rw") %>% 
  filter(season_f == "su") %>% 
  select(year, est, lwr, upr, log_est, se, model)
rbind(
  ind_rw,
  index_list[[1]] %>% mutate(model = "mvrfrw_only"),
  index_list[[2]] %>% mutate(model = "mvrfrw_year_fe")) %>% 
  ggplot() +
  geom_point(
    aes(x = year, y =  log_est, fill = model),
    shape = 21, position = position_dodge(0.6)
  ) +
  ggsidekick::theme_sleek()

sd(ind_rw$log_est)
sd(index_list[[1]]$log_est)
sd(index_list[[2]]$log_est)


# check diagnostics 
fit_mvrfrw3 <- fit_list[[3]]
sims <- simulate(fit_mvrfrw3, nsim = 50, newdata = dat)
pred_fixed <- fit_mvrfrw3$family$linkinv(
  predict(fit_mvrfrw3, newdata = dat)$est_non_rf
)
r_nb <- DHARMa::createDHARMa(
  simulatedResponse = sims,
  observedResponse = dat$n_juv,
  fittedPredictedResponse = pred_fixed
)
plot(r_nb)

(sum(dat$n_juv == 0) / length(dat$n_juv)) / (sum(sims == 0) / length(sims))

