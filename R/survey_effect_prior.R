### Experiment with prior on survey parameter
## Goal is to reduce bias in survey effects that may impact abundance indices



library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


# downscale data and predictive grid
dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  filter(species == "pink")


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
  max.edge = c(1, 5) * 500,
  cutoff = 20,
  offset = c(20, 200)
) 
spde <- make_mesh(
  dat,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 

# specify missing season-year levels
missing_surveys <- year_season_key$ys_index[!(year_season_key$ys_index %in% 
                                                dat$ys_index)]

# fitting anisotropic models to pink data leads to convergence issues; adjust
# accordingly
dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest() %>%
  mutate(
    anisotropy = ifelse(species == "pink", FALSE, TRUE)
  )


# use N(0, 1) prior for survey effect
beta_mu_priors <- rep(NA, times = 7)
beta_mu_priors[5] <- 0
beta_sd_priors <- rep(NA, times = 7)
beta_sd_priors[5] <- 1

fit <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth
  ,
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
  priors = sdmTMBpriors(
    # location = vector of means; scale = vector of standard deviations:
    b = normal(location = beta_mu_priors, scale = beta_sd_priors),
  ),
  silent = FALSE
)

# with ar1 term
fit_ar1 <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth
  ,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f + year_f,
  time_varying = ~ 1,
  time_varying_type = "ar1",
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


# constrained RW SD
fit_fix_tauV <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
    scale_depth
  ,
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
    start = list(
      ln_tau_V = matrix(log(0.1), nrow(year_season_key), 1)
    ),
    map = list(
      ln_tau_V = factor(rep(NA, times = nrow(year_season_key))),
      # 1 per season, fix all years to same value
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
      )
    )
  ),
  silent = FALSE
)


# original fit
dat_tbl_original <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_final.rds")
)
fit_orig <- dat_tbl_original$fit[[4]]

# par ests
fit_list <- list(fit, fit_orig, fit_ar1, fit_fix_tauV)
names_fits <- c("prior", "no prior", "no prior - AR1", "no prior - fix tauV")

purrr::map2(
  fit_list, names_fits,
  ~ tidy(.x, effects = "fixed", conf.int = T) %>% 
    mutate(model = .y)
) %>% 
  bind_rows() %>% 
  filter(term == "survey_fipes")

tidy(fit, effects = "fixed", conf.int = T)
