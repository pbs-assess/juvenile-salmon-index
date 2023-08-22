## Fit models w/ mvrw across different seasosn instead of RW RIs

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


# RI and RW
fit_ri <- sdmTMB(
  n_juv ~ 1,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  # spatial_varying = ~ 0 + season_f + year_f,
  time_varying = ~ 1,
  time_varying_type = "rw0",
  time = "ys_index",
  spatiotemporal = "off",
  anisotropy = FALSE,
  mvrw_category = NULL,
  extra_time = missing_surveys,
  control = sdmTMBcontrol(
    # map = list(
    #   # 1 per season, fix all years to same value
    #   ln_tau_Z = factor(
    #     c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
    #   )
    # )
  ),
  silent = FALSE
)

fit_mvrw <- sdmTMB(
  n_juv ~ 1,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  time = "year",
  spatiotemporal = "off",
  anisotropy = FALSE,
  mvrw_category = "season_f",
  silent = FALSE
)

fit_mvrw_sd <- sdmTMB(
  n_juv ~ 1,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  time = "year",
  spatiotemporal = "off",
  anisotropy = FALSE,
  mvrw_category = "season_f",
  control = sdmTMBcontrol(
    map = list(
      mvrw_logsds = factor(rep(1, times = length(unique(dat$season_f))))
    )
  ),
  silent = FALSE
)

fit_mvrw_fix_sumsd <- sdmTMB(
  n_juv ~ 1,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  time = "year",
  spatiotemporal = "off",
  anisotropy = FALSE,
  mvrw_category = "season_f",
  control = sdmTMBcontrol(
    start = list(
      mvrw_logsds = matrix(log(0.5), length(unique(dat$season_f)), 1)
    ),
    map = list(
      mvrw_logsds = factor(c(1, NA, 2))
    )
  ),
  silent = FALSE
)


# check pars
r <- fit_mvrw$tmb_obj$report(fit_mvrw$tmb_obj$env$last.par.best)
r_sd <- fit_mvrw_sd$tmb_obj$report(fit_mvrw_sd$tmb_obj$env$last.par.best)
r_fix <- fit_mvrw_fix_sumsd$tmb_obj$report(fit_mvrw_fix_sumsd$tmb_obj$env$last.par.best)
r$mvrw_Sigma
r_sd$mvrw_Sigma
r_fix$mvrw_Sigma

pars <- fit_mvrw$tmb_obj$env$parList()

pars$mvrw_rho
pars$mvrw_logsds
pars$mvrw_u

matplot(t(pars$mvrw_u), type = "p", pch = 21)

# make preds
nd <- year_season_key
preds_ri <- predict(fit_ri, newdata = nd) %>% 
  mutate(model = "rw")
preds_mvrw <- predict(
  fit_mvrw, 
  newdata = nd
) %>% 
  mutate(model = "mvrw")
preds_mvrw_sd <- predict(
  fit_mvrw_sd, 
  newdata = nd
) %>% 
  mutate(model = "mvrw_single_sd")
preds_mvrw_sd_fix <- predict(
  fit_mvrw_fix_sumsd, 
  newdata = nd
) %>% 
  mutate(model = "mvrw_summer_sd_fix")


preds <- list(preds_ri, preds_mvrw, preds_mvrw_sd, preds_mvrw_sd_fix) %>% 
  bind_rows() %>% 
  mutate(survey_missing = ifelse(ys_index %in% missing_surveys, "yes", "no"))
shape_pal <- c(1, 21)
names(shape_pal) <- c("yes", "no")

ggplot(preds) +
  geom_point(aes(x = year, y = est, colour = model, fill = model, shape = survey_missing)) +
  scale_shape_manual(values = shape_pal) +
  facet_wrap(~season_f) +
  ggsidekick::theme_sleek() +
  guides(
    shape = guide_legend(
      override.aes = list(fill = "grey")
    )
  )


## fit full model 
fit1 <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + #survey_f + 
    scale_dist +
    scale_depth,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f + year_f,
  time = "year",
  spatiotemporal = "off",
  anisotropy = TRUE,
  mvrw_category = "season_f",
  control = sdmTMBcontrol(
    map = list(
      # mvrw_logsds = factor(rep(1, times = length(unique(dat$season_f)))),
      # 1 per season, fix all years to same value
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
      )
    )
  ),
  silent = FALSE
)

mvrw_mu_priors <- rep(0, length = length(unique(dat$season_f)))
mvrw_sd_priors <- rep(3, length = length(unique(dat$season_f)))

fit1 <- sdmTMB(
  n_juv ~ 1,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f + year_f,
  time = "year",
  spatiotemporal = "off",
  # anisotropy = TRUE,
  mvrw_category = "season_f",
  # priors = sdmTMBpriors(
  #   # location = vector of means; scale = vector of standard deviations:
  #   mvrw_rho = normal(location = mvrw_mu_priors, scale = mvrw_sd_priors),
  # ),
  control = sdmTMBcontrol(
    map = list(
      # mvrw_logsds = factor(rep(1, times = length(unique(dat$season_f)))),
      # 1 per season, fix all years to same value
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
      )
    )
  ),
  silent = FALSE
)

fit2 <- sdmTMB(
  n_juv ~ 0 + season_f + day_night + #survey_f + 
    scale_dist +
    scale_depth,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f + year_f,
  time = "year",
  spatiotemporal = "off",
  anisotropy = TRUE,
  mvrw_category = "season_f",
  control = sdmTMBcontrol(
    map = list(
      mvrw_logsds = factor(rep(1, times = length(unique(dat$season_f)))),
      # 1 per season, fix all years to same value
      ln_tau_Z = factor(
        c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
      )
    )
  ),
  silent = FALSE
)

beta_mu_priors <- rep(NA, times = 7)
beta_mu_priors[5] <- 0
beta_sd_priors <- rep(NA, times = 7)
beta_sd_priors[5] <- 1
priors = sdmTMBpriors(
  # location = vector of means; scale = vector of standard deviations:
  b = normal(location = beta_mu_priors, scale = beta_sd_priors),
)


sims <- simulate(fit_mvrw_full, nsim = 50, newdata = dat)
pred_fixed <- fit_mvrw_full$family$linkinv(
  predict(fit_mvrw_full, newdata = dat)$est_non_rf
)
r_nb <- DHARMa::createDHARMa(
  simulatedResponse = sims,
  observedResponse = dat$n_juv,
  fittedPredictedResponse = pred_fixed
)
plot(r_nb)

(sum(dat$n_juv == 0) / length(dat$n_juv)) / (sum(sims == 0) / length(sims))


r_full <- fit_mvrw_full$tmb_obj$report(fit_mvrw_full$tmb_obj$env$last.par.best)
r_full$mvrw_Sigma

pars <- fit1$tmb_obj$env$parList()

pars$mvrw_rho
pars$mvrw_logsds
pars$mvrw_u
