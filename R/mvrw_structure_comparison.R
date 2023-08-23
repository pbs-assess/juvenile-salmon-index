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


# RI and RW
fit_ri <- sdmTMB(
  n_juv ~ 0 + season_f,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  # spatial_varying = ~ 0 + season_f + year_f,
  # time_varying = ~ 1,
  # time_varying_type = "rw0",
  # time = "ys_index",
  # spatiotemporal = "off",
  anisotropy = FALSE,
  extra_time = missing_surveys,
  # control = sdmTMBcontrol(
  #   map = list(
  #     # 1 per season, fix all years to same value
  #     ln_tau_Z = factor(
  #       c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
  #     )
  #   )
  # ),
  silent = FALSE
)

fit_mvrfrw <- sdmTMB(
  n_juv ~ 0 + season_f,
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
  silent = FALSE
)

fit_mvrfrw2 <- sdmTMB(
  n_juv ~ 0 + season_f,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time_varying =  ~ 1,
  time_varying_type = "rw",
  time = "year",
  spatiotemporal = "rw",
  anisotropy = FALSE,
  mvrw_category = "season_f",
  silent = FALSE
)



# check pars
pars <- fit_mvrfrw$tmb_obj$env$parList()

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
  list(fit_mvrfrw, fit_mvrfrw2),
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

rbind(
  index_list[[1]] %>% mutate(model = "1"),
  index_list[[2]] %>% mutate(model = "2")) %>% 
  ggplot() +
  geom_pointrange(
    aes(x = year, y =  est, ymin = lwr, ymax = upr, fill = model),
    shape = 21, position = "jitter"
  ) +
  ggsidekick::theme_sleek()





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
