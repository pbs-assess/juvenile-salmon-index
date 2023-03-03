
library(tidyverse)
library(sdmTMB)
library(ggplot2)

# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3)
  ) %>% 
  filter(!species == "steelhead") 


sox_dat <- dat %>% 
  filter(species == "sockeye") %>% 
  mutate(cpue = n_juv / volume_km3) %>% 
  droplevels()



## coordinates should be common to all species
dat_coords <- sox_dat %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()

## use INLA mesh based on SA recommendations and model selection (see notes)
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(1, 5) * 500,
  cutoff = 20,
  offset = c(20, 200)
) 
spde <- make_mesh(
  sox_dat,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 


## fits 
fit_nb <- sdmTMB(
  n_juv ~ 0 +
    as.factor(year) +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth,
  offset = sox_dat$effort,
  data = sox_dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = TRUE,
  # share_range = FALSE,
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1
  ),
  silent = FALSE
)

fit_tw <- sdmTMB(
  cpue ~ 0 +
    as.factor(year) +
    dist_to_coast_km +
    s(week, bs = "cc", k = 5) +
    target_depth,
  # offset = sox_dat$effort,
  data = sox_dat,
  mesh = spde,
  family = sdmTMB::tweedie(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year",
  anisotropy = TRUE,
  # share_range = FALSE,
  knots = list(
    week = c(0, 52)
  ),
  control = sdmTMBcontrol(
    newton_loops = 1
  ),
  silent = FALSE
)


## verbose predictive grid (overkill for here)
grid_list <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  purrr::map(
    .,
    ~ {.x %>% 
        mutate(utm_x_1000 = X / 1000,
               utm_y_1000 = Y / 1000,
               dist_to_coast_km = shore_dist / 1000)}
  )

summer_grid <- grid_list$ipes_grid  

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(sox_dat$year),
  survey_f = "hss",
  week = 25
) %>%
  mutate(id = row_number()) %>%
  split(., .$id) %>%
  purrr::map(., function (x) {
    summer_grid %>% 
      mutate(
        year = x$year,
        survey_f = x$survey_f,
        target_depth = 0,
        week = x$week,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() 


# nb model index
sox_preds_nb <- predict(
  fit_nb,
  newdata = exp_grid,
  return_tmb_object = TRUE
)
sox_ind_nb <- get_index(sox_preds_nb, area = 4)


# tweedie model index
sox_preds_tw <- predict(
  fit_tw,
  newdata = exp_grid,
  return_tmb_object = TRUE
)
sox_ind_tw <- get_index(sox_preds_tw, area = 4)
sox_ind_tw2 <- sox_ind_tw %>% 
  mutate(
    est_catch = est * median(sox_dat$volume_km3),
    lwr_catch = lwr * median(sox_dat$volume_km3),
    upr_catch = upr * median(sox_dat$volume_km3)
  )

p1 <- ggplot(sox_ind_tw2, 
             aes(year, est_catch)) +
  # geom_pointrange(aes(ymin = lwr_catch, ymax = upr_catch), 
  #                 shape = 21) +
  geom_point() +
  labs(x = "Year", y = "Tweedie Abundance") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
p2<- ggplot(sox_ind_nb, aes(year, est)) +
  geom_point() +
  # geom_pointrange(aes(ymin = lwr, ymax = upr), 
  #                 shape = 21) +
  labs(x = "Year", y = "NegBin Abundance") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")
cowplot::plot_grid(p1, p2, nrow = 2)
