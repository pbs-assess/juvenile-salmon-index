### Chinook only fit
## Experiment w/ spatiotemporal models allowing for spatially varying seasonal 
## effects
## Apr 6, 2023



library(tidyverse)
library(sdmTMB)
library(ggplot2)


# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    scale_season = scale(as.numeric(season_f))[ , 1]
  ) %>% 
  filter(species == "chinook") %>% 
  droplevels()


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
  spatial = "on",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "ar1",
  anisotropy = FALSE,  
  share_range = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)
test2 <- update(test1, spatiotemporal = "iid")
test3 <- update(
  test2, 
  priors = sdmTMBpriors(
    phi = halfnormal(0, 10),
    matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
  )
)
test4 <- update(test1, anisotropy = TRUE)

dat_out <- tibble(
  name = seq(1, 4, by = 1),
  fit = list(test1, test2, test3, test4)
)
saveRDS(dat_out, here::here("data", "fits", "chinook_season_sp.rds"))


# alternative to season-specific random fields, fit using multidimensional 
# group level-smoother
test5 <- sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night +
    t2(utm_x_1000, utm_y_1000, season_f, m=2, bs = c("tp", "tp", "re")),
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  # spatial = "on",
  time = "year",
  spatiotemporal = "iid",
  anisotropy = FALSE,  
  share_range = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)




## cross validation of top models

library(future)
plan(multisession)

set.seed(123)
cv5 <- sdmTMB_cv(
  n_juv ~ target_depth + day_night + survey_f,
  offset = "effort",
  data = dum,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  time = "year",
  anisotropy = FALSE,  
  share_range = TRUE,
  extra_time = extra_time,
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE,
  k_folds = 5,
  use_initial_fit = TRUE
)

set.seed(123)
cv6 <- sdmTMB_cv(
  n_juv ~ target_depth + day_night + survey_f,
  offset = "effort",
  data = dum,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "rw",
  time = "year",
  anisotropy = FALSE,  
  share_range = TRUE,
  extra_time = extra_time,
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE,
  k_folds = 5,
  use_initial_fit = TRUE
)

set.seed(123)
cv14 <- sdmTMB_cv(
  n_juv ~ 0 + year_f + target_depth + day_night + survey_f + s(week, k = 4),
  offset = "effort",
  data = dum,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  anisotropy = TRUE,  
  share_range = TRUE,
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE,
  k_folds = 5,
  use_initial_fit = TRUE
)

purrr::map(
  list(cv5, cv6, cv14), ~ .x$elpd
)
purrr::map(
  list(cv5, cv6, cv14), ~ mean(.x$fold_elpd)
)


## indices from plausible models -----------------------------------------------

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

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(dat$year),
  survey_f = unique(dum$survey_f)
) %>%
  mutate(id = row_number()) %>%
  split(., .$id) %>%
  purrr::map(., function (x) {
    summer_grid %>% 
      mutate(
        year = x$year,
        year_f = as.factor(x$year),
        survey_f = x$survey_f,
        target_depth = 0,
        week = 25,
        month_f = "7",
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  mutate(
    fake_survey = case_when(
      year > 2016 & survey_f == "hss" ~ "1",
      year < 2017 & survey_f == "ipes" ~ "1",
      TRUE ~ "0")
  )
sp_scalar <- 1000^2 * 13



pred5 <- predict(test5, 
                 newdata = exp_grid %>% 
                   filter(survey_f == "hss"), 
                 se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
pred6 <- predict(test6, 
                 newdata = exp_grid %>% 
                   filter(survey_f == "hss"), 
                 se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
pred13 <- predict(test13, 
                  newdata =  exp_grid %>% 
                    filter(survey_f == "hss",
                           year %in% dum$year), 
                  se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)

ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)

index_list <- purrr::map(
  list(pred5, pred6), 
  get_index, 
  area = sp_scalar, 
  bias_correct = TRUE
)
pp <- get_index(pred13, area = sp_scalar,
                bias_correct = TRUE)

indices <- list(
  index_list[[1]] %>% 
    mutate(model = "st_ar1"),
  index_list[[2]] %>% 
    mutate(model = "st_rw")#,
  # index_list[[3]] %>% 
  #   mutate(model = "sp")
) %>% 
  bind_rows() %>% 
  filter(!year %in% extra_time)


ggplot(indices, aes(year, log_est)) +
  geom_point(aes(fill = model),
             shape = 21, position = position_dodge(width=0.3)) +
  # geom_pointrange(aes(ymin = lwr, ymax = upr, fill = model),
  #                 shape = 21, position = position_dodge(width=0.3)) +
  labs(x = "Year", y = "Log Abundance") +
  ggsidekick::theme_sleek() 

