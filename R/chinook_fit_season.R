### Chinook only fit
## Experiment w/ spatiotemporal models allowing for spatially varying seasonal 
## effects
## Apr 6, 2023

library(tidyverse)
library(sdmTMB)


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    scale_season = scale(as.numeric(season_f))[ , 1],
    scale_week = scale(as.numeric(week))[ , 1],
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
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
test0 <-  sdmTMB(
  n_juv ~ 0 + year_f + season_f + target_depth + day_night,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  # spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "iid",
  anisotropy = TRUE,  
  share_range = TRUE, 
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
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
  anisotropy = TRUE,  
  share_range = TRUE, 
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)
# test2 <- update(test1, spatiotemporal = "iid")
# test3 <- update(test1, anisotropy = FALSE)
test4 <- sdmTMB(
  n_juv ~ 0 + year_f + target_depth + day_night,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatial_varying = ~ 0 + scale_week,
  time_varying = ~ 0 + poly(scale_week, 2),
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
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)
# test5 <- sdmTMB(
#   n_juv ~ 0 + year_f +
#     # s(scale_week, k = 4, bs = "cc", m = 2) +
#     s(week, k = 4, bs = "cc", by = year_f, m = 2) +
#     target_depth + day_night,
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatial_varying = ~ 0 + season_f,
#   time = "year",
#   spatiotemporal = "iid",
#   anisotropy = FALSE,
#   share_range = TRUE,
#   priors = sdmTMBpriors(
#     phi = halfnormal(0, 10),
#     matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
#     matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
#   ),
#   control = sdmTMBcontrol(
#     newton_loops = 1#,
#     # nlminb_loops = 2
#   ),
#   silent = FALSE
# )
# test6 <- sdmTMB(
#   n_juv ~ 0 + year_f + season_f + target_depth + day_night + 
#     (1 | season_year_f),
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatial_varying = ~ 0 + season_f,
#   time = "year",
#   spatiotemporal = "ar1",
#   anisotropy = FALSE,  
#   share_range = TRUE,
#   priors = sdmTMBpriors(
#     phi = halfnormal(0, 10),
#     matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
#     matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
#   ),
#   control = sdmTMBcontrol(
#     newton_loops = 1#,
#     # nlminb_loops = 2
#   ),
#   silent = FALSE
# )
# # as test1 but include time- and spatial-varying effect for seasons
# test7 <-  sdmTMB(
#   n_juv ~ 0 + year_f + season_f + target_depth + day_night,
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatial_varying = ~ season_f,
#   time_varying = ~ season_f,
#   time = "year",
#   spatiotemporal = "ar1",
#   anisotropy = FALSE,  
#   share_range = TRUE, 
#   priors = sdmTMBpriors(
#     phi = halfnormal(0, 10),
#     matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
#     matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
#   ),
#   control = sdmTMBcontrol(
#     newton_loops = 1#,
#     # nlminb_loops = 2
#   ),
#   silent = FALSE
# )

test8 <- sdmTMB(
  n_juv ~ target_depth + day_night + season_f +
    (1 | season_year_f),
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatial_varying = ~ 0 + season_f,
  # time_varying = ~ 1,
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
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)

# alternative to season-specific random fields, fit using multidimensional 
# group level-smoother
# DOES NOT CONVERGE
# test5 <- sdmTMB(
#   n_juv ~ 0 + year_f + season_f + target_depth + day_night +
#     t2(utm_x_1000, utm_y_1000, season_f, m=2, bs = c("tp", "tp", "re")),
#   offset = dat$effort,
#   data = dat,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   # spatial = "on",
#   time = "year",
#   spatiotemporal = "iid",
#   anisotropy = FALSE,  
#   share_range = TRUE,
#   control = sdmTMBcontrol(
#     newton_loops = 1#,
#     # nlminb_loops = 2
#   ),
#   silent = FALSE
# )


dat_tbl <- tibble(
  model = c("spatial season", "spatial poly", "spatial RI"),
  fits = list(test1, test4, test8)
)
saveRDS(dat_tbl, here::here("data", "fits", "chinook_sp_varying.rds"))

dat_tbl <- readRDS(here::here("data", "fits", "chinook_sp_varying.rds"))


dat_tbl$sims <- purrr::map(dat_tbl$fits, simulate, nsim = 50)
qq_list <- purrr::map2(dat_tbl$sims, dat_tbl$fits, sdmTMBextra::dharma_residuals)


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
  season_f = c("su", "wi")
) %>%
  filter(
    # remove fall ipes surveys (doesn't meet definition)
    !(season_f %in% c("wi") & survey_f == "ipes")
  ) %>% 
  mutate(id = row_number(),
         scale_week = ifelse(season_f == "wi", 1, -0.25)) %>%
  split(., .$id) %>%
  purrr::map(., function (x) {
    dum_grid <- if (x$season_f == "su") summer_grid else fall_grid
    
    dum_grid %>% 
      mutate(
        year = x$year,
        year_f = as.factor(x$year),
        survey_f = x$survey_f,
        season_f = x$season_f,
        scale_week = x$scale_week,
        target_depth = 0,
        day_night = "DAY",
        season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
      )
  }) %>%
  bind_rows() %>% 
  mutate(
    fake_survey = case_when(
      year > 2016 & survey_f == "hss" & season_f == "su" ~ "1",
      year < 2017 & survey_f == "ipes" & season_f == "su" ~ "1",
      TRUE ~ "0")
  ) %>% 
  filter(season_year_f %in% unique(dat$season_year_f)) %>% 
  select(-c(depth, slope, shore_dist)) %>% 
  droplevels()
exp_grid_hss <- exp_grid %>% filter(survey_f == "hss")

# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m)
sp_scalar <- 1000^2 * 13


dd <- exp_grid_hss %>% filter(year == "2017") %>% droplevels()

dum <-  sdmTMB(
  n_juv ~ 0 + season_f,# + target_depth + day_night,
  offset = dat$effort,
  data = dat,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatial_varying = ~ 0 + season_f,
  # time = "year",
  # spatiotemporal = "ar1",
  anisotropy = FALSE,  
  share_range = TRUE, 
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)

sp_preds1 <- predict(dum, 
                     newdata = dd,
                     se_fit = FALSE,
                     re_form = NULL)




# make spatial 
spatial_preds <- purrr::map(
  list(test1, test4), function (x) {
    predict(x, newdata = exp_grid, se_fit = FALSE, re_form = NULL)
  }
)
dat_tbl$sp_preds <- spatial_preds






## check FE predictions --------------------------------------------------------

nd <- expand.grid(
  year = unique(dat$year),
  season_f = c("su"),
  scale_week = seq(
    min(dat$scale_week),
    max(dat$scale_week),
    length = 30
  )
) %>% 
  mutate(
    year_f = as.factor(year),
    target_depth = 0,
    day_night = "DAY"
  )

p <- predict(dat_tbl$fits[[3]], newdata = nd, se_fit = FALSE, re_form = NA)

ggplot(p) +
  geom_line(aes(x = scale_week, y = exp(est), colour = year_f))



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

