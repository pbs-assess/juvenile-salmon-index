### Juvenile all species model comparison 
## Based on FE structure in all_species_fit, compete different ST model
## structures for each species using cross validation 
## Jan 21, 2022


library(tidyverse)
library(sdmTMB)
library(ggplot2)

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    week = lubridate::week(date),
    vessel = as.factor(vessel)
  ) %>% 
  # exclude 1995 and 1997 because sampling sparse
  filter(!year %in% c("1995", "1996", "1997", "2022"),
         !bath_depth_mean_m < 0) %>% 
  droplevels() 


## TOY EXAMPLE FOR BUG ---------------------------------------------------------

dd <- dat %>%
  filter(species == "CHINOOK",
         year > 2010) %>% 
  select(n_juv, volume_m3, utm_x_1000, utm_y_1000, effort)

mesh <- make_mesh(dd, c("utm_x_1000", "utm_y_1000"), cutoff = 20) # excercise

dum <- sdmTMB_cv(
  n_juv ~ 1,
  data = dd,
  mesh = mesh,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  k_folds = 5
)

dum2 <- sdmTMB_cv(
  n_juv ~ 1,
  data = dd,
  mesh = mesh,
  offset = "effort",
  family = sdmTMB::nbinom2(),
  spatial = "on",
  k_folds = 5
)


## -----------------------------------------------------------------------------

cv_tbl <- dat  %>% 
  group_by(species) %>% 
  group_nest() 

# # specify different spatiotemporal structures
# cv_tbl <- purrr::map(
#   c("off", "iid", "ar1"),
#   ~ {
#     cv_tbl %>% 
#       mutate(
#         time_varying = .x
#       )
#   }
#   ) %>% 
#   bind_rows()

## coordinates should be common to all species
dat_coords <- dat %>% 
  filter(species == "PINK") %>% 
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
  dat %>% filter(species == "PINK"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 


# spatial model  
set.seed(2022)
cv_spatial <- furrr::future_map(
  cv_tbl$data,
  ~ sdmTMB_cv(
      n_juv ~ 1 +
        as.factor(year) +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = "effort",#dat$effort,
      data = .x,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      anisotropy = TRUE,
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      use_initial_fit = TRUE,
      k_folds = 5
  ),
  .options = furrr::furrr_options(seed = TRUE)
)

# IID spatiotemporal model
set.seed(2022)
cv_st_iid <- furrr::future_map(
  cv_tbl$data,
  ~ sdmTMB_cv(
    n_juv ~ 1 +
      as.factor(year) +
      dist_to_coast_km +
      s(week, bs = "cc", k = 5) +
      target_depth +
      day_night +
      survey_f,
    offset = "effort",
    data = .x,
    mesh = spde,
    family = sdmTMB::nbinom2(),
    spatial = "on",
    spatiotemporal = "iid",
    anisotropy = TRUE,
    share_range = TRUE,
    time = "year",
    knots = list(
      week = c(0, 52)
    ),
    control = sdmTMBcontrol(
      newton_loops = 1
    ),
    use_initial_fit = TRUE,
    k_folds = 5
  ),
  .options = furrr::furrr_options(seed = TRUE)
)
    
# AR1 spatiotemporal model
set.seed(2022)
cv_st_ar1 <- furrr::future_map(
  cv_tbl$data,
  ~ sdmTMB_cv(
    n_juv ~ 1 +
      as.factor(year) +
      dist_to_coast_km +
      s(week, bs = "cc", k = 5) +
      target_depth +
      day_night +
      survey_f,
    offset = "effort",
    data = .x,
    mesh = spde,
    family = sdmTMB::nbinom2(),
    spatial = "on",
    spatiotemporal = "ar1",
    anisotropy = TRUE,
    share_range = TRUE,
    time = "year",
    knots = list(
      week = c(0, 52)
    ),
    control = sdmTMBcontrol(
      newton_loops = 1
    ),
    use_initial_fit = TRUE,
    k_folds = 5
  ),
  .options = furrr::furrr_options(seed = TRUE)
)








