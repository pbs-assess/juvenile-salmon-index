### Chinook model structure test
## Split summer and fall each into testing/training; fit full and season-
## constrained models, then compare with RMSE
## Use Chinook initially, then fit to all species
## Mar 13, 2022


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
  filter(!species == "steelhead") %>% 
  droplevels()


# identify summer and fall testing data
set.seed(123)
summer_test <- dat %>% 
  filter(season_f == "su") %>%
  pull(unique_event) %>% 
  unique() %>% 
  sample(size = 0.2 * length(.))
fall_test <- dat %>% 
  filter(season_f == "wi")  %>%
  pull(unique_event) %>% 
  unique() %>% 
  sample(size = 0.2 * length(.))
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


dat$fall_train <- ifelse(dat$unique_event %in% fall_test, 0, 1)
dat$summer_train <- ifelse(dat$unique_event %in% summer_test, 0, 1)


# fall_tbl <- dat %>% 
#   filter(season_f == "wi",
#          year %in% fall_years) %>% 
#   mutate(dataset = "fall") %>% 
#   droplevels()
# summer_tbl <- dat %>% 
#   filter(season_f == "su",
#          year %in% summer_years) %>% 
#   mutate(dataset = "summer") %>% 
#   droplevels()

# make tbl with two models for each season (necessary because testing/training
# data differs between full summer and fall data)
# dat_tbl <- dat  %>% 
#   mutate(dataset = "all_summer") %>% 
#   rbind(., dat  %>% 
#           mutate(dataset = "all_fall")) %>% 
#   rbind(., summer_tbl) %>% 
#   rbind(., fall_tbl) %>% 
#   group_by(species, dataset) %>% 
#   group_nest() %>% 
#   mutate(
#     anisotropy = ifelse(species == "pink", FALSE, TRUE)
#   )
# dat_tbl$fits <- vector(mode = "list", length = nrow(dat_tbl))
# dat_tbl$preds <- vector(mode = "list", length = nrow(dat_tbl))
dat_tbl$cv_fits <- vector(mode = "list", length = nrow(dat_tbl))

dat_tbl <- readRDS(here::here("data", "fits", "model_structure_fits.rds"))

library(future)
plan(multisession)

for (i in 1:4) {
  dd <- dat_tbl$data[[i]] %>% droplevels()
  
  dat_coords <- dd %>% 
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
    dd,
    c("utm_x_1000", "utm_y_1000"),
    mesh = inla_mesh_raw
  ) 
  
  bs <- ifelse(grepl("all_", dat_tbl$dataset[i]), "cc", "tp")
  knts <- NULL
  if (grepl("all_", dat_tbl$dataset[i])) {
    knts <- list(
      week = c(0, 52)
    )
  }
  extra_time <- NULL
  if (dat_tbl$dataset[i] == "fall") {
    extra_time <-  c(2018, 2020) 
  }
  if (dat_tbl$dataset[i] == "summer") {
    extra_time <-  c(2016, 2020, 2021) 
  }
  
  if (dat_tbl$dataset[[i]] == "fall") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ #dist_to_coast_km +
        week +
        target_depth +
        day_night ,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = TRUE,  
      share_range = TRUE,
      knots = knts,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  } 
  if (dat_tbl$dataset[[i]] == "summer" & dat_tbl$species[[i]] == "pink") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ #dist_to_coast_km +
        #as.factor(month) +
        target_depth +
        day_night +
        survey_f,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = TRUE, 
      share_range = TRUE,
      knots = knts,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      parallel = TRUE,
      use_initial_fit = FALSE,
      k_folds = 5
    )
  }
  
  if (dat_tbl$dataset[[i]] == "summer") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ #dist_to_coast_km +
        as.factor(month) +
        target_depth +
        day_night +
        survey_f,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = TRUE, 
      share_range = TRUE,
      knots = knts,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      parallel = TRUE,
      use_initial_fit = FALSE,
      k_folds = 5
    )
  }
  if (dat_tbl$dataset[[i]] == "all_fall") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ 0 + year_f + 
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night ,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = dat_tbl$anisotropy[i], 
      share_range = TRUE,
      knots = knts,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      parallel = TRUE,
      use_initial_fit = FALSE,
      k_folds = 5
    )
  }
  if (dat_tbl$dataset[[i]] == "all_summer") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ 0 + year_f + 
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = dat_tbl$anisotropy[i], 
      share_range = TRUE,
      knots = knts,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      parallel = TRUE,
      use_initial_fit = FALSE,
      k_folds = 5
    )
  }
}


dd <- purrr::pmap(
  list(dat_tbl$fits[1:4], dat_tbl$data[1:4], dat_tbl$dataset[1:4]), 
  function (x, y, z) {
    if (grepl("summer", z)) {
      y$remove <- ifelse(y$unique_event %in% summer_test, 1, 0)
    }
    
    if (grepl("fall", z)) {
      y$remove <- ifelse(y$unique_event %in% fall_test, 1, 0)
    }
    dd_test <- y %>% filter(remove == "1")
    predict(x, newdata = dd_test, nsim = 10)
  }
)


for (i in 16:20) {
  dum <- dat_tbl$data[[i]] %>% droplevels()
  
  if (grepl("summer", dat_tbl$dataset[[i]])) {
    dum$remove <- ifelse(dum$unique_event %in% summer_test, 1, 0)
  }
  
  if (grepl("fall", dat_tbl$dataset[[i]])) {
    dum$remove <- ifelse(dum$unique_event %in% fall_test, 1, 0)
  }
  
  dd <- dum %>% filter(!remove == "1")
  
  dat_coords <- dd %>% 
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
    dd,
    c("utm_x_1000", "utm_y_1000"),
    mesh = inla_mesh_raw
  ) 
  
  bs <- ifelse(grepl("all_", dat_tbl$dataset[i]), "cc", "tp")
  kk <- ifelse(grepl("all_", dat_tbl$dataset[i]), 5, 3)
  knts <- NULL
  if (grepl("all_", dat_tbl$dataset[i])) {
    knts <- list(
      week = c(0, 52)
    )
  }
  extra_time <- NULL
  if (dat_tbl$dataset[i] == "fall") {
    extra_time <-  c(2018, 2020) 
  }
  if (dat_tbl$dataset[i] == "summer") {
    extra_time <-  c(2016, 2020, 2021) 
  }
  
  if (dat_tbl$dataset[[i]] == "fall") {
    dat_tbl$fits[[i]] <- sdmTMB(
      n_juv ~ #dist_to_coast_km +
        week +
        target_depth +
        day_night ,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = TRUE,  
      share_range = TRUE,
      knots = knts,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  } 
  if (dat_tbl$dataset[[i]] == "summer" & dat_tbl$species[[i]] == "pink") {
    dat_tbl$fits[[i]] <- sdmTMB(
      n_juv ~ #dist_to_coast_km +
        #as.factor(month) +
        target_depth +
        day_night +
        survey_f,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = TRUE, 
      share_range = TRUE,
      knots = knts,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  }
  
  if (dat_tbl$dataset[[i]] == "summer") {
    dat_tbl$fits[[i]] <- sdmTMB(
      n_juv ~ #dist_to_coast_km +
        as.factor(month) +
        target_depth +
        day_night +
        survey_f,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = TRUE, 
      share_range = TRUE,
      knots = knts,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  }
  if (dat_tbl$dataset[[i]] == "all_fall") {
    dat_tbl$fits[[i]] <- sdmTMB(
      n_juv ~ 0 + year_f + 
        dist_to_coast_km +
        s(week, bs = bs, k = kk) +
        target_depth +
        day_night ,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = dat_tbl$anisotropy[i], 
      share_range = TRUE,
      knots = knts,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  }
  if (dat_tbl$dataset[[i]] == "all_summer") {
    dat_tbl$fits[[i]] <- sdmTMB(
      n_juv ~ 0 + year_f + 
        dist_to_coast_km +
        s(week, bs = bs, k = kk) +
        target_depth +
        day_night +
        survey_f,
      offset = dd$effort,
      data = dd,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = dat_tbl$anisotropy[i], 
      share_range = TRUE,
      knots = knts,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  }

  
  dd_test <- dum %>% filter(remove == "1")
  dat_tbl$preds[[i]] <- predict(dat_tbl$fits[[i]], newdata = dd_test)
}

saveRDS(dat_tbl, here::here("data", "fits", "model_structure_fits.rds"))

purrr::map(dat_tbl$fits[1:4], 
           sanity)
dat_tbl[16, ]

# despite differences in predictions, very similar RMSE values;
purrr::map(dat_tbl$preds[5:12], 
           ~ Metrics::rmse(
             .x$n_juv,
             exp(.x$est)
           ))




