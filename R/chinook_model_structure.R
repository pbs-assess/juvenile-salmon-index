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
      knots = list(
        week = c(0, 52)
      ),
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
      knots = list(
        week = c(0, 52)
      ),
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
purrr::map(dat_tbl$preds[1:4], 
           ~ Metrics::rmse(
             .x$n_juv,
             exp(.x$est)
           ))


## simulate data and check RMSE
set.seed(345)
sims_list <- purrr::map(dat_tbl$fits, simulate, nsim = 25)
rmse_list <- purrr::map2(sims_list, dat_tbl$fits, function (z, y) {
  apply(z, 2, function (x) {
    Metrics::rmse(y$data$n_juv, x)
  }
  )
})
dat_tbl$train_rmse <- rmse_list

train_rmse_plot <- dat_tbl %>% 
  select(species, dataset, train_rmse) %>% 
  unnest(cols = c(train_rmse)) %>%
  ggplot() +
  geom_boxplot(aes(x = dataset, y = train_rmse)) +
  facet_wrap(~species, scales = "free_y")

png(here::here("figs", "diagnostics", "in_sample_rmse.png"), units = "in", 
    height = 6, width = 6, res = 150)
train_rmse_plot
dev.off()

  
## cross validation ------------------------------------------------------------

library(future)
plan(multisession)

for (i in 2:4) {
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
  
  extra_time <- NULL
  if (dat_tbl$dataset[i] == "fall") {
    extra_time <-  c(2018, 2020) 
  }
  if (dat_tbl$dataset[i] == "summer") {
    extra_time <-  c(2016, 2020, 2021) 
  }
  
  if (dat_tbl$dataset[[i]] == "fall") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ week +
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
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      k_folds = 4,
      use_initial_fit = TRUE
    )
  } 
  if (dat_tbl$dataset[[i]] == "summer" & dat_tbl$species[[i]] == "pink") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ target_depth +
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
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      k_folds = 4,
      use_initial_fit = TRUE
    )
  }
  if (dat_tbl$dataset[[i]] == "summer") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ as.factor(month) +
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
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      k_folds = 4,
      use_initial_fit = TRUE
    )
  }
  if (dat_tbl$dataset[[i]] == "all_fall") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ s(week, bs = "cc", k = 5) +
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
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      k_folds = 4,
      use_initial_fit = TRUE
    )
  }
  if (dat_tbl$dataset[[i]] == "all_summer") {
    dat_tbl$cv_fits[[i]] <- sdmTMB_cv(
      n_juv ~ s(week, bs = "cc", k = 5) +
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
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE,
      k_folds = 4,
      use_initial_fit = TRUE
    )
  }
}

saveRDS(dat_tbl, here::here("data", "fits", "model_structure_fits_cv.rds"))



## spatial predictions and index -----------------------------------------------

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
               utm_y_1000 = Y / 1000,
               dist_to_coast_km = shore_dist / 1000)}
  )

summer_grid <- grid_list$ipes_grid
fall_grid <- grid_list$wcvi_grid

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(dat$year),
  survey_f = unique(dat$survey_f),
  week = c(25, 42)
) %>%
  filter(
    # remove fall ipes surveys (doesn't meet definition)
    !(week == "42" & survey_f == "ipes")
  ) %>% 
  mutate(id = row_number()) %>%
  split(., .$id) %>%
  purrr::map(., function (x) {
    dum_grid <- if (x$week == "42") fall_grid else summer_grid
    
    dum_grid %>% 
      mutate(
        year = x$year,
        year_f = as.factor(x$year),
        survey_f = x$survey_f,
        target_depth = 0,
        week = x$week,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  mutate(
    fake_survey = case_when(
      year > 2016 & survey_f == "hss" & week == "25" ~ "1",
      year < 2017 & survey_f == "ipes" & week == "25" ~ "1",
      TRUE ~ "0")
  )
sp_scalar <- 1000^2 * 13


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


spatial_preds <- furrr::future_map2(
  dat_tbl$fits[1:4], dat_tbl$dataset[1:4], function (x, y) {
    gg <- if (grepl("fall", y)) {
      exp_grid %>% filter(week == "42")
    } else {
        exp_grid %>% filter(week == "25")
      }
    predict(x, newdata = gg, se_fit = FALSE, re_form = NULL)
  },
  .options = furrr::furrr_options(seed = TRUE)
)

# summer high seas
ind_preds <- furrr::future_map2(
  dat_tbl$fits[1:4], dat_tbl$dataset[1:4], function (x, y) {
    gg <- if (grepl("fall", y)) {
      exp_grid %>% filter(week == "42")
    } else {
      exp_grid %>% filter(week == "25")
    }
    predict(x,
            newdata = gg %>% filter(survey_f == "hss"),
            return_tmb_object = TRUE)
  }
)
index_list <- furrr::future_map(
  ind_preds, 
  get_index, 
  area = sp_scalar, 
  bias_correct = TRUE
)

index_fall <- rbind(
  index_list[[1]] %>% 
    mutate(dataset = "all"),
  index_list[[3]] %>% 
    mutate(dataset = "subset")) %>% 
  mutate(
    season = "fall"
  )
index_summer <- rbind(
  index_list[[2]] %>% 
    mutate(dataset = "all"),
  index_list[[4]] %>% 
    mutate(dataset = "subset")) %>% 
  mutate(
    season = "summer"
  )

index_all <- rbind(index_fall, index_summer)

ggplot(index_all, aes(year, log_est)) +
  geom_point(aes(fill = dataset), shape = 21) +
  # geom_pointrange(aes(ymin = lwr, ymax = upr, fill = species),
  #                 shape = 21) +
  labs(x = "Year", y = "Absolute Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(dataset~season) +
  theme(legend.position = "none")
