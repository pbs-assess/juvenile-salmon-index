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
  sample_n(size = 0.2 * nrow(.)) %>% 
  pull(unique_event)
fall_test <- dat %>% 
  filter(season_f == "wi") %>%
  sample_n(size = 0.2 * nrow(.)) %>% 
  pull(unique_event)


dat$fall_train <- ifelse(dat$unique_event %in% fall_test, 0, 1)
dat$summer_train <- ifelse(dat$unique_event %in% summer_test, 0, 1)

fall_tbl <- dat %>% 
  filter(season_f == "wi") %>% 
  mutate(dataset = "fall")
summer_tbl <- dat %>% 
  filter(season_f == "su") %>% 
  mutate(dataset = "summer")

# make tbl with two models for each season (necessary because testing/training
# data differs between full summer and fall data)
dat_tbl <- dat  %>% 
  mutate(dataset = "all_summer") %>% 
  rbind(., dat  %>% 
          mutate(dataset = "all_fall")) %>% 
  rbind(., summer_tbl) %>% 
  rbind(., fall_tbl) %>% 
  group_by(species, dataset) %>% 
  group_nest() %>% 
  mutate(
    anisotropy = ifelse(species == "pink", FALSE, TRUE)
  )
dat_tbl$fits <- vector(mode = "list", length = nrow(dat_tbl))
dat_tbl$preds <- vector(mode = "list", length = nrow(dat_tbl))

for (i in 2:4) {
  dum <- dat_tbl$data[[i]]
  
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
  
  bs <- ifelse(grepl("all", dat_tbl$dataset[i]), "cc", "tp")
  knts <- ifelse(grepl("all", dat_tbl$dataset[i]), 
                 list(
                   week = c(0, 52)
                 ),
                 NULL)
  
  dat_tbl$fits[[i]] <- sdmTMB(
    n_juv ~ 0 +
      year_f +
      dist_to_coast_km +
      s(week, bs = bs, k = 5) +
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
    share_range = FALSE,
    knots = knts,
    control = sdmTMBcontrol(
      newton_loops = 1
    ),
    silent = FALSE
  )
  
  dd_test <- dum %>% filter(remove == "1")
  dat_tbl$preds[[i]] <- predict(dat_tbl$fits[[i]], newdata = dd_test)
}

saveRDS(dat_tbl, here::here("data", "fits", "model_structure_fits.rds"))

purrr::map(dat_tbl$preds[1:4], 
           ~ Metrics::rmse(
             .x$n_juv,
             exp(.x$est)
           ))

Metrics::rmse(dat_tbl$preds[[i]]$n_juv, exp(dat_tbl$preds[[i]]$est))
