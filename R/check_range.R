# check whether share_range can be changed to false

library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)

ncores <- parallel::detectCores()
future::plan(future::multisession, workers = 5L)


# downscale data and predictive grid
dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds"))


# make year key representing pink/sockeye cycle lines
yr_key <- data.frame(
  year = unique(dat_in$year)
) %>%
  arrange(year) %>%
  mutate(
    sox_cycle = rep(1:4, length.out = 25L) %>%
      as.factor(),
    pink_cycle = rep(1:2, length.out = 25L) %>%
      as.factor()
  )

# downscale data and predictive grid
dat <- dat_in %>%
  mutate(
    year_f = as.factor(year),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_km3 * 500),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    day_night = as.factor(day_night)) %>%
  filter(!species == "steelhead") %>%
  droplevels() %>%
  left_join(., yr_key, by = "year")


## mesh shared among species
dat_coords <- dat %>%
  filter(species == "chinook") %>%
  select(utm_x_1000, utm_y_1000) %>%
  as.matrix()
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(2, 10) * 500,
  cutoff = 30,
  offset = c(10, 50)
)
spde <- make_mesh(
  dat %>%
    filter(species == "chinook"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
)


dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest()


index_grid_hss <- readRDS(here::here("data", "index_hss_grid.rds")) %>%
  mutate(day_night = as.factor(day_night),
         trim = ifelse(
           season_f == "wi" & utm_y_1000 < 5551, "yes", "no"
         )) %>%
  #subset to northern domain
  filter(trim == "no") %>%
  left_join(., yr_key, by = "year")

sp_scalar <- 1 * (13 / 1000) * 500

# fitted model from all_species_fit.R
all_fit_tbl <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_final.rds")
)


## FIT MODEL TO EACH SPECIES ---------------------------------------------------

dat_tbl$fit <- furrr::future_map2(
  dat_tbl$data, dat_tbl$species,
  function(dat_in, sp) {
    if (sp == "pink") {
      sdmTMB(
        n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
          scale_depth + pink_cycle,
        offset = dat_in$effort,
        data = dat_in,
        mesh = spde,
        family = sdmTMB::nbinom2(),
        spatial = "off",
        spatial_varying = ~ 0 + season_f,
        time = "year",
        spatiotemporal = "rw",
        anisotropy = TRUE,
        groups = "season_f",
        share_range = FALSE,
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE
      )
    } else if (sp == "sockeye") {
      sdmTMB(
        n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
          scale_depth + sox_cycle,
        offset = dat_in$effort,
        data = dat_in,
        mesh = spde,
        family = sdmTMB::nbinom2(),
        spatial = "off",
        spatial_varying = ~ 0 + season_f,
        time = "year",
        spatiotemporal = "rw",
        anisotropy = TRUE,
        groups = "season_f",
        share_range = FALSE,
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE
      )
    } else {
      sdmTMB(
        n_juv ~ 0 + season_f + day_night + survey_f + scale_dist +
          scale_depth,
        offset = dat_in$effort,
        data = dat_in,
        mesh = spde,
        family = sdmTMB::nbinom2(),
        spatial = "off",
        spatial_varying = ~ 0 + season_f,
        time = "year",
        spatiotemporal = "rw",
        anisotropy = TRUE,
        groups = "season_f",
        share_range = FALSE,
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE
      )
    }
  }
)
