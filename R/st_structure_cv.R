### Juvenile all species model comparison 
## Based on FE structure in all_species_fit, compete different ST model
## structures for each species using cross validation 
## Jan 21, 2022


library(tidyverse)
library(sdmTMB)
library(ggplot2)

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
  filter(!year %in% c("1995", "1997", "2022"),
         !bath_depth_mean_m < 0) %>% 
  droplevels() 


## TOY EXAMPLE FOR BUG ---------------------------------------------------------

dd <- dat %>%
  filter(species == "CHINOOK",
         year > 2010) %>% 
  select(n_juv, volume_m3, utm_x_1000, utm_y_1000)

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
  offset = log(dd$volume_m3),
  family = sdmTMB::nbinom2(),
  spatial = "on",
  k_folds = 5
)


## -----------------------------------------------------------------------------
cv_tbl <- dat  %>% 
  group_by(species) %>% 
  group_nest() %>% 
  mutate(
    # assign share_range based on AIC (see notes for details)
    share_range = ifelse(species %in% c("CHINOOK", "SOCKEYE"), FALSE, TRUE)
  )

# specify different spatiotemporal structures
cv_tbl <- purrr::map(
  c("off", "iid", "ar1"),
  # c(NA, "year", "year"),
  ~ {
    cv_tbl %>% 
      mutate(
        time_varying = .x
      )
  }
  ) %>% 
  bind_rows()

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


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)



# CROSS-VALIDATION -------------------------------------------------------------

cv_mods <- furrr::future_pmap(
  list(cv_tbl$data, cv_tbl$share_range, cv_tbl$time_varying),
  function(dat, sr, tv) {
    set.seed(123)
    
    if (tv == "off") {
      fit_cv <- sdmTMB_cv(
        n_juv ~ 1 +
          as.factor(year) +
          dist_to_coast_km +
          s(week, bs = "cc", k = 5) +
          target_depth +
          day_night +
          survey_f,
        offset = dat$effort,
        data = dat,
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
        k_folds = 5
      )
    } else {
      fit_cv <- sdmTMB_cv(
        n_juv ~ 1 +
          as.factor(year) +
          dist_to_coast_km +
          s(week, bs = "cc", k = 5) +
          target_depth +
          day_night +
          survey_f,
        offset = dat$effort,
        data = dat,
        mesh = spde,
        family = sdmTMB::nbinom2(),
        spatial = "on",
        spatiotemporal = tv,
        anisotropy = TRUE,
        share_range = sr,
        time = "year",
        knots = list(
          week = c(0, 52)
        ),
        extra_time = c(1997),
        control = sdmTMBcontrol(
          newton_loops = 1
        ),
        k_folds = 5
      )
    }
    return(fit_cv)
  },
  .options = furrr::furrr_options(seed = TRUE)
)

purrr::map(cv_mods, sanity)







