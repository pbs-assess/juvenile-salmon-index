### Survey Effects CV
## Use out of sample predictive performance to evaluate changes in forecast
## between models that include/exclude survey effects
## Jan 25, 2024


# ensure mvrfrw branch installed
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "mvrfrw")


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


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
    day_night = as.factor(day_night)
    ) %>% 
  filter(!species == "steelhead") %>% 
  droplevels() %>% 
  left_join(., yr_key, by = "year")


# define training and testing data; use 30% of post-2016 (IPES period) data as
# testing
# test_dum <- dat %>% 
#   filter(year > 2016) %>% 
#   pull(unique_event) %>% 
#   unique()
# set.seed(1234)
# test_sets <- sample(
#   test_dum, size = round(0.4 * length(test_dum), digits = 0), replace = FALSE
# )
# set.seed(1234)
# sets <- unique(dat$unique_event)
# test_sets <- sample(
#   sets, size = round(0.1 * length(sets), digits = 0), replace = FALSE
# )
# 
# dat$test <- ifelse(
#   dat$unique_event %in% test_sets,
#   TRUE,
#   FALSE
# )
# 
# train_dat <- dat %>% filter(test == FALSE)
# test_dat <- dat %>% filter(test == TRUE)

fold_ids <- dat %>% 
  select(unique_event) %>% 
  distinct() %>% 
  mutate(
    fold_id = sample(seq(1, 5, by = 1), nrow(.), replace = T)
  )



## plotting color palette
col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = 5)


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

spde$mesh$n


# shape file for coastline
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -129, ymin = 48.25, xmax = -123.25, ymax = 51.3) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))


## FIT -------------------------------------------------------------------------

train_dat_tbl <- dat %>%
  left_join(., fold_ids, by = "unique_event") %>% 
  group_by(species) %>%
  group_nest()

set.seed(123)
train_fits_list <- purrr::map2(
  train_dat_tbl$data, train_dat_tbl$species,
  function(dat_in, sp) {
    if (sp == "pink") {
      sdmTMB_cv(
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
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE,
        fold_ids = dat_in$fold_id,
        use_initial_fit = TRUE
      )
    } else if (sp == "sockeye") {
      sdmTMB_cv(
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
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE,
        fold_ids = dat_in$fold_id,
        use_initial_fit = TRUE
      )
    } else {
      sdmTMB_cv(
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
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE,
        fold_ids = dat_in$fold_id,
        use_initial_fit = TRUE
      )
    }
  }
)


# as above but no survey effects
set.seed(123)
train_fits_list2 <- purrr::map2(
  train_dat_tbl$data, train_dat_tbl$species,
  function(dat_in, sp) {
    if (sp == "pink") {
      sdmTMB_cv(
        n_juv ~ 0 + season_f + scale_dist + scale_depth + pink_cycle,
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
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE,
        fold_ids = dat_in$fold_id,
        use_initial_fit = TRUE
      )
    } else if (sp == "sockeye") {
      sdmTMB_cv(
        n_juv ~ 0 + season_f + scale_dist + scale_depth + sox_cycle,
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
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE,
        fold_ids = dat_in$fold_id,
        use_initial_fit = TRUE
      )
    } else {
      sdmTMB_cv(
        n_juv ~ 0 + season_f + scale_dist + scale_depth,
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
        control = sdmTMBcontrol(
          map = list(
            ln_tau_Z = factor(
              rep(1, times = length(unique(dat$season_f)))
            )
          )
        ),
        silent = FALSE,
        fold_ids = dat_in$fold_id, 
        use_initial_fit = TRUE
      )
    }
  }
)



set.seed(123)
coho_refit1 <- sdmTMB_cv(
  n_juv ~ 0 + season_f + scale_dist + scale_depth + day_night + survey_f,
  offset = "effort",
  data = train_dat_tbl$data[[4]],
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat$season_f)))
      )
    )
  ),
  silent = FALSE,
  fold_ids = train_dat_tbl$data[[4]]$fold_id,
  use_initial_fit = TRUE
)
set.seed(123)
coho_refit2 <- sdmTMB_cv(
  n_juv ~ 0 + season_f + scale_dist + scale_depth,
  offset =  "effort",
  data = train_dat_tbl$data[[4]],
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatial_varying = ~ 0 + season_f,
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  groups = "season_f",
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(dat$season_f)))
      )
    )
  ),
  silent = FALSE,
  fold_ids = train_dat_tbl$data[[4]]$fold_id,
  use_initial_fit = TRUE
)

train_fits_list[[3]] <- coho_refit1
train_fits_list2[[3]] <- coho_refit2

train_dat_tbl$fit_survey <- train_fits_list
train_dat_tbl$fit_no_survey <- train_fits_list2


saveRDS(
  train_dat_tbl,
  here::here("data", "fits", "cv_survey_effect_fits2.RDS")
)


## EVALUTATE WITH LOGLIK -------------------------------------------------------

train_dat_tbl$fit_survey[[1]]$sum_loglik
train_dat_tbl$fit_no_survey[[1]]$sum_loglik

purrr::map(
  train_dat_tbl$fit_no_survey, ~ .x$converged
) %>%
  unlist()

surv_loglik <- purrr::map(
  train_dat_tbl$fit_survey, ~ .x$sum_loglik
) %>% 
  unlist()
no_surv_loglik <- purrr::map(
  train_dat_tbl$fit_no_survey, ~ .x$sum_loglik
) %>% 
  unlist()

loglik_out <- data.frame(
  species = rep(train_dat_tbl$species, times = 2),
  model = rep(c("no survey effects", "survey effects"), each = 5),
  kfolds_cv_sum_loglik = c(no_surv_loglik, surv_loglik)
)


# import AIC model comparison table
aic_tab <- read.csv(
  here::here("data", "survey_eff_model_comp.csv"))

tab_out <- aic_tab %>% 
  left_join(., loglik_out, by = c("species", "model"))

write.csv(tab_out,
          here::here("data", "survey_eff_model_comp_kfolds.csv"),
          row.names = FALSE)
