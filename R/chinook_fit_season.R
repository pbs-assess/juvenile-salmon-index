### Juvenile all species fit - season specific 
## Use all_species_fit.R and survey_model_structure.R as templates to fit 
## species- and survey-specific models
## 1) Fit spatial model for each species to ensure reasonable convergence and
## to test for effects of survey domain on fixed effect estimates 
## 2) Fit saturated spatiotemporal model to each species
## 3) Calculate index for summer and fall (assuming surface and day tow)
## 4) Calculate fixed effects for spatial covariates
## 5) Calculate spatiotemporal effects (maps by year)
## Mar 28, 2023


library(tidyverse)
library(sdmTMB)
library(ggplot2)


# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    month_f = as.factor(month)
  ) %>% 
  filter(!species == "steelhead") %>% 
  droplevels()


## plotting color palette
col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


# separate seasons and join
summer_years <- dat %>% 
  filter(season_f == "su"#, 
         # remove 2021 since only partial survey
         # !year_f == "2021"
         ) %>%
  pull(year) %>% 
  unique()
fall_years <- dat %>% 
  filter(season_f == "wi"#,
         # remove 2020 since most of survey was outside of grid
         # !year_f == "2020"
         ) %>%
  pull(year) %>% 
  unique()


dum <- dat %>% filter(season_f == "wi", species == "chinook")
dat_coords <- dum %>% 
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
  dum,
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 

test1 <-  sdmTMB(
  n_juv ~ #month_f + 
    target_depth + day_night ,
  offset = dum$effort,
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
    # newton_loops = 1#,
    nlminb_loops = 2
  ),
  silent = FALSE
)
test2 <- update(test1, spatiotemporal = "rw")
# test3 <- update(test2, time_varying = ~1)
# test7 <-  sdmTMB(
#   n_juv ~ 0 + month_f + target_depth + day_night + survey_f + year_f,
#   offset = dum$effort,
#   data = dum,
#   mesh = spde,
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatiotemporal = "AR1",
#   time = "year",
#   anisotropy = FALSE,  
#   share_range = TRUE,
#   # extra_time = extra_time,
#   control = sdmTMBcontrol(
#     newton_loops = 1#,
#     # nlminb_loops = 2
#   ),
#   silent = FALSE
# )
test5 <- sdmTMB(
  n_juv ~ s(week, k = 4) + target_depth + day_night,
  offset = dum$effort,
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
  silent = FALSE
)
test6 <- update(test5, spatiotemporal = "rw")
test8 <- update(test3, anisotropy = FALSE)
test9 <-  sdmTMB(
  n_juv ~ 0 + target_depth + day_night + survey_f + year_f,
  offset = dum$effort,
  data = dum,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "off",
  spatiotemporal = "AR1",
  time = "year",
  anisotropy = FALSE,  
  share_range = TRUE,
  # extra_time = extra_time,
  control = sdmTMBcontrol(
    # newton_loops = 1#,
    nlminb_loops = 2
  ),
  silent = FALSE
)
test10 <-  sdmTMB(
  n_juv ~ target_depth + day_night + survey_f + month_f,
  offset = dum$effort,
  data = dum,
  mesh = spde,
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "AR1",
  time = "year",
  anisotropy = FALSE,  
  share_range = TRUE,
  # extra_time = extra_time,
  control = sdmTMBcontrol(
    newton_loops = 1#,
    # nlminb_loops = 2
  ),
  silent = FALSE
)
test11 <- update(test10, spatiotemporal = "rw")
test12 <- update(test9, spatial = "off")
test13 <- sdmTMB(
  n_juv ~ 0 + year_f + target_depth + day_night + survey_f + s(week, k = 4),
  offset = dum$effort,
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
  silent = FALSE
)

# top models
AIC(test5, test6, test13)
sanity(test12)


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

