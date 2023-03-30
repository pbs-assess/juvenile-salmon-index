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
library(sdmTMBextra)


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
fall_tbl <- dat %>%
  filter(season_f == "wi",
         year %in% fall_years) %>%
  mutate(dataset = "fall") %>%
  droplevels()
summer_tbl <- dat %>%
  filter(season_f == "su",
         year %in% summer_years) %>%
  mutate(dataset = "summer") %>%
  droplevels()

dat_tbl <- rbind(summer_tbl, fall_tbl) %>%
  group_by(species, dataset) %>%
  group_nest() %>% 
  mutate(
    anisotropy = ifelse(species == "chinook" & dataset == "summer", FALSE, TRUE),
    time_model = ifelse(species == "sockeye" & dataset == "fall", "rw", "ar1")
  )


dat_tbl$fits <- purrr::pmap(
  list(dat_tbl$dataset, dat_tbl$time_model, dat_tbl$anisotropy, dat_tbl$data),
  function(x, y, z, dat_in) {
    dum <- dat_in %>% droplevels()
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
    
    extra_time <- NULL
    if (x == "fall") {
      extra_time <-  c(2018#, 2020
      )
      formula_in <- as.formula(
        paste("n_juv ~ target_depth + day_night")
      )
    }
    if (x == "summer") {
      extra_time <-  c(2016, 2020
                       #, 2021
      )
      formula_in <- as.formula(
        paste("n_juv ~ survey_f + target_depth + day_night")
      )
    }
    sdmTMB(
      formula_in,
      offset = dum$effort,
      data = dum,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = y,
      time = "year",
      anisotropy = z,  
      share_range = TRUE,
      extra_time = extra_time,
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  } 
)

dat_tbl$fits <- fits_list

# saveRDS(dat_tbl, here::here("data", "fits", "st_mod_all_sp_season.rds"))

# sub_tbl <- dat_tbl %>%
#   filter(
#     (species == "sockeye" & dataset == "fall" & model == "rw") |
#       (!(species == "sockeye" & dataset == "fall") & model == "ar1")
#   )
# # saveRDS(sub_tbl, here::here("data", "fits", "top_st_mod_all_sp_season.rds"))


sub_tbl %>% 
  filter(dataset == "summer") %>% 
  pull(fits) %>% 
  purrr::map(., sanity)


# check residuals
sub_tbl$sims <- purrr::map(sub_tbl$fits, simulate, nsim = 50)
qq_list <- purrr::map2(sub_tbl$sims, sub_tbl$fits, dharma_residuals)


## INDEX -----------------------------------------------------------------------


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


# summer high seas
ind_preds_sum <- purrr::map(
  sub_tbl$fits,
  ~ {
    predict(.x,
            newdata = exp_grid %>% filter(survey_f == "hss", week == "25"),
            return_tmb_object = TRUE)
  }
)
index_list_sum <- furrr::future_map(
  ind_preds_sum, 
  get_index, 
  area = sp_scalar, 
  bias_correct = TRUE
)

# fall high seas
ind_preds_fall <- purrr::map(
  sub_tbl$fits,
  ~ {
    predict(.x,
            newdata = exp_grid %>% filter(survey_f == "hss", week == "42"),
            return_tmb_object = TRUE)
  }
)
index_list_fall <- furrr::future_map(
  ind_preds_fall, get_index, 
  area = sp_scalar, 
  bias_correct = TRUE
)

index_sum <- purrr::map2(
  index_lists[1:5], dat_tbl$species, function (x, sp) {
    x$species <- sp
    x$season <- "su"
    return(x)
  }) %>% 
  bind_rows()
index_fall <- purrr::map2(
  index_lists[6:10], dat_tbl$species, function (x, sp) {
    x$species <- sp
    x$season <- "fa"
    return(x)
  }) %>% 
  bind_rows()


index_dat <- rbind(index_sum %>% filter(year %in% summer_years), 
                   index_fall %>% filter(year %in% fall_years)) %>% 
  mutate(
    season = factor(season, levels = c("su", "fa"), 
                    labels = c("summer", "fall"))
  ) %>% 
  group_by(species) %>% 
  mutate(
    max_est = max(est),
    scale_est = est / max_est,
    scale_lwr = lwr / max_est,
    scale_upr = upr / max_est
  ) %>% 
  ungroup()

## scaled abundance with no intervals
index <- ggplot(index_dat, aes(year, est)) +
  geom_point(aes(fill = species), shape = 21) +
  # geom_pointrange(aes(ymin = lwr, ymax = upr, fill = species),
  #                 shape = 21) +
  labs(x = "Year", y = "Absolute Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none")


# SPATIAL PREDS ----------------------------------------------------------------

##NOTE HASN'T BEEN TWEAKED TO MATCH NEW RUNS


# make spatial 
spatial_preds <- furrr::future_map(
  dat_tbl$st_mod, function (x) {
    predict(x, newdata = exp_grid, se_fit = FALSE, re_form = NULL)
  },
  .options = furrr::furrr_options(seed = TRUE)
)


# make new tibble of predictions
spatial_pred_tbl <- tibble(
  species = rep(tolower(unique(dat$species)), each = 2),
  season = rep(c("summer", "fall"), times = 5),
  week = rep(unique(exp_grid$week), times = 5)
)
# separate seasons into their own list
spatial_pred_tbl$spatial_preds <- purrr::map(
  spatial_preds,
  ~ {
    split(.x, .x$week)
  }
) %>% 
  do.call(c, .)


# which years are viable?
# index from summer
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


# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells
sp_scalar <- 1000^2 * 13


# spatiotemporal and total predsrandom fields for subset of years
year_seq <- seq(1999, 2019, by = 5)
sub_spatial <- spatial_pred_tbl %>% 
  filter(season == "fall") %>% 
  select(-week) %>% 
  unnest(cols = "spatial_preds") %>%
  group_by(species) %>% 
  mutate(
    grid_est = sp_scalar * exp(est),
    scale_est = grid_est / max(grid_est),
    # fix large scaled values
    scale_est2 = ifelse(scale_est > 0.1, 
                        0.1, 
                        scale_est)
  ) %>% 
  ungroup() %>% 
  filter(year %in% year_seq)
eps_max <- quantile(sub_spatial$epsilon_st, 0.999)#max(abs(epsilon_dat$epsilon_st))
max_est <- quantile(sub_spatial$scale_est, 0.99)

png(here::here("figs", "ms_figs", 
               "scaled_sp_preds.png"), 
    height = 8, width = 8, units = "in", res = 200)
ggplot() + 
  geom_raster(data = sub_spatial, aes(X, Y, fill = scale_est2)) +
  coord_fixed() +
  # geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Scaled\nAbundance",
    labels = c("0", "0.025", "0.05", "0.075", ">0.1")
  ) +
  facet_grid(year~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
dev.off()


png(here::here("figs", "ms_figs", "st_rf.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = sub_spatial, aes(X, Y, fill = epsilon_st)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(-1 * eps_max, eps_max),
                       name = "Spatial\nField") +
  facet_grid(year~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
dev.off()


# spatial random effects by species
omega_dat <- spatial_pred_tbl %>% 
  select(-week) %>% 
  unnest(cols = "spatial_preds") %>% 
  select(-c(year, fake_survey, est, est_non_rf, est_rf, epsilon_st)) %>%
  # remove duplicated summer data
  filter(survey_f == "hss",
         season == "fall"
  ) %>% 
  distinct()
omega_max <- max(abs(omega_dat$omega_s))

png(here::here("figs", "ms_figs", "spatial_rf.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = omega_dat, aes(X, Y, fill = omega_s)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(-1 * omega_max, omega_max),
                       name = "Spatial\nField") +
  facet_wrap(~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
dev.off()