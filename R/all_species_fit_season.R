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
    adj_year = ifelse(season_f == "sp", year - 1, year),
    year_f = as.factor(adj_year),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    month_f = as.factor(month),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    scale_season = scale(as.numeric(season_f))[ , 1],
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night)
  ) %>% 
  filter(!species == "steelhead") %>% 
  droplevels()


## plotting color palette
col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')



# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


## mesh shared among species
dat_coords <- dat %>% 
  filter(species == "chinook") %>% 
  select(utm_x_1000, utm_y_1000) %>% 
  as.matrix()
inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = dat_coords,
  max.edge = c(1, 5) * 500,
  cutoff = 20,
  offset = c(20, 200)
) 
spde <- make_mesh(
  dat %>% 
    filter(species == "chinook"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 

dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest() #%>%
  # mutate(
    # anisotropy = ifelse(species == "chinook" & dataset == "summer", FALSE, TRUE),
    # time_model = ifelse(species == "sockeye" & dataset == "fall", "rw", "ar1")
  # )


fits_list <- furrr::future_map(
  dat_tbl$data,
  function(dat_in) {
    # fit_ri <- sdmTMB(
    #   n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f +
    #      dist_to_coast_km + (1 | season_year_f),
    #   offset = dat_in$effort,
    #   data = dat_in,
    #   mesh = spde,
    #   family = sdmTMB::nbinom2(),
    #   spatial = "off",
    #   spatial_varying = ~ 0 + season_f,
    #   # time = "year",
    #   # spatiotemporal = "ar1",
    #   anisotropy = FALSE,
    #   share_range = TRUE,
    #   priors = sdmTMBpriors(
    #     phi = halfnormal(0, 10),
    #     matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
    #     matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
    #   ),
    #   control = sdmTMBcontrol(
    #     newton_loops = 1
    #   ),
    #   silent = FALSE
    # )
    # fit_tv <- 
      sdmTMB(
        n_juv ~ 0 + year_f + season_f + target_depth + day_night + survey_f + 
          dist_to_coast_km,
        offset = dat_in$effort,
        data = dat_in,
        mesh = spde,
        family = sdmTMB::nbinom2(),
        spatial = "off",
        spatial_varying = ~ 0 + season_f,
        time_varying = ~ 1 + season_f,
        time_varying_type = "rw0",
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
          newton_loops = 1,
          start = list(
            ln_tau_V = matrix(log(c(0.25, 0.25, 0.25)), nrow = 3, ncol = 1)
          ),
          map = list(
            ln_tau_V = factor(c(NA, NA, NA))
          )
        ),
        silent = FALSE
      )
    # list(fit_ri, fit_tv)
  }
)



dat_tbl$model <- "tv"
dat_tbl$fit <- fits_list

saveRDS(dat_tbl, here::here("data", "fits", "all_spatial_varying_eps.rds"))
dat_tbl <- readRDS(here::here("data", "fits", "all_spatial_varying_no_eps.rds"))


purrr::map(dat_tbl$fit, sanity)
dat_tbl$AIC <- purrr::map(dat_tbl$fit, AIC) %>% unlist()


# check residuals
dat_tbl$sims <- purrr::map(dat_tbl$fit, simulate, nsim = 50)
qq_list <- purrr::map2(dat_tbl$sims, dat_tbl$fit, dharma_residuals)


## PARAMETER ESTIMATES ---------------------------------------------------------

# fixed parameter estimates
fix_pars <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species, dat_tbl$model),
  function(x, y, z) {
    tidy(x, effects = "fixed", conf.int = T) %>% 
      mutate(
        species = y,
        model = z
      ) 
  }
) %>% 
  bind_rows() %>% 
  filter(term %in% c("season_fsu", "season_fwi", "survey_fipes", 
                     "day_nightNIGHT")) %>% 
  mutate(
    term = fct_recode(as.factor(term), 
                      "Summer" = "season_fsu", "Fall" = "season_fwi",
                      "IPES Survey" = "survey_fipes", 
                      "Nocturnal Sampling" = "day_nightNIGHT")
  )

png(here::here("figs", "ms_figs_season", "fix_ints.png"), height = 3, width = 4,
    units = "in", res = 200)
ggplot(
  fix_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high, fill = species,
      shape = model)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(23, 24)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  facet_wrap(~term, nrow = 2) +
  labs(y = "Parameter Estimate") +
  guides(fill = "none") +
  theme(axis.title.x = element_blank())
dev.off()


## random parameter estimates
ran_pars <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species, dat_tbl$model),
  function(x, y, z) {
    tidy(x, effects = "ran_par", conf.int = T) %>% 
      mutate(
        species = y,
        model = z
      ) 
  }
) %>% 
  bind_rows() 

png(here::here("figs", "ms_figs_season", "ran_pars.png"), height = 4, width = 8,
    units = "in", res = 200)
ggplot(
  ran_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high,
      fill = species, shape = model)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(23, 24)) +
  ylab("Parameter Estimate") +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  facet_wrap(~term, scales = "free_y") +
  guides(fill = "none")
dev.off()


# headrope depth
target_dat <- expand.grid(
  day_night = unique(dat$day_night),
  season_f = unique(dat$season_f),
  target_depth = seq(min(dat$target_depth), 
                     max(dat$target_depth), length.out = 100),
  dist_to_coast_km = median(dat$dist_to_coast_km),
  survey_f = unique(dat$survey_f),
  year = unique(dat$year)
) %>% 
  mutate(
    year_f = as.factor(year),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    day_night == "DAY", season_f == "su", survey_f == "hss", year_f == "2012"
  )

depth_preds <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species, dat_tbl$model),
  function(x, y, z) {
    predict(x, newdata = target_dat, se_fit = T, re_form = NA, 
            re_form_iid = NA) %>% 
      mutate(species = y,
             model = z)
  }
) %>%
  bind_rows()

depth_preds2 <- depth_preds %>% 
  group_by(species, model) %>% 
  mutate(
    exp_est = exp(est),
    max_est = max(exp_est),
    scale_est = exp_est / max_est,
    up = (est + 1.96 * est_se),
    lo = (est - 1.96 * est_se),
    exp_up = exp(up),
    exp_lo = exp(lo),
    scale_up = exp(up) / max_est,
    scale_lo = exp(lo) / max_est
  )

png(here::here("figs", "ms_figs_season", "depth_preds.png"), 
    height = 4, width = 8, units = "in", res = 200)
ggplot(
  data = depth_preds2,
  aes(x = target_depth, y = scale_est, ymin = scale_lo, ymax = scale_up, 
      fill = species)
) +
  ggsidekick::theme_sleek() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(model~species) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Target Headrope Depth (m)") +
  guides(fill = "none") +
  ylab("Abundance Index")
dev.off()



# distance to coast
dist_dat <- expand.grid(
  dist_to_coast_km = seq(min(dat$dist_to_coast_km), 
                         max(dat$dist_to_coast_km), length.out = 100),
  day_night = unique(dat$day_night),
  season_f = unique(dat$season_f),
  target_depth = 0,
  survey_f = unique(dat$survey_f),
  year = unique(dat$year)
  ) %>% 
  mutate(
    year_f = as.factor(year),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    day_night == "DAY", season_f == "su", survey_f == "hss", year_f == "2012"
  )

dist_preds <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species, dat_tbl$model),
  function(x, y, z) {
    predict(x, newdata = dist_dat, se_fit = T, re_form = NA, 
            re_form_iid = NA) %>% 
      mutate(species = y,
             model = z)
  }
) %>%
  bind_rows()

dist_preds2 <- dist_preds %>% 
  group_by(species, model) %>% 
  mutate(
    exp_est = exp(est),
    max_est = max(exp_est),
    scale_est = exp_est / max_est,
    up = (est + 1.96 * est_se),
    lo = (est - 1.96 * est_se),
    exp_up = exp(up),
    exp_lo = exp(lo),
    scale_up = exp(up) / max_est,
    scale_lo = exp(lo) / max_est
  )

ggplot(
  data = dist_preds2,
  aes(x = dist_to_coast_km, y = scale_est, ymin = scale_lo, ymax = scale_up, 
      fill = species)
) +
  ggsidekick::theme_sleek() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(model~species) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Target Headrope Depth (m)") +
  guides(fill = "none") +
  ylab("Abundance Index")



## SPATIAL GRID  ---------------------------------------------------------------

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

pred_grid_list <- expand.grid(
  year = unique(dat$year),
  survey_f = unique(dat$survey_f),
  season_f = unique(dat$season_f)
) %>%
  mutate(
    year_f = as.factor(year),
    id = row_number(),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    # remove fall ipes surveys (doesn't meet definition)
    !(season_f %in% c("wi") & survey_f == "ipes"),
    !season_f == "sp"#,
    # season_year_f %in% unique(dat$season_year_f)
  ) %>% 
  split(., .$id)


## INDEX -----------------------------------------------------------------------

# add unique years and seasons
index_grid <- pred_grid_list %>%
  purrr::map(., function (x) {
    dum_grid <- if (x$season_f == "su") summer_grid else fall_grid
    
    dum_grid %>% 
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        season_f = x$season_f,
        season_year_f = x$season_year_f,
        target_depth = 0,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  mutate(
    fake_survey = case_when(
      year > 2016 & survey_f == "hss" & season_f == "su" ~ "1",
      year < 2017 & survey_f == "ipes" & season_f == "su" ~ "1",
      TRUE ~ "0")
  ) %>% 
  select(-c(depth, slope, shore_dist)) 
# remove season-year factor levels not present in original dataframe
levels(index_grid$season_year_f)[!levels(index_grid$season_year_f) %in%
                                 levels(dat$season_year_f)] <- NA


## crashes with RI model so use TV only
sub_tbl <- dat_tbl# %>% 
  # filter(model == "tv")


# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m)
sp_scalar <- 1000^2 * 13


# summer high seas
ind_preds_sum <- purrr::map(
  sub_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid %>% 
              filter(survey_f == "hss", season_f == "su"),
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
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
  sub_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid %>% 
              filter(survey_f == "hss", season_f == "wi"),
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)
index_list_fall <- furrr::future_map(
  ind_preds_fall, 
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)

index_sum <- purrr::map2(
  index_list_sum,
  sub_tbl$species,
  function (x, sp) {
    x$species <- sp
    x$season <- "su"
    return(x)
  }) %>%
  bind_rows()
index_fall <- purrr::map2(
  index_list_fall,
  sub_tbl$species,
  function (x, sp) {
    x$species <- sp
    x$season <- "fa"
    return(x)
  }) %>%
  bind_rows()


# identify viable years by season
summer_years <- dat %>%
  filter(season_f == "su",
         # remove 2021 since few tows
         !year == "2021") %>%
  pull(year) %>%
  unique() %>% 
  sort()
fall_years <- dat %>%
  filter(season_f == "wi",
         # remove 2020 since majority of tows in north
         !year == "2020") %>%
  pull(year) %>%
  unique() %>% 
  sort()

index_dat <- rbind(index_sum, #%>% filter(year %in% summer_years),
                   index_fall #%>% filter(year %in% fall_years)
                   ) %>%
  mutate(
    season = factor(season, levels = c("su", "fa"),
                    labels = c("summer", "fall")),
    log_lwr = log_est - (1.96 * se),
    log_upr = log_est + (1.96 * se),
    survey = case_when(
      season == "fall" & year %in% fall_years ~ "sampled",
      season == "summer" & year %in% summer_years ~ "sampled",
      TRUE ~ "no survey"
    )
  ) %>%
  group_by(species, season) %>%
  mutate(
    mean_log_est = mean(log_est)
    # max_est = max(est),
    # scale_est = est / max_est,
    # scale_lwr = lwr / max_est,
    # scale_upr = upr / max_est
  ) %>%
  ungroup()

saveRDS(index_dat, here::here("data", "fits", "season_index.rds"))

index_dat <- readRDS(here::here("data", "fits", "season_index.rds"))


## scaled abundance with no intervals
shape_pal <- c(21, 23)
names(shape_pal) <- unique(index_dat$survey)
log_index <- ggplot(index_dat, aes(year, log_est)) +
  geom_pointrange(aes(ymin = log_lwr, ymax = log_upr, fill = species, 
                      shape = survey)) +
  geom_hline(aes(yintercept = mean_log_est), lty = 2) +
  labs(x = "Year", y = "Log Abundance Index") +
  scale_shape_manual(values = shape_pal) +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none")



scaled_index <- ggplot(index_dat, aes(year, scale_est)) +
  geom_pointrange(aes(ymin = scale_lwr, ymax = scale_upr, fill = species),
                  shape = 21) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(x = "Year", y = "Log Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none")



# SPATIAL PREDS ----------------------------------------------------------------

# similar to index grid except only HSS and fall survey domain to highlight 
# spatial contrasts
spatial_grid <- pred_grid_list %>%
  purrr::map(., function (x) {
    fall_grid %>% 
      mutate(
        year = x$year,
        year_f = x$year_f,
        survey_f = x$survey_f,
        season_f = x$season_f,
        season_year_f = x$season_year_f,
        target_depth = 0,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  filter(
    !survey_f == "ipes"
  ) %>% 
  select(-c(depth, slope, shore_dist)) 
# remove season-year factor levels not present in original dataframe
levels(index_grid$season_year_f)[!levels(index_grid$season_year_f) %in%
                                   levels(dat$season_year_f)] <- NA


# make spatial 
spatial_preds <- furrr::future_map2(
  dat_tbl$fit, dat_tbl$species, function (x, y) {
    predict(x, newdata = spatial_grid, 
            se_fit = FALSE, re_form = NULL) %>% 
      mutate(species = y)
  },
  .options = furrr::furrr_options(seed = TRUE)
)
dat_tbl$sp_preds <- spatial_preds


# spatiotemporal and total preds random fields for subset of years
year_seq <- seq(1999, 2019, by = 5)
sub_spatial <- dat_tbl %>% 
  unnest(cols = "sp_preds") %>%
  group_by(species, dataset) %>% 
  mutate(
    grid_est = sp_scalar * exp(est),
    scale_est = grid_est / max(grid_est),
    # fix large scaled values
    scale_est2 = ifelse(scale_est > 0.25, 
                        0.25, 
                        scale_est)
  ) %>% 
  ungroup() %>% 
  filter(year %in% year_seq)
eps_max <- quantile(sub_spatial$epsilon_st, 0.9999)

png(here::here("figs", "ms_figs_season", 
               "scaled_sp_preds.png"), 
    height = 8, width = 8, units = "in", res = 200)
ggplot() + 
  geom_raster(data = sub_spatial %>% filter(dataset == "summer"),
              aes(X, Y, fill = scale_est2)) +
  coord_fixed() +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Scaled\nAbundance",
    breaks = c(0, 0.05, 0.15, 0.25),
    labels = c("0", "0.05", "0.15", ">0.25")
  ) +
  facet_grid(year~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
dev.off()


png(here::here("figs", "ms_figs_season", "st_rf.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = sub_spatial %>% filter(dataset == "summer"),
              aes(X, Y, fill = epsilon_st)) +
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
omega_dat <- dat_tbl %>% 
  unnest(cols = "sp_preds") %>% 
  filter(survey_f == "hss",
         year == unique(dat$year)[1]) 
omega_max <- max(abs(omega_dat$omega_s))

png(here::here("figs", "ms_figs_season", "spatial_rf.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = omega_dat, aes(X, Y, fill = omega_s)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(-1 * omega_max, omega_max),
                       name = "Spatial\nField") +
  facet_grid(dataset~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
dev.off()


# FIXED EFFECTS ----------------------------------------------------------------

