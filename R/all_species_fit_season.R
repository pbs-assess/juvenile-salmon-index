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


# alternative branch
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "zeta-intercept")


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)


# downscale data and predictive grid
dat_in <- readRDS(here::here("data", "catch_survey_sbc.rds")) 


# make key so that missing year-season combinations can be indexed
year_season_key <- expand.grid(
  year = unique(dat_in$year) %>% sort(),
  season_f = unique(dat_in$season_f) %>% sort()
) %>% 
  arrange(year, season_f) %>% 
  mutate(
    # make consecutive index for model fitting
    ys_index = seq(1, nrow(.), by = 1),
    year_season_f = paste(year, season_f, sep = "_") %>% as.factor()
  )
saveRDS(year_season_key, here::here("data", "year_season_key.rds"))


# downscale data and predictive grid
dat <- dat_in %>% 
  mutate(
    year_f = as.factor(year),
    yday = lubridate::yday(date),
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    scale_dist = scale(as.numeric(dist_to_coast_km))[ , 1],
    scale_depth = scale(as.numeric(target_depth))[ , 1],
    year_season_f = paste(year_f, season_f, sep = "_") %>% as.factor(),
    day_night = as.factor(day_night)
  ) %>% 
  left_join(., year_season_key %>% select(ys_index, year_season_f)) %>% 
  filter(!species == "steelhead") %>% 
  droplevels()


# ggplot(dat %>% filter(!n_juv == "0")) +
#   geom_boxplot(aes(x = season_f, y = log(n_juv))) +
#   facet_wrap(~species)

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

# specify missing season-year levels
missing_surveys <- year_season_key$ys_index[!(year_season_key$ys_index %in% 
                                                dat$ys_index)]


dat_tbl <- readRDS(here::here("data", "fits", 
                              "all_spatial_varying_new_scale_iso.rds"))

# dat_tbl2 <- dat %>%
#   group_by(species) %>%
#   group_nest()

# fit
# fits_list_iso <- furrr::future_map(
#   dat_tbl2$data,
#   function(dat_in) {
#    sdmTMB(
#      n_juv ~ 0 + season_f + day_night + survey_f +
#        scale_depth + scale_dist,
#        # target_depth + dist_to_coast_km,
#      offset = dat_in$effort,
#      data = dat_in,
#      mesh = spde,
#      family = sdmTMB::nbinom2(),
#      spatial = "off",
#      spatial_varying = ~ 0 + season_f + year_f,
#      time_varying = ~ 1,
#      time_varying_type = "rw0",
#      time = "ys_index",
#      spatiotemporal = "off",
#      anisotropy = FALSE,
#      share_range = TRUE,
#      extra_time = missing_surveys,
#      priors = sdmTMBpriors(
#        phi = halfnormal(0, 10),
#        matern_s = pc_matern(range_gt = 25, sigma_lt = 10),
#        matern_st = pc_matern(range_gt = 25, sigma_lt = 10)
#      ),
#      control = sdmTMBcontrol(
#        newton_loops = 1,
#        map = list(
#          # 1 per season, fix all years to same value
#          ln_tau_Z = factor(
#            c(1, 2, 3, rep(4, times = length(unique(dat$year)) - 1))
#          )
#        )
#      ),
#      silent = FALSE
#    )
#   }
# )
# 
# dat_tbl2$model <- "tv_iso"
# dat_tbl2$fit <- fits_list_iso
# 
# saveRDS(dat_tbl2, here::here("data", "fits", "all_spatial_varying_new_scale_iso.rds"))

purrr::map(dat_tbl$fit, sanity)

# check residuals
dat_tbl$sims <- purrr::map(dat_tbl$fit, simulate, nsim = 50)
qq_list <- purrr::map2(dat_tbl$sims, dat_tbl$fit, dharma_residuals)


## PARAMETER ESTIMATES ---------------------------------------------------------

# fixed parameter estimates
fix_pars <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species),
  function(x, y) {
    tidy(x, effects = "fixed", conf.int = T) %>% 
      mutate(
        species = y
      ) 
  }
) %>% 
  bind_rows() %>% 
  filter(term %in% c("season_fsp", "season_fsu", "season_fwi", "survey_fipes", 
                     "day_nightNIGHT")) %>% 
  mutate(
    term2 = ifelse(grepl("season", term), "Season", term),
    term = fct_recode(as.factor(term), "Spring" = "season_fsp", 
                      "Summer" = "season_fsu", "Fall" = "season_fwi",
                      "IPES Survey" = "survey_fipes", 
                      "Nocturnal Sampling" = "day_nightNIGHT")
  )


shape_pal <- c(23, 24, 25)
names(shape_pal) <- c("Spring", "Summer", "Fall")

alpha_1 <- ggplot(
  fix_pars %>% filter(!term %in% c("IPES Survey", "Nocturnal Sampling")),
  aes(species, estimate, ymin = conf.low, ymax = conf.high, fill = species,
      shape = term)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = shape_pal, name = "") +
  geom_pointrange(position = position_dodge(width = 0.35)) +
  guides(fill = "none") +
  theme(axis.title = element_blank()) +
  facet_wrap(~term2)

alpha_2 <- ggplot(
  fix_pars %>% filter(term %in% c("IPES Survey", "Nocturnal Sampling")),
  aes(species, estimate, ymin = conf.low, ymax = conf.high, fill = species)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  geom_pointrange(position = position_dodge(width = 0.2),
                  shape = 21) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  facet_wrap(~term) +
  guides(fill = "none") +
  theme(axis.title = element_blank())

plot_g <- cowplot::plot_grid(alpha_1, alpha_2, ncol = 1)
y_grob <- grid::textGrob("Parameter Estimate", rot = 90)

png(here::here("figs", "ms_figs_season", "fix_ints.png"), height = 4, width = 6,
    units = "in", res = 200)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(plot_g, left = y_grob)
)
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
  headrope_depth = seq(0, 60, length.out = 100),
  scale_dist = 0,
  survey_f = unique(dat$survey_f),
  year = unique(dat$year)
) %>% 
  mutate(
    scale_depth = (headrope_depth - mean(dat$target_depth)) / 
      sd(dat$target_depth),
    year_f = as.factor(year),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    day_night == "DAY", season_f == "su", survey_f == "hss", year_f == "2012"
  ) %>% 
  left_join(
    ., year_season_key, by = c("season_f", "year")
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

depth_plot <- ggplot(
  data = depth_preds2,
  aes(x = headrope_depth, y = scale_est, ymin = scale_lo, ymax = scale_up, 
      fill = species)
) +
  ggsidekick::theme_sleek() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, nrow = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Target Headrope Depth (m)") +
  guides(fill = "none") +
  theme(axis.title.y = element_blank()) +
  ylab("Scaled Abundance Index") 


# distance to coast
dist_dat <- expand.grid(
  dist_to_coast_km = seq(0.1, max(dat$dist_to_coast_km), length.out = 100),
  scale_depth = 0,
  day_night = unique(dat$day_night),
  season_f = unique(dat$season_f),
  target_depth = 0,
  survey_f = unique(dat$survey_f),
  year = unique(dat$year)
  ) %>% 
  mutate(
    scale_dist = (dist_to_coast_km - mean(dat$dist_to_coast_km)) / 
      sd(dat$dist_to_coast_km),
    year_f = as.factor(year),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    day_night == "DAY", season_f == "su", survey_f == "hss", year_f == "2012"
  ) %>% 
  left_join(
    ., year_season_key, by = c("season_f", "year")
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

dist_plot <- ggplot(
  data = dist_preds2,
  aes(x = dist_to_coast_km, y = scale_est, ymin = scale_lo, ymax = scale_up, 
      fill = species)
) +
  ggsidekick::theme_sleek() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, nrow = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Distance to Coastline (km)") +
  theme(axis.title.y = element_blank()) +
  guides(fill = "none") 

plot_g2 <- cowplot::plot_grid(depth_plot, dist_plot, ncol = 1)
y_grob <- grid::textGrob("Scaled Abundance Index", rot = 90)

png(here::here("figs", "ms_figs_season", "smooth_effs.png"), height = 5, 
    width = 8,
    units = "in", res = 200)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(plot_g2, left = y_grob)
)
dev.off()



## SPATIAL GRID  ---------------------------------------------------------------

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

pred_grid_list <- rbind(
  year_season_key %>% mutate(survey_f = "hss"),
  year_season_key %>% mutate(survey_f = "ipes")
) %>% 
  mutate(
    year_f = as.factor(year),
    id = row_number(),
    survey_f = survey_f %>% as.factor()
  ) %>% 
  split(., .$id)


## INDICES ---------------------------------------------------------------------

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
        ys_index = x$ys_index,
        target_depth = 0,
        day_night = "DAY",
        scale_depth = (target_depth - mean(dat$target_depth)) / 
          sd(dat$target_depth),
        scale_dist = (dist_to_coast_km - mean(dat$dist_to_coast_km)) / 
          sd(dat$dist_to_coast_km)
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


# predictions for index (set to HSS reference values)
index_grid_hss <- index_grid %>% filter(survey_f == "hss")

ind_preds <- purrr::map(
  dat_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid_hss,
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)



# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m)
sp_scalar <- 1000^2 * 13

index_list <- purrr::map(
  ind_preds,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
saveRDS(index_list, here::here("data", "fits", "season_index_list.rds"))


index_list <- readRDS(
  here::here("data", "fits", "season_index_list.rds")
)

# define years where survey occurred
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


index_dat <- purrr::map2(
  dat_tbl$species, index_list,
  ~ .y %>% 
    mutate(species = .x) 
) %>%
  bind_rows() %>% 
  left_join(., year_season_key, by = "ys_index") %>% 
  filter(!season_f == "sp") %>% 
  mutate(
    season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
    survey = case_when(
      season == "fall" & year %in% fall_years ~ "sampled",
      season == "summer" & year %in% summer_years ~ "sampled",
      TRUE ~ "no survey"
    ),
    log_lwr = log_est - (1.96 * se),
    log_upr = log_est + (1.96 * se)
  ) 

## index plots
shape_pal <- c(21, 23)
names(shape_pal) <- unique(index_dat$survey)
# log abundance 
log_index <- ggplot(index_dat, aes(year, log_est)) +
  geom_pointrange(aes(ymin = log_lwr, ymax = log_upr, fill = species, 
                      shape = survey)) +
  geom_hline(aes(yintercept = mean_log_est), lty = 2) +
  labs(x = "Year", y = "Log Abundance Index") +
  scale_shape_manual(values = shape_pal) +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

png(here::here("figs", "ms_figs_season", "log_index.png"), 
    height = 8, width = 8, units = "in", res = 200)
log_index
dev.off()


## scaled abundance with no intervals
scaled_index <- index_dat %>% 
  filter(survey == "sampled")%>%
  group_by(species, season) %>%
  mutate(
    mean_log_est = mean(log_est),
    max_est = max(est),
    scale_est = est / max_est,
    scale_lwr = lwr / max_est,
    scale_upr = upr / max_est
  ) %>%
  ungroup() %>% 
  ggplot(., aes(year, scale_est)) +
  geom_pointrange(aes(ymin = scale_lwr, ymax = scale_upr, fill = species,
                      shape = survey)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(x = "Year", y = "Scaled Abundance Index") +
  scale_shape_manual(values = shape_pal) +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

png(here::here("figs", "ms_figs_season", "scaled_index.png"), 
    height = 8, width = 8, units = "in", res = 200)
scaled_index
dev.off()


# check among season correlations
index_dat %>% 
  split(., .$species) %>% 
  purrr::map(., function (x) {
    fall_dat <- x %>% filter(season == "fall") %>% pull(log_est)
    summ_dat <- x %>% filter(season == "summer") %>% pull(log_est)
    cor(fall_dat, summ_dat)
  })


# compare HSS index to one where survey design is not accounted for
index_grid_true <- index_grid %>% filter(fake_survey == "0")
ind_preds_true <- purrr::map(
  dat_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid_true,
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)
index_true_survey_list <- purrr::map(
  ind_preds_true,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
saveRDS(index_true_survey_list, 
        here::here("data", "fits", "season_index_true_survey_list.rds"))

index_dat2 <- purrr::map2(
  dat_tbl$species, index_true_survey_list,
  ~ .y %>% 
    mutate(species = .x) 
) %>%
  bind_rows() %>% 
  left_join(., year_season_key, by = "ys_index") %>% 
  mutate(preds = "obs")

index_dat %>% 
  mutate(preds = "hss") %>% 
  select(colnames(index_dat2)) %>% 
  rbind(., 
        index_dat2 ) %>% 
  filter(season_f == "su",
         year > 2017) %>% 
  group_by(preds, species) %>% 
  summarize(mean_est = mean(est)) %>% 
  pivot_wider(names_from = preds, values_from = mean_est) %>% 
  mutate(diff = obs / hss)


# SPATIAL PREDS ----------------------------------------------------------------

# shape file for coastline
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., xmin = -129, ymin = 48.25, xmax = -124, ymax = 51.2) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))


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
        ys_index = x$ys_index,
        target_depth = 0,
        day_night = "DAY"
      )
  }) %>%
  bind_rows() %>% 
  filter(
    !survey_f == "ipes"
  ) %>% 
  select(-c(depth, slope, shore_dist)) 

spatial_preds_list <- purrr::map2(
  dat_tbl$fit,
  dat_tbl$species,
  ~ {
    predict(.x,
            newdata = spatial_grid,
            se_fit = FALSE, re_form = NULL) %>% 
      mutate(species = .y)
  }
)
spatial_preds_list <- spatial_preds
spatial_preds <- spatial_preds_list %>%
  bind_rows() %>% 
  filter(!season_f == "sp")
saveRDS(spatial_preds,
        here::here("data", "fits", "all_spatial_varying_new_sp_preds.rds"))


plot_map_raster <- function(dat, column = est) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis_c()
}


# spatiotemporal and total preds random fields for subset of years
year_seq <- seq(1999, 2019, by = 5)
sub_spatial <- spatial_preds %>% 
  filter(year %in% year_seq) %>% 
  group_by(species) %>% 
  mutate(
    grid_est = sp_scalar * exp(est),
    scale_est = grid_est / max(grid_est)
  ) %>% 
  ungroup() 
scale_max <- quantile(sub_spatial$scale_est, 0.99)
sub_spatial$scale_est2 <- ifelse(
  sub_spatial$scale_est > scale_max, 
  scale_max, 
  sub_spatial$scale_est
)

png(here::here("figs", "ms_figs_season", 
               "scaled_sp_preds_su.png"), 
    height = 8, width = 8, units = "in", res = 200)
ggplot() + 
  geom_raster(data = sub_spatial %>% filter(season_f == "su"),
              aes(X, Y, fill = scale_est2)) +
  coord_fixed() +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Scaled\nAbundance"
  ) +
  facet_grid(year~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
dev.off()


png(here::here("figs", "ms_figs_season", 
               "scaled_sp_preds_fall.png"), 
    height = 8, width = 8, units = "in", res = 200)
ggplot() + 
  geom_raster(data = sub_spatial %>% filter(season_f == "wi"),
              aes(X, Y, fill = scale_est2)) +
  coord_fixed() +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Scaled\nAbundance"
  ) +
  facet_grid(year~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
dev.off()

year_field_seq <- paste("zeta_s_year_f", year_seq, sep = "")
omega_yr <- sub_spatial %>% 
  # values will be duplicated across years, seasons and surveys; select one
  filter(year_f == year_seq[1], season_f == "su") %>% 
  select(X, Y, species, year_field_seq) %>% 
  pivot_longer(cols = starts_with("zeta"), values_to = "omega_est", 
             names_to = "year", names_prefix = "zeta_s_year_f")

png(here::here("figs", "ms_figs_season", "year_rf.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = omega_yr,
              aes(X, Y, fill = omega_est)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_distiller(palette = "Spectral", 
                       # limits = c(-1 * eps_max, eps_max),
                       name = "Spatial\nField") +
  facet_grid(year~species) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
dev.off()


# spatial random effects by season
omega_season <- sub_spatial %>% 
  # values will be duplicated across years, seasons and surveys; select one
  filter(year_f == year_seq[1], season_f == "wi") %>% 
  select(X, Y, species, starts_with("zeta_s_season_")) %>%
  pivot_longer(cols = starts_with("zeta"), values_to = "omega_est", 
               names_to = "season", names_prefix = "zeta_s_season_f") %>% 
  filter(!season == "sp")

png(here::here("figs", "ms_figs_season", "seasonal_rf.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = omega_season, aes(X, Y, fill = omega_est)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(season~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
dev.off()
