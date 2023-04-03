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

# identify viable years by season
summer_years <- dat %>%
  filter(season_f == "su") %>%
  pull(year) %>%
  unique() %>% 
  sort()
fall_years <- dat %>%
  filter(season_f == "wi") %>%
  pull(year) %>%
  unique() %>% 
  sort()

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


dat_tbl <- readRDS(here::here("data", "fits", "top_st_mod_all_sp_season.rds"))


# separate seasons and join
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
# 
# dat_tbl <- rbind(summer_tbl, fall_tbl) %>%
#   group_by(species, dataset) %>%
#   group_nest() %>% 
#   mutate(
#     anisotropy = ifelse(species == "chinook" & dataset == "summer", FALSE, TRUE),
#     time_model = ifelse(species == "sockeye" & dataset == "fall", "rw", "ar1")
#   )


# dat_tbl$fits <- purrr::pmap(
#   list(dat_tbl$dataset, dat_tbl$time_model, dat_tbl$anisotropy, dat_tbl$data),
#   function(x, y, z, dat_in) {
#     dum <- dat_in %>% droplevels()
#     dat_coords <- dum %>% 
#       select(utm_x_1000, utm_y_1000) %>% 
#       as.matrix()
#     
#     ## use INLA mesh based on SA recommendations and model selection (see notes)
#     inla_mesh_raw <- INLA::inla.mesh.2d(
#       loc = dat_coords,
#       max.edge = c(1, 5) * 500,
#       cutoff = 20,
#       offset = c(20, 200)
#     ) 
#     spde <- make_mesh(
#       dum,
#       c("utm_x_1000", "utm_y_1000"),
#       mesh = inla_mesh_raw
#     ) 
#     
#     extra_time <- NULL
#     if (x == "fall") {
#       extra_time <-  c(2018#, 2020
#       )
#       formula_in <- as.formula(
#         paste("n_juv ~ target_depth + day_night")
#       )
#     }
#     if (x == "summer") {
#       extra_time <-  c(2016, 2020
#                        #, 2021
#       )
#       formula_in <- as.formula(
#         paste("n_juv ~ survey_f + target_depth + day_night")
#       )
#     }
#     sdmTMB(
#       formula_in,
#       offset = dum$effort,
#       data = dum,
#       mesh = spde,
#       family = sdmTMB::nbinom2(),
#       spatial = "on",
#       spatiotemporal = y,
#       time = "year",
#       anisotropy = z,  
#       share_range = TRUE,
#       extra_time = extra_time,
#       control = sdmTMBcontrol(
#         newton_loops = 1
#       ),
#       silent = FALSE
#     )
#   } 
# )
# 
# dat_tbl$fits <- fits_list

# sub_tbl <- dat_tbl %>%
#   filter(
#     (species == "sockeye" & dataset == "fall" & model == "rw") |
#       (!(species == "sockeye" & dataset == "fall") & model == "ar1")
#   )
# # saveRDS(sub_tbl, here::here("data", "fits", "top_st_mod_all_sp_season.rds"))


dat_tbl %>% 
  filter(dataset == "summer") %>% 
  pull(fits) %>% 
  purrr::map(., sanity)


# check residuals
sub_tbl$sims <- purrr::map(sub_tbl$fits, simulate, nsim = 50)
qq_list <- purrr::map2(sub_tbl$sims, sub_tbl$fits, dharma_residuals)


## PARAMETER ESTIMATES ---------------------------------------------------------

# fixed parameter estimates
fix_pars <- purrr::pmap(
  list(dat_tbl$fits, dat_tbl$species, dat_tbl$dataset),
  function(x, y, z) {
    tidy(x, effects = "fixed", conf.int = T) %>% 
      mutate(
        species = y,
        season = z
      ) 
  }
) %>% 
  bind_rows() %>% 
  filter(term %in% c("survey_fipes", "day_nightNIGHT")) %>% 
  mutate(
    term = fct_recode(as.factor(term), "IPES Survey" = "survey_fipes", 
                      "Nocturnal Sampling" = "day_nightNIGHT"),
    season = fct_relevel(season, "summer", "fall")
  )

png(here::here("figs", "ms_figs_season", "fix_ints.png"), height = 3, width = 4,
    units = "in", res = 200)
ggplot(
  fix_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high, fill = species,
      shape = season)
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
  list(dat_tbl$fits, dat_tbl$species, dat_tbl$dataset),
  function(x, y, z) {
    tidy(x, effects = "ran_par", conf.int = T) %>% 
      mutate(
        species = y,
        season = z
      ) 
  }
) %>% 
  bind_rows() %>% 
  mutate(
    season = fct_relevel(season, "summer", "fall")
  )

png(here::here("figs", "ms_figs_season", "ran_pars.png"), height = 4, width = 8,
    units = "in", res = 200)
ggplot(
  ran_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high,
      fill = species, shape = season)
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
target_dat <- data.frame(
  day_night = "DAY",
  target_depth = seq(min(dat$target_depth), 
                     max(dat$target_depth), length.out = 100),
  survey_f = "hss",
  year = 2012
) 
depth_preds <- purrr::pmap(
  list(dat_tbl$fits, dat_tbl$species, dat_tbl$dataset),
  function(x, y, z) {
    predict(x, newdata = target_dat, se_fit = T, re_form = NA) %>% 
      mutate(species = y,
             season = z)
  }
) %>%
  bind_rows()



depth_preds2 <- depth_preds %>% 
  mutate(season = fct_relevel(season, "summer", "fall")) %>% 
  group_by(species, season) %>% 
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

png(here::here("figs", "ms_figs_season", "depth_preds.png"), height = 4, width = 8,
    units = "in", res = 200)
ggplot(
  data = depth_preds2,
  aes(x = target_depth, y = scale_est, ymin = scale_lo, 
  ymax = scale_up, fill = species)
) +
  ggsidekick::theme_sleek() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(season~species) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Target Headrope Depth (m)") +
  guides(fill = "none") +
  ylab("Abundance Index")
dev.off()




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

# scalar for spatial predictions; since preds are in m3, multiply by
# (1000 * 1000 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m)
sp_scalar <- 1000^2 * 13


## INDEX -----------------------------------------------------------------------

index_dat <- readRDS(here::here("data", "fits", "season_index.rds")) %>% 
  filter(!(year == "2021" & season == "summer")) %>% 
  mutate(
    log_lwr = log_est + (qnorm(0.025) * se),
    log_upr = log_est + (qnorm(0.975) * se)
  ) %>% 
  group_by(species, season) %>% 
  mutate(
    log_mean = mean(log_est)
  ) %>% 
  ungroup()


# summer high seas
# ind_preds_sum <- purrr::map(
#   sub_tbl %>% filter(dataset == "summer") %>% pull(fits),
#   ~ {
#     predict(.x,
#             newdata = exp_grid %>% filter(survey_f == "hss", week == "25"),
#             return_tmb_object = TRUE)
#   }
# )
# index_list_sum <- furrr::future_map(
#   ind_preds_sum, 
#   get_index, 
#   area = sp_scalar, 
#   bias_correct = TRUE
# )
# index_list_sum_trim <- index_list_sum[c(2,4,6,8,9)]
# 
# # fall high seas
# ind_preds_fall <- purrr::map(
#   sub_tbl %>% filter(dataset == "fall") %>% pull(fits),
#   ~ {
#     predict(.x,
#             newdata = exp_grid %>% filter(survey_f == "hss", week == "42"),
#             return_tmb_object = TRUE)
#   }
# )
# index_list_fall <- furrr::future_map(
#   ind_preds_fall, get_index, 
#   area = sp_scalar, 
#   bias_correct = TRUE
# )
# 
# index_sum <- purrr::map2(
#   index_list_sum_trim, 
#   sub_tbl %>% filter(dataset == "summer") %>% pull(species),
#   function (x, sp) {
#     x$species <- sp
#     x$season <- "su"
#     return(x)
#   }) %>% 
#   bind_rows()
# index_fall <- purrr::map2(
#   index_list_fall,   
#   sub_tbl %>% filter(dataset == "fall") %>% pull(species),
#   function (x, sp) {
#     x$species <- sp
#     x$season <- "fa"
#     return(x)
#   }) %>% 
#   bind_rows()
# 
# 
# index_dat <- rbind(index_sum %>% filter(year %in% summer_years), 
#                    index_fall %>% filter(year %in% fall_years)) %>% 
#   mutate(
#     season = factor(season, levels = c("su", "fa"), 
#                     labels = c("summer", "fall"))
#   ) %>% 
#   group_by(species) %>% 
#   mutate(
#     max_est = max(est),
#     scale_est = est / max_est,
#     scale_lwr = lwr / max_est,
#     scale_upr = upr / max_est
#   ) %>% 
#   ungroup()
# 
# saveRDS(index_dat, here::here("data", "fits", "season_index.rds"))


## scaled abundance with no intervals
log_index <- ggplot(index_dat, aes(year, log_est)) +
  geom_pointrange(aes(ymin = log_lwr, ymax = log_upr, fill = species),
                  shape = 21) +
  geom_hline(aes(yintercept = log_mean), lty = 2) +
  labs(x = "Year", y = "Log Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none")


pdf(here::here("figs", "sopo_figs", "sp_preds_2022_v2.pdf"))
log_index
dev.off()


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

dat_tbl$exp_grid <- purrr::map(
  dat_tbl$dataset, 
  ~ {
    if (.x == "summer") {
      exp_grid %>% filter(week == "25")
      } else {
        exp_grid %>% filter(week == "42")  
      }
  }
)

# make spatial 
spatial_preds <- furrr::future_map2(
  dat_tbl$fits, dat_tbl$exp_grid, function (x, y) {
    predict(x, newdata = y, se_fit = FALSE, re_form = NULL)
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

