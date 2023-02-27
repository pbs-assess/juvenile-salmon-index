### Juvenile all species fit 
## Use chinook_fit as a template to fit different species models
## 1) Fit spatial model for each species to ensure reasonable convergence and
## to test for effects of survey domain on fixed effect estimates 
## (all_species_fit_spatial.R)
## 2) Fit saturated spatiotemporal model to each species
## 3) Calculate index for summer and fall (assuming surface and day tow)
## 4) Calculate fixed effects for spatial covariates
## 5) Calculate spatiotemporal effects (maps by year)
## Oct 13, 2022


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

dat_tbl <- dat  %>% 
  group_by(species) %>% 
  group_nest() %>% 
  mutate(
    anisotropy = ifelse(species == "pink", FALSE, TRUE)
  )

## coordinates should be common to all species
dat_coords <- dat %>% 
  filter(species == "pink") %>% 
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
  dat %>% filter(species == "pink"),
  c("utm_x_1000", "utm_y_1000"),
  mesh = inla_mesh_raw
) 

## plotting color palette
col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')

# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)



### EXP PLOTS ------------------------------------------------------------------

ggplot(dat) +
  geom_boxplot(aes(x = target_depth_bin, y = log(n_juv))) +
  facet_wrap(~species, scales = "free_x")

ggplot(dat) +
  geom_point(aes(x = bath_depth_mean_m, y = log(n_juv))) +
  facet_wrap(~species)

ggplot(dat) +
  geom_point(aes(x = dist_to_coast_km, y = log(n_juv))) +
  facet_wrap(~species)

ggplot(dat) +
  geom_point(aes(x = week, y = log(n_juv))) +
  facet_wrap(~species)

ggplot(dat) +
  geom_boxplot(aes(x = day_night, y = log(n_juv))) +
  facet_wrap(~species)


### HISTOGRAMS OF CATCH --------------------------------------------------------

catch_hist <- dat %>% 
  filter(n_juv > 0) %>% 
  group_by(species) %>% 
  mutate(median_catch = mean(n_juv)) %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_histogram(aes(x = n_juv, fill = species), colour = "black",
                 bins = 30) +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(trans='log10') +
  geom_vline(aes(xintercept = median_catch), lty = 2, colour = "red") +
  facet_wrap(~species, ncol = 1) +
  labs(x = "Number of Individuals", y = "Number of Tows") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none")


png(here::here("figs", "ms_figs", "catch_histogram.png"), height = 8.5, 
    width = 4, units = "in", res = 200)
catch_hist
dev.off()


### FIT SATURATED --------------------------------------------------------------

st_mod_ar1 <- furrr::future_map2(
  dat_tbl$data, dat_tbl$anisotropy,
  ~ {
    sdmTMB(
      n_juv ~ 0 +
        year_f +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = .x$effort,
      data = .x,
      mesh = spde,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "ar1",
      time = "year",
      anisotropy = .y, #TRUE,
      share_range = FALSE,
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      ),
      silent = FALSE
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)

purrr::map(st_mod_ar1, sanity)
## all look good


dat_tbl$st_mod <- st_mod_ar1
saveRDS(dat_tbl, here::here("data", "fits", "st_mod_all_sp_ar1.rds"))


dat_tbl <- readRDS(here::here("data", "fits", "st_mod_all_sp_ar1.rds"))

purrr::map(
  dat_tbl$st_mod, sanity
)


## check residuals
dat_tbl$sims <- purrr::map(dat_tbl$st_mod, simulate, nsim = 100)
qq_list <- purrr::map2(dat_tbl$sims, dat_tbl$st_mod, dharma_residuals)

samps <- sample.int(250, size = 9)
hist_list <- purrr::pmap(
  list(dat_tbl$data, dat_tbl$st_mod, dat_tbl$sims), function (x, y, z) {
    max_n <- 250#quantile(exp(x$n_juv), 0.8) %>% as.numeric
    breaks_vec <- c(seq(0, max_n, by = 10))
    obs_dat <- x %>% filter(n_juv < max_n & n_juv > 0) %>% pull(n_juv)
    par(mfrow = c(3, 3))
    for (i in seq_along(samps)) {
      dum <- z[, i ]
      hist(obs_dat, col="green", pch=20, cex=4, breaks=breaks_vec)
      hist(dum[dum < max_n & dum > 0], pch=20, cex=4, breaks=breaks_vec, 
           col=rgb(1,0,0,0.5), add=TRUE)
    }
  }
)

pred_obs_list <- purrr::pmap(
  list(dat_tbl$data, dat_tbl$sims), function (x, y) {
    for (i in seq_along(samps)) {
      dum <- y[, i ]
      plot(log(dum) ~ log(x$n_juv))
      abline(0, 1, col = "red")
    }
  }
)

resid_list <- purrr::pmap(
  list(dat_tbl$data, dat_tbl$species, dat_tbl$st_mod), function (x, y, z) {
    x$resids <- residuals(z)
    ggplot(x, aes(utm_x_1000, utm_y_1000, col = resids)) +
      scale_colour_gradient2() +
      geom_point() +
      facet_wrap(~year_f) +
      coord_fixed() +
      labs(title = y) +
      ggsidekick::theme_sleek()
  }
)

pdf(here::here("figs", "diagnostics", "spatial_resids.pdf"))
resid_list
dev.off()


## MAKE FE PREDICTIONS ---------------------------------------------------------

# make conditional predictive dataframes
week_dat <- data.frame(
 week = seq(0, 52, length.out = 100),
 dist_to_coast_km = median(dat$dist_to_coast_km),
 day_night = "DAY",
 target_depth = 0,
 survey_f = "hss",
 year = 2012
 ) %>% 
  mutate(year_f = as.factor(year))
dist_dat <- data.frame(
  week = median(dat$week),
  dist_to_coast_km = seq(min(dat$dist_to_coast_km), 
                         max(dat$dist_to_coast_km), length.out = 100),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year = 2012L
  ) %>% 
  mutate(year_f = as.factor(year))
target_dat <- data.frame(
  week = median(dat$week),
  dist_to_coast_km = median(dat$dist_to_coast_km),
  day_night = "DAY",
  target_depth = seq(min(dat$target_depth), 
                         max(dat$target_depth), length.out = 100),
  survey_f = "hss",
  year = 2012L
  ) %>% 
  mutate(year_f = as.factor(year))
# day_dat <- data.frame(
#   week = median(dat$week),
#   dist_to_coast_km = median(dat$dist_to_coast_km),
#   day_night = c("DAY", "NIGHT"),
#   target_depth = 0,
#   survey_f = "hss",
#   year = 2012L
#   ) %>% 
#   mutate(year_f = as.factor(year))
# survey_dat <- data.frame(
#   week = median(dat$week),
#   dist_to_coast_km = median(dat$dist_to_coast_km),
#   day_night = "DAY",
#   target_depth = 0,
#   survey_f = c("hss", "ipes"),
#   year = 2012L
#   ) %>% 
#   mutate(year_f = as.factor(year))

pred_tbl <- tibble(
  var = c("week",
          "dist_to_coast_km",
          "target_depth"#,
          # "day_night",
          # "survey_f"
          ),
  data = list(week_dat, dist_dat, target_dat#, day_dat, survey_dat
              ),
  plot = c("line", "line", "line"#, "dot", "dot"
           )
)

# predict
pred_list <- vector(length = nrow(pred_tbl), mode = "list")
for (i in seq_along(pred_tbl$var)) {
  pred_list[[i]] <- furrr::future_map2(
    dat_tbl$st_mod, dat_tbl$species, function(x , sp) {
      predict(x, newdata = pred_tbl$data[[i]], se_fit = T, re_form = NA) %>% 
        mutate(species = sp)
      }
  ) %>%
    bind_rows()
}

pred_tbl$pred_dat <- pred_list


# takes a long time so save
saveRDS(pred_tbl, here::here("data", "preds", "fe_preds.rds"))
pred_tbl <- readRDS(here::here("data", "preds", "fe_preds.rds"))


pred_tbl$pred_dat <- purrr::map(
  pred_tbl$pred_dat,
  ~ {
    .x %>% 
      group_by(species) %>% 
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
  }
)

# make individual plots for ms
p_dat <- pred_tbl %>% 
  select(var, pred_dat) %>% 
  unnest(cols = pred_dat) %>% 
  mutate(species = tolower(species),
         day_night = tolower(day_night),
         survey_f = toupper(survey_f)) 


pp_foo <- function(dat, ...) {
  ggplot(
    dat,
    mapping = aes(!!!ensyms(...))
  ) +
    ggsidekick::theme_sleek() +
    scale_fill_manual(values = col_pal) +
    ylab("Abundance Index")
}

png(here::here("figs", "ms_figs", "week_preds.png"), height = 8.5, width = 4,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "week"),
  x = "week", y = "scale_est", ymin = "scale_lo", 
  ymax = "scale_up", fill = "species"
) +
  ggsidekick::theme_sleek() +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, ncol = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Week") 
dev.off()


png(here::here("figs", "ms_figs", "depth_preds.png"), height = 8.5, width = 4,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "target_depth"),
  x = "target_depth", y = "scale_est", ymin = "scale_lo", 
  ymax = "scale_up", fill = "species"
) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, ncol = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Target Headrope Depth (m)") 
dev.off()


png(here::here("figs", "ms_figs", "dist_preds.png"), height = 8.5, width = 4,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "dist_to_coast_km"),
  x = "dist_to_coast_km", y = "scale_est", ymin = "scale_lo", 
  ymax = "scale_up", fill = "species"
) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~species, ncol = 1) +
  coord_cartesian(y = c(0, 2.5)) +
  xlab("Distance to Nearest Coastline (km)") 
dev.off()



## PARAMETER ESTIMATES ---------------------------------------------------------

## survey and DN estimates
# fixed parameter estimates
fix_pars <- purrr::map2(
  dat_tbl$st_mod, dat_tbl$species,
  ~ tidy(.x, effects = "fixed", conf.int = T) %>% 
    mutate(
      species = .y
    )
) %>% 
  bind_rows() %>% 
  filter(term %in% c("survey_fipes", "day_nightNIGHT")) %>% 
  mutate(
    term = fct_recode(as.factor(term), "IPES Survey" = "survey_fipes", 
                                "Nocturnal Sampling" = "day_nightNIGHT")
    )

png(here::here("figs", "ms_figs", "fix_ints.png"), height = 3, width = 4,
    units = "in", res = 200)
ggplot(
  fix_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high, fill = species)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  geom_pointrange(shape = 21) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  facet_wrap(~term, nrow = 2) +
  labs(y = "Parameter Estimate") +
  theme(legend.position = "none",
        axis.title.x = element_blank())
dev.off()


## random parameter estimates
ran_pars <- purrr::map2(
  dat_tbl$st_mod, dat_tbl$species,
  ~ tidy(.x, effects = "ran_par", conf.int = T) %>% 
    mutate(
      species = .y
    )
) %>% 
  bind_rows() %>% 
  # add unique identifier for second range term
  group_by(species, term) %>% 
  mutate(par_id = row_number(),
         term = ifelse(par_id > 1, paste(term, par_id, sep = "_"), term),
         term = case_when(
           term == "range" ~ "range_O",
           term == "range_2" ~ "range_E",
           TRUE ~ term
         )) %>% 
  ungroup()

png(here::here("figs", "ms_figs", "ran_pars.png"), height = 3, width = 7,
    units = "in", res = 200)
ggplot(
  ran_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high,
      fill = species)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  ylab("Parameter Estimate") +
  geom_pointrange(shape = 21) +
  facet_wrap(~term, scales = "free_y") +
  theme(
    legend.position = "none"
  )
dev.off()


## MAKE SPATIAL PREDICTIONS ----------------------------------------------------

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
  filter(season_f == "wi") %>%
  pull(year) %>% 
  unique()


# fixed effects plots
fall_tbl <- spatial_pred_tbl %>%
  filter(season == "fall")
for (i in seq_along(fall_tbl$species)) {
  .x <- fall_tbl$spatial_preds[[i]] %>% 
    filter(
      year %in% fall_years
    ) %>% 
    mutate(
      scale_est = exp(est) / max(exp(est))
    )
  .y <- fall_tbl$species[[i]]
  
  # total effects
  max_est <- quantile(.x$scale_est, 0.999)
  p <-  ggplot() + 
    geom_raster(data = .x, aes(X, Y, fill = scale_est)) +
    coord_fixed() +
    geom_sf(data = coast, color = "black", fill = "white") +
    ggsidekick::theme_sleek() +
    scale_fill_viridis_c(
      trans = "sqrt",
      limits = c(0, max_est)
    ) +
    facet_wrap(~year) +
    theme(axis.title = element_blank(),
          axis.text = element_blank())
  png(here::here("figs", "ms_figs", "fall_fe_preds",
                 paste(.y, "fall_pred.png", sep = "_")), 
      height = 8, width = 8, units = "in", res = 200)
  print(p)
  dev.off() 
  
  # epsilon effects
  epsilon_max <- max(abs(.x$epsilon_st))
  q <- ggplot() +
    geom_raster(data = .x, aes(X, Y, fill = epsilon_st)) +
    coord_fixed() +
    geom_sf(data = coast, color = "black", fill = "white") +
    ggsidekick::theme_sleek() +
    scale_fill_distiller(palette = "Spectral", 
                         limits = c(-1 * epsilon_max, epsilon_max)) +
    facet_wrap(~year) +
    theme(axis.title = element_blank(),
          axis.text = element_blank())
  png(here::here("figs", "ms_figs", "fall_fe_preds",
                 paste(.y, "fall_eps.png", sep = "_")), 
      height = 8, width = 8, units = "in", res = 200)
  print(q)
  dev.off() 
}


# spatial random effects by species and season
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
  # facet_grid(season~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1, 'cm')
  )
dev.off()


## SOPO MAPS -------------------------------------------------------------------

# use subset of predictions above to make maps of 2022 preds (FEs) and random 
# effects

fall_22 <- fall_tbl %>% 
  select(-week) %>% 
  unnest(cols = spatial_preds) %>% 
  filter(
    year %in% fall_years
  ) %>% 
  group_by(species) %>% 
  mutate(
    scale_est = exp(est) / max(exp(est))
  ) %>% 
  filter(
    year == "2022"
  )

max_est <- quantile(fall_22$scale_est, 0.999)

png(here::here("figs", "ms_figs", "fixed_preds_2022.png"), 
    height = 6, width = 7.5, units = "in", res = 200)
ggplot() + 
  geom_raster(data = fall_22, aes(X, Y, fill = scale_est)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    limits = c(0, max_est),
    name = "Scaled\nAbundance"
  ) +
  facet_wrap(~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.key.size = unit(1, 'cm'))
dev.off()


## INDICES ---------------------------------------------------------------------


# fix to HSS survey (can't combine because predictions shouldn't be passed
# duplicates, but require tmb_object stored in preds)
ind_preds_sum <- purrr::map(
  dat_tbl$st_mod,
  ~ {
    predict(.x,
            newdata = exp_grid %>% filter(survey_f == "hss", week == "25"),
            return_tmb_object = TRUE)
  }
)
index_list_sum <- purrr::map(ind_preds_sum, get_index, bias_correct = TRUE)

ind_preds_fall <- purrr::map(
  dat_tbl$st_mod,
  ~ {
    predict(.x,
            newdata = exp_grid %>% filter(survey_f == "hss", week == "42"),
            return_tmb_object = TRUE)
  }
)
index_list_fall <- purrr::map(ind_preds_fall, get_index, bias_correct = TRUE)

index_lists <- c(index_list_sum,
                 index_list_fall)
saveRDS(index_lists, here::here("data", "fits", "index_list.rds"))
# index_lists <- readRDS(here::here("data", "fits", "index_list.rds")) 


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
  group_by(species, season) %>% 
  mutate(
    max_est = max(est),
    scale_est = est / max_est,
    scale_lwr = lwr / max_est,
    scale_upr = upr / max_est
    ) %>% 
  ungroup()


index_scaled <- ggplot(index_dat, 
                     aes(year, scale_est)) +
  geom_pointrange(aes(ymin = scale_lwr, ymax = scale_upr, fill = species), 
                  shape = 21) +
  labs(x = "Year", y = "Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  coord_cartesian(y = c(0, 2.5)) +
  theme(legend.position = "none")

  
log_index_plot <- ggplot(index_dat, 
       aes(year, log_est)) +
  geom_pointrange(aes(ymin = log(lwr), ymax = log(upr), fill = species), 
                  shape = 21) +
  labs(x = "Year", y = "Log Abundance Index") +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  facet_grid(species~season, scales = "free_y") +
  theme(legend.position = "none")


png(here::here("figs", "ms_figs", "log_hss_index.png"))
log_index_plot
dev.off()

png(here::here("figs", "ms_figs", "hss_index_scaled.png"))
index_scaled
dev.off()



## correlations in indices
summer_cor <- index_dat %>% 
  filter(season == "summer") %>% 
  select(species, year, est) %>% 
  pivot_wider(names_from = species, 
              values_from = est) %>%
  select(-year) %>% 
  cor(., use = "complete.obs")  %>% 
  ggcorrplot::ggcorrplot(., 
                         # hc.order = TRUE, 
                         type = "lower",
                       lab = TRUE,
                       title = "summer") + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank())

fall_cor <- index_dat %>% 
  filter(season == "fall") %>% 
  select(species, year, est) %>% 
  pivot_wider(names_from = species, 
              values_from = est) %>%
  select(-year) %>% 
  cor(., use = "complete.obs")  %>% 
  ggcorrplot::ggcorrplot(., 
                         # hc.order = TRUE, 
                         type = "lower",
                         lab = TRUE,
                         title = "fall") + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank())

png(here::here("figs", "ms_figs", "cor_plot.png"), height = 4, width = 7.5,
    res = 250, units = "in")
cowplot::plot_grid(summer_cor,
                   fall_cor,
                   ncol = 2)
dev.off()


# separate survey effects (only summer since IPES not a fall survey)
ind_preds2 <- purrr::map(dat_tbl$st_mod, function (x) {
  predict(x, newdata = exp_grid %>% filter(fake_survey == "0", week == "25"),
          return_tmb_object = TRUE)
})
index_list2 <- purrr::map(ind_preds2, get_index, bias_correct = TRUE)
index_df2 <- purrr::map2(index_list2, dat_tbl$species, function (x, sp) {
  x$species <- sp
  return(x)
}) %>%
  bind_rows()


index1 <- index_dat %>% 
  filter(season == "summer") %>%
  # remove so DFs can be bound, then recalculate
  select(-c(season, max_est:scale_upr)) %>% 
  mutate(survey = "hss") 
index_combined <- index_df2 %>%
  mutate(survey = "ipes") %>%
  filter(year > 2016) %>%
  rbind(., index1) %>%
  filter(year %in% summer_years) %>% 
  group_by(species) %>% 
  mutate(
    max_est = max(est),
    log_lwr = log_est + (qnorm(0.025) * se),
    log_upr = log_est + (qnorm(0.975) * se),
    scale_est = est / max_est,
    scale_lwr = lwr / max_est,
    scale_upr = upr / max_est
  ) %>% 
  ungroup()

comb_index_scaled <- ggplot(index_combined, aes(year, scale_est)) +
  geom_pointrange(aes(ymin = scale_lwr, ymax = scale_upr, fill = survey),
                  shape = 21) +
  labs(x = "Year", y = "Scaled Count") +
  ggsidekick::theme_sleek() +
  facet_wrap(~species, scales = "free_y")

comb_index_log <- ggplot(index_combined, aes(year, log_est)) +
  geom_pointrange(aes(ymin = log_lwr, ymax = log_upr, fill = survey),
                  shape = 21) +
  labs(x = "Year", y = "Log Count") +
  ggsidekick::theme_sleek() +
  facet_wrap(~species, scales = "free_y")

comb_index <- ggplot(index_combined, aes(year, est)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr, fill = survey),
                  shape = 21) +
  labs(x = "Year", y = "Count") +
  ggsidekick::theme_sleek() +
  facet_wrap(~species, scales = "free_y")


pdf(here::here("figs", "ms_figs", "st_index_surv_all_sp.pdf"), height = 7,
    width = 9)
comb_index_scaled
comb_index_log
comb_index
dev.off()