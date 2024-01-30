### Juvenile all species fit - season specific 
## Use all_species_fit_season as template to fit species-specific models that
## model seasonal abundance as MVN random walk
## 1) Fit spatiotemporal model to each species
## 2) Calculate index for summer and fall (assuming surface and day tow)
## 3) Compare indices including/excluding survey effects
## 4) Calculate fixed effects for spatial covariates
## 5) Calculate spatiotemporal effects (maps by year)
## Aug 22, 223


# ensure mvrfrw branch installed
devtools::install_github("https://github.com/pbs-assess/sdmTMB",
                         ref = "mvrfrw")


library(tidyverse)
library(sdmTMB)
library(ggplot2)
library(sdmTMBextra)

# make subdirectories for storage
dir.create("data/fits", recursive = TRUE, showWarnings = FALSE)
dir.create("data/preds", recursive = TRUE, showWarnings = FALSE)

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
    day_night = as.factor(day_night),
    cpue = n_juv / (volume_km3 * 500)) %>% 
  filter(!species == "steelhead") %>% 
  droplevels() %>% 
  left_join(., yr_key, by = "year")


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

dat_tbl <- dat %>%
  group_by(species) %>%
  group_nest()


# mvrfrw only
fits_list_nb2 <- furrr::future_map(
  dat_tbl$data,
  function(dat_in) {
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
)

# mvrfrw only + FE by year or cycle line
fits_list_nb2_fe <- furrr::future_map2(
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
          scale_depth + year_f,
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
        silent = FALSE
      )
    }
  }
)

dat_tbl_mvrfrw <- dat_tbl %>%
  mutate(fit = fits_list_nb2)
saveRDS(
  dat_tbl_mvrfrw,
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_only.rds")
)

dat_tbl_fe <- dat_tbl %>%
  mutate(fit = fits_list_nb2_fe)
saveRDS(
  dat_tbl_fe, here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_fe.rds")
)

dat_tbl_mvrfrw <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_only.rds")
)
dat_tbl_fe <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_fe.rds")
)

# # make AIC table
dat_tbl <- rbind(
  dat_tbl_mvrfrw %>% mutate(model = "mvrfrw"),
  dat_tbl_fe %>% mutate(model = "mvrfrw_fe")
) %>%
  mutate(
    AIC = purrr::map(.$fit, AIC) %>%
      unlist(),
    log_likelihood = purrr::map(
      .$fit,
      ~ logLik(.x) %>%
        # attr(., "df") %>%
        as.numeric()
    ) %>%
      unlist(),
    n_FE = purrr::map(
      .$fit,
      ~ tidy(.x, effects = "fixed", conf.int = T) %>%
        nrow()
    ) %>%
      unlist()
  )

dat_tbl %>%
  select(species, model, AIC, log_likelihood, n_FE) %>%
  arrange(species, AIC)


# use FEs for species with cycles, otherwise don't
dat_tbl <- dat_tbl %>%
  filter((species %in% c("pink", "sockeye") & model == "mvrfrw_fe") |
           (species %in% c("coho", "chum", "chinook") & model == "mvrfrw"))
saveRDS(
  dat_tbl,
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_final.rds")
)


## CHECKS ----------------------------------------------------------------------

purrr::map(dat_tbl$fit, sanity)

# check residuals
dat_tbl$sims <- purrr::map2(
  dat_tbl$fit, dat_tbl$data, 
  ~ simulate(.x, newdata = .y, nsim = 50)
)
fixed_pred_list <- purrr::map2(
  dat_tbl$fit, dat_tbl$data,
  ~ predict(.x, newdata = .y)$est_non_rf %>% 
    .x$family$linkinv(.)
)
qq_list <- purrr::pmap(
  list(dat_tbl$data, dat_tbl$sims, fixed_pred_list), 
  function(x, y, z) {
    dharma_res <- DHARMa::createDHARMa(
      simulatedResponse = y,
      observedResponse = x$n_juv,
      fittedPredictedResponse = z
    )
    plot(dharma_res)
  }
)

# proportion zeros
purrr::map2(dat_tbl$data, dat_tbl$sims, function(x, y) {
  (sum(x$n_juv == 0) / length(x$n_juv)) / (sum(y == 0) / length(y))
})

# residual spatial dist
resid_list_sum <- purrr::pmap(
  list(dat_tbl$data, dat_tbl$fit, dat_tbl$species),
  function (x, y, z) {
    x$resid <- resid(y)
    ggplot(x %>% filter(season_f == "su")) +
      geom_point(aes(x = utm_x, y = utm_y, color = resid)) +
      facet_wrap(~year) +
      scale_colour_gradient2() +
      ggsidekick::theme_sleek() +
      labs(title = paste(z, unique(x$season_f), sep = "_")) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "top")
  }
)
resid_list_fall <- purrr::pmap(
  list(dat_tbl$data, dat_tbl$fit, dat_tbl$species),
  function (x, y, z) {
    x$resid <- resid(y)
    ggplot(x %>% filter(season_f == "wi")) +
      geom_point(aes(x = utm_x, y = utm_y, color = resid)) +
      facet_wrap(~year) +
      scale_colour_gradient2() +
      ggsidekick::theme_sleek() +
      labs(title = paste(z, unique(x$season_f), sep = "_")) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "top")
  }
)

pdf(here::here("figs", "diagnostics", "spatial_resids.pdf"))
resid_list_sum
resid_list_fall
dev.off()


## PARAMETER ESTIMATES ---------------------------------------------------------

dat_tbl <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_final.rds")
)

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
y_grob <- grid::textGrob(
  "Parameter Estimate", rot = 90,
  gp = grid::gpar(col = "grey30", fontsize = 12, fontface = "plain")
)

png(here::here("figs", "ms_figs_season_mvrw", "fix_ints.png"), 
    height = 4, width = 6,  units = "in", res = 250)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(plot_g, left = y_grob)
)
dev.off()


## random parameter estimates
ran_pars <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species),
  function(x, y) {
    pars1 <- tidy(x, effects = "ran_par", conf.int = T) 
    
    # pull upsilon estimate separately (not currently generated by predict)
    est <- as.list(x$sd_report, "Estimate", report = TRUE)
    se <- as.list(x$sd_report, "Std. Error", report = TRUE)
    upsilon <- data.frame(
      term = "sigma_upsilon",
      log_est = est$log_sigma_U,
      log_se = se$log_sigma_U
    ) %>% 
      mutate(
        estimate = exp(log_est),
        std.error = exp(log_se),
        conf.low = exp(log_est + (stats::qnorm(0.025) * log_se)),
        conf.high = exp(log_est + (stats::qnorm(0.975) * log_se)) 
      ) %>% 
      select(term, estimate, std.error, conf.low, conf.high)
    
    rbind(pars1, upsilon) %>% 
      mutate(
        species = y
      )
  }
) %>% 
  bind_rows()
ran_pars$term <- fct_recode(ran_pars$term, "sigma_omega" = "sigma_Z", 
                            "kappa" = "range")

png(here::here("figs", "ms_figs_season_mvrw", "ran_pars.png"), height = 4,
    width = 8,
    units = "in", res = 250)
ggplot(
  ran_pars,
  aes(species, estimate, ymin = conf.low, ymax = conf.high,
      fill = species)
) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = col_pal) +
  ylab("Parameter Estimate") +
  geom_pointrange(position = position_dodge(width = 0.4), shape = 21) +
  facet_wrap(~term, scales = "free_y") +
  guides(fill = "none") +
  theme(axis.title.x = element_blank())
dev.off()


# headrope depth
target_dat <- expand.grid(
  day_night = unique(dat$day_night),
  season_f = unique(dat$season_f),
  # sequence from surface to reasonable maximum headrope depth
  headrope_depth = seq(0, 60, length.out = 100),
  scale_dist = 0,
  survey_f = unique(dat$survey_f),
  sox_cycle = "1",
  pink_cycle = "1",
  year = unique(dat$year)
) %>% 
  mutate(
    scale_depth = (headrope_depth - mean(dat$target_depth)) / 
      sd(dat$target_depth),
    year_f = as.factor(year),
    season_year_f = paste(season_f, year_f, sep = "_") %>% as.factor()
  ) %>% 
  filter(
    day_night == "DAY", survey_f == "hss", year_f == "2012"
  ) %>% 
  distinct()

depth_preds <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species),
  function(x, y, z) {
    pp <- predict(x, newdata = target_dat, se_fit = FALSE, re_form = NA, 
            re_form_iid = NA, nsim = 100) 
    target_dat %>% 
      mutate(
        est = apply(pp, 1, mean),
        est_se = apply(pp, 1, sd),
        species = y
        )
  }
) %>%
  bind_rows() 

depth_preds2 <- depth_preds %>% 
  filter(season_f == "su") %>% 
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
  theme(axis.title.y = element_blank(),
        panel.spacing = unit(0.8, "lines"),
        axis.title.x = element_text(size = 12)) +
  ylab("Scaled Abundance Index") 


# distance to coast 
dist_dat <- expand.grid(
  dist_to_coast_km = seq(0.1, max(dat$dist_to_coast_km), length.out = 100),
  scale_depth = 0,
  day_night = unique(dat$day_night),
  season_f = unique(dat$season_f),
  target_depth = 0,
  sox_cycle = "1",
  pink_cycle = "1",
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
    day_night == "DAY", survey_f == "hss", year_f == "2012"
  ) 

dist_preds <- purrr::pmap(
  list(dat_tbl$fit, dat_tbl$species),
  function(x, y, z) {
    pp <- predict(x, newdata = dist_dat, se_fit = FALSE, re_form = NA, 
                  re_form_iid = NA, nsim = 100) 
    dist_dat %>% 
      mutate(
        est = apply(pp, 1, mean),
        est_se = apply(pp, 1, sd),
        species = y
      )
  }
) %>%
  bind_rows()

dist_preds2 <- dist_preds %>%
  filter(season_f == "su") %>% 
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
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)) +
  guides(fill = "none")

plot_g2 <- cowplot::plot_grid(depth_plot, dist_plot, ncol = 1)
y_grob <- grid::textGrob(
  "Scaled Abundance Index", rot = 90,
  gp = grid::gpar(col = "grey30", fontsize = 14, fontface = "plain")
  )

png(here::here("figs", "ms_figs_season_mvrw", "smooth_effs.png"), height = 5,
    width = 8,    units = "in", res = 250)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(plot_g2, left = y_grob)
)
dev.off()



## SPATIAL GRID  ---------------------------------------------------------------

# scalar for spatial predictions; since preds are in m3, multiply by
# (1 * 1 * 13) because effort in m but using 1x1 km grid cells and 
# assuming mean net opening (13 m); multiply by 500 to match input effort that
# was scaled to be centered near 1 unit for offsetting purposes
sp_scalar <- 1 * (13 / 1000) * 500


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

year_season_key <- expand.grid(
  year = unique(dat_in$year) %>% sort(),
  season_f = unique(dat_in$season_f) %>% sort()
) %>% 
  arrange(year, season_f) 

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

#add unique years and seasons
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
  select(-c(depth, slope, shore_dist)) %>% 
  left_join(., yr_key, by = "year")


# predictions for index (set to HSS reference values and constrain in summer)
index_grid_hss <- index_grid %>% 
  filter(survey_f == "hss") %>% 
  mutate(day_night = as.factor(day_night),
         trim = ifelse(
           season_f == "wi" & utm_y_1000 < 5551, "yes", "no"
         )) %>%
  #subset to northern domain that's well sampled in recent years
  filter(trim == "no")

saveRDS(index_grid_hss, here::here("data", "index_hss_grid.rds"))


# plot grid domains
index_grid_plot <- index_grid_hss %>% 
  filter(!season_f == "sp",
         year == "2020") %>% 
  mutate(season_f2 = factor(season_f, labels = c("summer", "fall")))
index_grids <- ggplot() + 
  geom_raster(
    data = index_grid_plot,
    aes(X, Y, fill = season_f2)
  ) +
  coord_fixed() +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_brewer(
    type = "qual", palette = 5
  ) +
  facet_wrap(~season_f2) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))

png(here::here("figs", "ms_figs_season_mvrw", "index_grids.png"),
    height = 4.5, width = 8, units = "in", res = 200)
index_grids
dev.off()

 
ind_preds_summer <- purrr::map(
  dat_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid_hss %>%
              filter(season_f == "su"),
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)
ind_preds_fall <- purrr::map(
  dat_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid_hss %>%
              filter(season_f == "wi"),
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)

index_list_summer <- furrr::future_map(
  ind_preds_summer,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
index_list_summer_out <- purrr::map2(
  index_list_summer, dat_tbl$species,
  function (x, y) {
    x %>%
      mutate(season_f = "su",
             species = y)
  }
)


index_list_fall <- furrr::future_map(
  ind_preds_fall,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
index_list_fall_out <- purrr::map2(
  index_list_fall, dat_tbl$species,
  function (x, y) {
    x %>%
      mutate(season_f = "wi",
             species = y)
  }
)

saveRDS(c(index_list_summer_out, index_list_fall_out),
        here::here("data", "season_index_list_mvrfrw.rds"))


index_list <- readRDS(
  here::here("data", "season_index_list_mvrfrw.rds")
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


index_dat <- index_list %>%
  bind_rows() %>%
  mutate(
    season = fct_recode(season_f, "summer" = "su", "fall" = "wi"),
    survey = case_when(
      season == "fall" & year %in% fall_years ~ "sampled",
      season == "summer" & year %in% summer_years ~ "sampled",
      TRUE ~ "no survey"
    ),
    # identify surveys with relatively fewer samples
    sparse = case_when(
      survey == "no survey" ~ "yes",
      season == "summer" & 
        year %in% c("1998", "2000", "2003", "2013", "2014") ~ "yes",
      season == "fall" & 
        year %in% c("2016", "2017", "2019") ~ "yes",
      TRUE ~ "no"
    )
    ) %>% 
  group_by(season, species) %>% 
  mutate(mean_est = mean(est),
         geo_mean_est = exp(mean(log_est))) %>% 
  ungroup() %>% 
  mutate(
    scale_geo_est = est / geo_mean_est,
    scale_geo_lwr = lwr / geo_mean_est,
    scale_geo_upr = upr / geo_mean_est
    # scale_geo_lwr = ifelse(survey == "sampled", lwr / geo_mean_est, scale_geo_est),
    # scale_geo_upr = ifelse(survey == "sampled", upr / geo_mean_est, scale_geo_est)
  )


## index plots
shape_pal <- c(21, 1)
names(shape_pal) <- unique(index_dat$survey)

# log abundance 
log_index <- ggplot(index_dat, aes(year, est)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr, fill = species,
                      shape = survey)) +
  geom_hline(aes(yintercept = mean_est), lty = 2) +
  labs(x = "Year", y = "Abundance Index") +
  scale_shape_manual(values = shape_pal) +
  scale_y_log10() +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13))

png(here::here("figs", "ms_figs_season_mvrw", "log_index.png"),
    height = 8, width = 8, units = "in", res = 250)
log_index
dev.off()


# scaled geometric mean abundance 
geo_mean_index <- ggplot(index_dat, aes(year, scale_geo_est)) +
  geom_pointrange(aes(ymin = scale_geo_lwr, ymax = scale_geo_upr, 
                      fill = species,
                      shape = survey)) +
  labs(x = "Year", y = "Abundance Index (Geometric Mean)") +
  scale_shape_manual(values = shape_pal) +
  scale_y_log10(labels = scales::comma) +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)) 

png(here::here("figs", "ms_figs_season_mvrw", "gm_index.png"), 
    height = 8, width = 8, units = "in", res = 250)
geo_mean_index
dev.off()


# log abundance coded by survey coverage (looks good)
# shape_pal2 <- c(1, 21)
# names(shape_pal2) <- unique(index_dat$sparse)
# log_index_sparse <- ggplot(index_dat, aes(year, est)) +
#   geom_pointrange(aes(ymin = lwr, ymax = upr, fill = species,
#                       shape = sparse)) +
#   scale_y_log10() +
#   # geom_hline(aes(yintercept = mean_log_est), lty = 2) +
#   labs(x = "Year", y = "Log Abundance Index") +
#   scale_shape_manual(values = shape_pal2) +
#   ggsidekick::theme_sleek() +
#   facet_grid(species~season, scales = "free_y") +
#   scale_fill_manual(values = col_pal) +
#   theme(axis.title.x = element_blank())
# 
# png(here::here("figs", "diagnostics", "log_index_sparse.png"),
#     height = 8, width = 8, units = "in", res = 200)
# log_index_sparse
# dev.off()


## scaled abundance 
scaled_index <- index_dat %>% 
  # filter(survey == "sampled")%>%
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
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13))

png(here::here("figs", "ms_figs_season_mvrw", "scaled_index.png"), 
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


## SURVEY EFFECTS --------------------------------------------------------------


## fit alternative models without survey or day/night effects
naive_fit <- furrr::future_map2(
  dat_tbl$data, dat_tbl$species,
  function(dat_in, sp) {
    if (sp == "pink") {
      sdmTMB(
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
        silent = FALSE
      )
    } else if (sp == "sockeye") {
      sdmTMB(
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
        silent = FALSE
      )
    } else {
      sdmTMB(
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
        silent = FALSE
      )
    }
  }
)

dat_tbl_naive <- dat_tbl %>%
  mutate(fit = naive_fit)
saveRDS(
  dat_tbl_naive,
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_naive.rds")
)

dat_tbl_naive <- readRDS(
  here::here("data", "fits", "all_spatial_varying_nb2_mvrfrw_naive.rds")
)


# AIC table
survey_eff_aic <- purrr::map(
  dat_tbl$fit,
  AIC
) %>% 
  unlist()
no_survey_eff_aic <- purrr::map(
  dat_tbl_naive$fit,
  AIC
) %>% 
  unlist()

survey_eff_aic - no_survey_eff_aic

both_fits <- c(dat_tbl$fit, dat_tbl_naive$fit)
survey_eff_out <- data.frame(
  species = rep(dat_tbl$species, times = 2),
  model = rep(c("survey effects", "no survey effects"), each = 5),
  AIC = purrr::map(both_fits, AIC) %>%
    unlist(),
  log_likelihood = purrr::map(
    both_fits,
    ~ logLik(.x) %>%
      as.numeric()
  ) %>%
    unlist(),
  n_FE = purrr::map(
    both_fits,
    ~ tidy(.x, effects = "fixed", conf.int = T) %>%
      nrow()
  ) %>%
    unlist()
) %>% 
  arrange(
    species, model
  )
write.csv(survey_eff_out,
          here::here("data", "survey_eff_model_comp.csv"),
          row.names = FALSE)
  

# compare HSS index to one where survey design is not accounted for
ind_preds_naive <- purrr::map(
  dat_tbl_naive %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid_hss %>%
              filter(season_f == "su"),
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)

index_list_naive <- furrr::future_map(
  ind_preds_naive,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
index_list_naive_out <- purrr::map2(
  index_list_naive, dat_tbl_naive$species,
  function (x, y) {
    x %>%
      mutate(season_f = "su",
             species = y)
  }
)

saveRDS(index_list_naive_out,
        here::here("data", "fits", "season_index_naive_survey_list_mvrfrw.rds"))
index_list_naive_out <- readRDS(
  here::here("data", "fits", "season_index_naive_survey_list_mvrfrw.rds"))

index_dat2 <- index_list_naive_out %>% 
  bind_rows() %>% 
  mutate(model = "no survey effects") %>% 
  filter(season_f == "su")

index_by_model <- index_dat %>% 
  mutate(
    model = "survey effects"
  ) %>% 
  filter(season_f == "su") %>% 
  select(colnames(index_dat2)) %>% 
  bind_rows(., index_dat2)

shape_pal2 <- c(24, 21)
names(shape_pal2) <-  c("no survey effects", "survey effects")

png(here::here("figs", "ms_figs_season_mvrw", "log_index_both_surveys.png"), 
    height = 8, width = 8, units = "in", res = 250)
ggplot(index_by_model, aes(year, est)) +
  geom_pointrange(
    aes(ymin = lwr, ymax = upr, fill = species, shape = model),
    position = position_dodge(width = .7)
  ) +
  labs(x = "Year", y = "Abundance Index") +
  scale_shape_manual(values = shape_pal2) +
  scale_y_log10() +
  ggsidekick::theme_sleek() +
  facet_wrap(~species, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = col_pal, guide = "none") +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13))
dev.off()

# sd of index in summer vs fall
index_dat %>% 
  group_by(species, season_f) %>% 
  summarize(sd_index = sd(log_est)) %>% 
  pivot_wider(names_from = season_f, values_from = sd_index) %>% 
  mutate(
    ppn = wi / su)


# ESTIMATE DECLINE -------------------------------------------------------------

index_lm_dat <- rbind(
  index_dat %>% 
    filter(season_f == "su",
           year > 2010) %>% 
    dplyr::select(log_est, se, species, year, season_f) %>% 
    mutate(survey_eff = "yes",
           weight = (1 / se^2)),
  index_dat2 %>% 
    dplyr::select(log_est, se, species, year, season_f) %>% 
    filter(year > 2010) %>% 
    mutate(survey_eff = "no",
           weight = (1 / se^2))
) %>% 
  group_by(species, survey_eff) %>% 
  group_nest()

# fit
index_lm_dat$fit <- purrr::map(
  index_lm_dat$data, 
  ~ lm(log_est ~ year, data = .x, weights = weight)
)

# predict
index_lm_dat$pred <- purrr::map2(
  index_lm_dat$fit, index_lm_dat$data, 
  function(x , y) {
    pp <- predict(x, se.fit = TRUE) 
    y %>% 
      mutate(
        pred_mean = as.numeric(pp$fit),
        pred_se = as.numeric(pp$se.fit)
      )
  }
)

# plot trends
index_lm_dat %>% 
  dplyr::select(species, survey_eff, pred) %>% 
  unnest(cols = c(pred)) %>%
  mutate(upr = pred_mean + 1.96 * pred_se,
         lwr = pred_mean - 1.96 * pred_se) %>% 
  ggplot(., aes(x = year, y = pred_mean, colour = species)) +
  geom_line(aes(linetype = survey_eff)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species, 
                  linetype = survey_eff),
              alpha = 0.2) +
  facet_wrap(~species) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ggsidekick::theme_sleek() 


# slope estimates 
index_lm_dat$coefs <- purrr::map(
   index_lm_dat$fit,
   ~ tidy(.x, conf.int = TRUE)
 ) 
 
index_lm_plot_dat <- index_lm_dat %>% 
  dplyr::select(species, survey_eff, coefs) %>% 
  unnest(cols = c(coefs)) %>% 
  filter(term == "year") %>% 
  mutate(
    survey_eff = fct_relevel(as.factor(survey_eff), "no", after = Inf)
  )
  
png(here::here("figs", "ms_figs_season_mvrw", "decline_survey_effect.png"), 
    height = 3, width = 8, units = "in", res = 200)
ggplot(index_lm_plot_dat, 
       aes(x = survey_eff, y = (exp(estimate) - 1) * 100)) +
  geom_pointrange(
    aes(ymin = (exp(conf.low) - 1) * 100, ymax = (exp(conf.high) - 1) * 100,
        fill = species),
    shape = 21
  ) +
  geom_hline(aes(yintercept = 0), linetype = 2, colour = "red") +
  facet_wrap(~species, nrow = 1) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  ggsidekick::theme_sleek() +
  labs(x = "Survey/Diel Effects Accounted For", y = "Estimated Decline") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)) +
  scale_y_continuous(
    breaks = seq(-60, 20, by = 20),
    labels = paste(seq(-60, 20, by = 20), "%", sep = "")
  )
dev.off()

 
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
        ys_index = x$ys_index,
        target_depth = 0,
        day_night = "DAY",
        scale_depth = (target_depth - mean(dat$target_depth)) /
          sd(dat$target_depth),
        scale_dist = (dist_to_coast_km  - mean(dat$dist_to_coast_km )) /
          sd(dat$dist_to_coast_km )
      )
  }) %>%
  bind_rows() %>%
  filter(
    !survey_f == "ipes"
  ) %>%
  dplyr::select(-c(depth, slope)) %>% 
  left_join(
    ., yr_key, by = "year"
  )

spatial_preds_list <- furrr::future_map2(
  dat_tbl$fit,
  dat_tbl$species,
  ~ {
    # separate predictions to generate upsilon estimates
    pp <- predict(.x,
                  newdata = spatial_grid,
                  return_tmb_report = TRUE)
    p_out <- predict(.x,
                     newdata = spatial_grid,
                     se_fit = FALSE, re_form = NULL)

    p_out %>%
      mutate(species = .y,
             upsilon_stc = as.numeric(pp$proj_upsilon_st_A_vec))
  }
)

spatial_preds <- spatial_preds_list %>%
  bind_rows() %>%
  filter(!season_f == "sp")
saveRDS(spatial_preds,
        here::here("data", "fits", "sp_preds_mvrfrw.rds"))

spatial_preds <- readRDS(here::here("data", "fits", "sp_preds_mvrfrw.rds"))

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

png(here::here("figs", "ms_figs_season_mvrw", "scaled_sp_preds_su.png"), 
    height = 8.5, width = 8.5, units = "in", res = 200)
ggplot() + 
  geom_raster(
    data = sub_spatial %>% 
      filter(season_f == "su",
        # species == "coho",     
        year %in% year_seq),
    aes(X, Y, fill = scale_est2)
  ) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Scaled\nAbundance"
  ) +
  facet_grid(year~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.key.size = unit(1.1, 'cm')) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))
dev.off()


png(here::here("figs", "ms_figs_season_mvrw", "scaled_sp_preds_fall.png"), 
    height = 8.5, width = 8.5, units = "in", res = 200)
ggplot() + 
  geom_raster(
    data = sub_spatial %>% 
      filter(season_f == "wi",
             year %in% year_seq),
    aes(X, Y, fill = scale_est2)
  ) +
  coord_fixed() +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_c(
    trans = "sqrt",
    name = "Scaled\nAbundance"
  ) +
  facet_grid(year~species) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.key.size = unit(1.1, 'cm')) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))
dev.off()


# all yearly preds for each species
full_spatial <- spatial_preds %>% 
  group_by(species) %>% 
  mutate(
    grid_est = sp_scalar * exp(est),
    scale_est = grid_est / max(grid_est)
  ) %>% 
  ungroup() 
scale_max <- quantile(full_spatial$scale_est, 0.99)
full_spatial$scale_est2 <- ifelse(
  full_spatial$scale_est > scale_max, 
  scale_max, 
  full_spatial$scale_est
)
summer_list <- purrr::map(
  unique(dat$species) %>% sort(),
  ~ ggplot() +
      geom_raster(
        data = full_spatial %>%
          filter(season_f == "su", species == .x),
        aes(X, Y, fill = scale_est2)
      ) +
      coord_fixed() +
      geom_sf(data = coast, color = "black", fill = "white") +
      ggsidekick::theme_sleek() +
      scale_fill_viridis_c(
        trans = "sqrt",
        name = "Scaled\nAbundance"
      ) + 
      facet_wrap(~year) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "top",
            legend.key.size = unit(1.1, 'cm')) +
    labs(title = .x)
)
pdf(here::here("figs", "ms_figs_season_mvrw", "all_year_preds",
               "summer_year_rf.pdf"))
summer_list
dev.off()

fall_list <- purrr::map(
  unique(dat$species) %>% sort(),
  ~ ggplot() +
    geom_raster(data = full_spatial %>%
                  filter(season_f == "wi", species == .x),
                aes(X, Y, fill = scale_est2)) +
    coord_fixed() +
    geom_sf(data = coast, color = "black", fill = "white") +
    ggsidekick::theme_sleek() +
    scale_fill_viridis_c(
      trans = "sqrt",
      name = "Scaled\nAbundance"
    ) + 
    facet_wrap(~year) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position = "top",
          legend.key.size = unit(1.1, 'cm')) +
    labs(title = .x)
)
pdf(here::here("figs", "ms_figs_season_mvrw", "all_year_preds",
               "fall_year_rf.pdf"))
fall_list
dev.off()



png(here::here("figs", "ms_figs_season_mvrw", "su_year_rf.png"), 
    height = 8.5, width = 8.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = sub_spatial %>% 
                filter(season_f == "su",
                       year %in% year_seq),
              aes(X, Y, fill = upsilon_stc)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(year~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1.1, 'cm')
  ) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))
dev.off()


png(here::here("figs", "ms_figs_season_mvrw", "fa_year_rf.png"), 
    height = 8.5, width = 8.5, units = "in", res = 200)
ggplot() +
  geom_raster(data = sub_spatial %>% 
                filter(season_f == "wi",
                       year %in% year_seq),
              aes(X, Y, fill = upsilon_stc)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(year~species) +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    axis.text = element_blank(),
    legend.key.size = unit(1.1, 'cm')
  ) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))
dev.off()


# spatial random effects by season
omega_season <- sub_spatial %>% 
  # values will be duplicated across years, seasons and surveys; select one
  filter(year_f == year_seq[1], season_f == "wi") %>% 
  dplyr::select(X, Y, species, starts_with("zeta_s_season_")) %>%
  pivot_longer(cols = starts_with("zeta"), values_to = "omega_est", 
               names_to = "season", names_prefix = "zeta_s_season_f") %>% 
  filter(!season == "sp") %>% 
  mutate(season = fct_recode(season, "summer" = "su", "fall" = "wi"))

png(here::here("figs", "ms_figs_season_mvrw", "seasonal_rf.png"), 
    height = 7, width = 6, units = "in", res = 200)
ggplot() +
  geom_raster(data = omega_season, aes(X, Y, fill = omega_est)) +
  geom_sf(data = coast, color = "black", fill = "white") +
  ggsidekick::theme_sleek() +
  scale_fill_gradient2(name = "Spatial\nField") +
  facet_grid(species~season) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(.95, 'cm')
  ) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0))
dev.off()


# alternative formulations - average abundance
mean_spatial <- spatial_preds %>% 
  mutate(location = paste(X, Y, sep = "_"),
         grid_est = sp_scalar * exp(est)) %>% 
  group_by(X, Y, species, season_f, location) %>% 
  summarize(
    mean_grid_est = mean(grid_est),
    median_grid_est = median(grid_est)
  ) %>% 
  group_by(species, season_f) %>% 
  mutate(
    scale_mean_grid_est = mean_grid_est / max(mean_grid_est),
    scale_median_grid_est = median_grid_est / max(median_grid_est)
  ) %>% 
  ungroup() 
scale_max <- quantile(mean_spatial$scale_median_grid_est, 0.99)
mean_spatial$scale_est2 <- ifelse(
  mean_spatial$scale_median_grid_est > scale_max, 
  scale_max, 
  mean_spatial$scale_median_grid_est
)

inset_rect <- data.frame(
  x1 = 645000, x2 = 720000, y1 = 5465000, y2 = 5545000
)

plot_mean_foo <- function (x) {
  ggplot() + 
    geom_raster(
      data = x,
      aes(X, Y, fill = scale_est2)
    ) +
    coord_fixed() +
    ggsidekick::theme_sleek() +
    scale_fill_viridis_c(
      trans = "sqrt",
      name = "Scaled\nAbundance"
    )  +
    facet_wrap(~species, ncol = 1) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.key.size = unit(1.1, 'cm')
    ) 
} 

summer_main <- mean_spatial %>% 
  filter(season_f == "su") %>% 
  plot_mean_foo(x = .) +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_rect(data = inset_rect,
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = NA, lty = 2, col = "red", linewidth = 0.7) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0)) 
summer_inset <-  mean_spatial %>% 
  filter(season_f == "su") %>% 
  plot_mean_foo(x = .) +
  scale_y_continuous(limits = c(5465000, 5535500), expand = c(0, 0)) +
  scale_x_continuous(limits = c(645000, 720000), expand = c(0, 0))
fall_main <- mean_spatial %>% 
  filter(season_f == "wi") %>% 
  plot_mean_foo(x = .) +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_rect(data = inset_rect,
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = NA, lty = 2, col = "red", linewidth = 0.7) +
  scale_x_continuous(limits = c(462700, 814800), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5350050, 5681850), expand = c(0, 0)) 
fall_inset <-  mean_spatial %>% 
  filter(season_f == "wi") %>% 
  plot_mean_foo(x = .) +
  scale_y_continuous(limits = c(5465000, 5535500), expand = c(0, 0)) +
  scale_x_continuous(limits = c(645000, 720000), expand = c(0, 0))

plot_legend <- cowplot::get_legend(summer_main)
p_cols <- cowplot::plot_grid(
  summer_main + theme(legend.position = "none"),
  summer_inset + theme(legend.position = "none"),
  fall_main + theme(legend.position = "none"),
  fall_inset + theme(legend.position = "none"),
  nrow = 1
)

png(here::here("figs", "ms_figs_season_mvrw", "mean_spatial_preds_season_scaled.png"), 
    height = 7, width = 7, units = "in", res = 200)
cowplot::plot_grid(
  p_cols, plot_legend, ncol = 2, rel_widths = c(1, 0.125)
  # plot_legend, p_cols, nrow = 2, rel_heights = c(0.1, 1)
)
dev.off()


## SUMMARY TABLE INFO ----------------------------------------------------------

#ntows
length(dat$unique_event)

#ntows w/ zero salmon
dat %>% 
  filter(n_juv == "0") %>% 
  group_by(species) %>% 
  tally() %>% 
  mutate(
    ppn = n / length(dat$unique_event)
  )

# effort
dum <- dat %>% filter(species == "chinook", n_juv < 200)
mean(dum$volume_km3)
sd(dum$volume_km3)

# depth
mean(dum$target_depth)
sd(dum$target_depth)


## IPES ONLY PREDICTIONS -------------------------------------------------------

index_grid_ipes <- index_grid %>% 
  filter(survey_f == "ipes") %>% 
  mutate(day_night = as.factor(day_night),
         trim = ifelse(
           season_f == "wi" & utm_y_1000 < 5551, "yes", "no"
         )) %>%
  #subset to northern domain that's well sampled in recent years
  filter(trim == "no")

ind_preds_summer_ipes <- purrr::map(
  dat_tbl %>% pull(fit),
  ~ {
    predict(.x,
            newdata = index_grid_ipes %>%
              filter(season_f == "su"),
            se_fit = FALSE, re_form = NULL, return_tmb_object = TRUE)
  }
)

index_list_summer_ipes <- furrr::future_map(
  ind_preds_summer_ipes,
  get_index,
  area = sp_scalar,
  bias_correct = TRUE
)
index_list_summer_ipes_out <- purrr::map2(
  index_list_summer_ipes, dat_tbl$species,
  function (x, y) {
    x %>%
      mutate(season_f = "su",
             species = y)
  }
)

saveRDS(index_list_summer_ipes_out,
        here::here("data", "season_index_list_mvrfrw_summer_ipes.rds"))
