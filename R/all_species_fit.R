### Juvenile all species fit 
## Use chinook_fit as a template to fit different species models
## 1) Fit spatial model for each species to ensure reasonable convergence
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
    effort = log(volume_m3),
    week = lubridate::week(date),
    vessel = as.factor(vessel),
    # redefine bins for all species except chinook
    target_depth_bin = case_when(
      species %in% c("SOCKEYE", "PINK") & 
        target_depth_bin %in% c("30", "45", "60") ~ ">30",
      species %in% c("COHO", "CHUM") & 
        target_depth_bin %in% c("45", "60") ~ ">45",
      species == "CHINOOK" & target_depth_bin == "60" ~ ">60",
      TRUE ~ as.character(target_depth_bin)
    )
  ) %>% 
  droplevels() 


# keep early years in dataset for now but exclude 2022 and consider dropping 
# years with v. limited summer data
dat_in <- dat %>% 
  filter(!year == "2022",
         !bath_depth_mean_m < 0)


coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -137, ymin = 47, xmax = -121.25, ymax = 57) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# make meshes (differ among groups because off pooling of depth strata)
dat_tbl <- dat_in %>% 
  group_by(species) %>% 
  group_nest() %>% 
  mutate(
    spde = purrr::map(data, make_mesh,  c("utm_x_1000", "utm_y_1000"), 
                      cutoff = 15, type = "kmeans"),
    bspde = purrr::map(spde, add_barrier_mesh,
                       coast_utm, range_fraction = 0.1,
                       # scaling = 1000 since UTMs were rescaled above
                       proj_scaling = 1000)
  )


# prep multisession
ncores <- parallel::detectCores() 
future::plan(future::multisession, workers = ncores - 3)


### EXP PLOTS ------------------------------------------------------------------

ggplot(dat_in) +
  geom_boxplot(aes(x = target_depth_bin, y = log(n_juv))) +
  facet_wrap(~species, scales = "free_x")

ggplot(dat_in) +
  geom_point(aes(x = bath_depth_mean_m, y = log(n_juv))) +
  facet_wrap(~species)

ggplot(dat_in) +
  geom_point(aes(x = dist_to_coast_km, y = log(n_juv))) +
  facet_wrap(~species)

ggplot(dat_in) +
  geom_point(aes(x = week, y = log(n_juv))) +
  facet_wrap(~species)

ggplot(dat_in) +
  geom_boxplot(aes(x = day_night, y = log(n_juv))) +
  facet_wrap(~species)


### FIT SPATIAL ----------------------------------------------------------------


spatial_mod <- furrr::future_map2(
  dat_tbl$data, dat_tbl$bspde, function (x, spde_in) {
    sdmTMB(
      n_juv ~ 1 +  
        dist_to_coast_km +
        # s(dist_to_coast_km, bs = "tp", k = 3) +
        s(week, bs = "cc", k = 4) +
        day_night +
        target_depth_bin +
        survey_f,
      offset = x$effort,
      data = x,
      mesh = spde_in,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      # spatiotemporal = "ar1",
      anisotropy = FALSE,
      priors = sdmTMBpriors(
        matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
      ),
      control = sdmTMBcontrol(
        nlminb_loops = 2,
        newton_loops = 1
      )
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)

purrr::map(spatial_mod, sanity)
purrr::map(spatial_mod, summary)

# check single spatial mod
sp_mod <- sdmTMB(
  n_juv ~ 1 +  
    dist_to_coast_km +
    s(week, bs = "cc", k = 4) +
    day_night +
    target_depth_bin +
    survey_f,
  offset = dat_tbl$data[[3]]$effort,
  data = dat_tbl$data[[3]],
  mesh = dat_tbl$bspde[[3]],
  family = sdmTMB::nbinom2(),
  spatial = "on",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
  ),
  control = sdmTMBcontrol(
    nlminb_loops = 2,
    newton_loops = 1
  )
)
spatial_preds <- predict(sp_mod, 
                         exp_grid %>% filter(year == "2017") %>% distinct())


### FIT SATURATED --------------------------------------------------------------

st_mod <- furrr::future_pmap(
  list(dat_tbl$data, dat_tbl$bspde, dat_tbl$species),
  function (x, spde_in, sp_in) {
    sdmTMB(
      n_juv ~ 1 +
        dist_to_coast_km +
        s(week, bs = "tp", k = 5) +
        s(target_depth, bs = "tp", k = 4) +
        day_night +
        survey_f,
      offset = x$effort,
      data = x,
      mesh = spde_in,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "AR1",
      time = "year",
      anisotropy = FALSE,
      priors = sdmTMBpriors(
        matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
      ),
      control = sdmTMBcontrol(
        nlminb_loops = 2,
        newton_loops = 1
      )
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)


purrr::map(st_mod, sanity)
purrr::map(st_mod, summary)
purrr::map(dat_tbl$data, function (x) range(x$n_juv))


saveRDS(dat_tbl, here::here("data", "fits", "st_mod_all_sp.rds"))


dat_tbl <- readRDS(here::here("data", "fits", "st_mod_all_sp.rds"))


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


## MAKE FE PREDICTIONS ---------------------------------------------------------


# make conditional predictive dataframes
week_dat <- data.frame(
 week = seq(5, 45, length.out = 100),
 dist_to_coast_km = median(dat_in$dist_to_coast_km),
 day_night = "DAY",
 target_depth = 0,
 survey_f = "hss",
 year = 2011L
)
dist_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = seq(min(dat_in$dist_to_coast_km), 
                         max(dat_in$dist_to_coast_km), length.out = 100),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year = 2011L
)
day_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = c("DAY", "NIGHT"),
  target_depth = 0,
  survey_f = "hss",
  year = 2011L
)
target_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = "DAY",
  target_depth = seq(min(dat_in$target_depth), 
                         max(dat_in$target_depth), length.out = 100),
  survey_f = "hss",
  year = 2011L
)
survey_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = c("hss", "ipes"),
  year = 2011L
)
pred_tbl <- tibble(
  var = c("week",
          "dist_to_coast_km",
          "day_night",
          "target_depth",
          "survey_f"),
  data = list(week_dat, dist_dat, day_dat, target_dat, survey_dat),
  plot = c("line", "line", "dot", "line", "dot")
)


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

plot_eff <- function (dat, x_var, type = c("line", "dot")) {
  p <- dat %>% 
    mutate(
      exp_est = exp(est),
      up = exp(est + 1.96 * est_se),
      lo = exp(est - 1.96 * est_se)
    ) %>% 
    ggplot(
    ., 
    aes_string(x_var, "exp_est", ymin = "lo", ymax = "up")
  ) +
    facet_wrap(~species, scales = "free_y") +
    labs(title = x_var) +
    ggsidekick::theme_sleek()
  if (type == "line") {
    p +
      geom_line() +
      geom_ribbon(alpha = 0.3)
  }
  if (type == "dot") {
    p + 
      geom_pointrange() 
  }
} 

purrr::pmap(list(pred_tbl$pred_dat, pred_tbl$var, pred_tbl$plot), plot_eff)
plot_eff(dat = pred_tbl$pred_dat[[1]],
         x_var = pred_tbl$var[[1]])

ggplot(pred_list[[1]], 
       aes(week, exp(est),
           ymin = exp(est - 1.96 * est_se), ymax = exp(est + 1.96 * est_se))) +
  geom_line() +
  geom_ribbon(alpha = 0.3) + 
  facet_wrap(~species, scales = "free_y") +
  ggsidekick::theme_sleek()

## MAKE SPATIAL PREDICTIONS ----------------------------------------------------

# spatial distribution of residuals
grid <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  mutate(utm_x_1000 = X / 1000,
         utm_y_1000 = Y / 1000,
         dist_to_coast_km = shore_dist / 1000)

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(dat_in$year),
  survey_f = unique(dat_in$survey_f)
) %>%
  mutate(id = row_number()) %>%
  split(., .$id) %>%
  purrr::map(., function (x) {
    grid %>% 
      mutate(
        year = x$year,
        survey_f = x$survey_f,
        target_depth_bin = "0",
        # summer week 
        week = 25,
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


dat_tbl$preds <- purrr::map(dat_tbl$st_mod, function (x) {
  predict(x, newdata = exp_grid)
})

stock_preds <- dat_tbl %>% 
  select(species, preds) %>% 
  unnest(cols = preds) %>% 
  select(-est, -epsilon_st, -year) %>% 
  distinct()


plot_map <- function(dat, column) {
  ggplot(dat, aes_string("utm_x_1000", "utm_y_1000", fill = column)) +
    geom_raster() +
    coord_fixed() +
    ggsidekick::theme_sleek()
}

fe_plot_list <- purrr::map(dat_tbl$preds, function (x) {
  plot_map(x, "exp(est)") +
    scale_fill_viridis_c(
      trans = "sqrt",
      limits = c(0, quantile(exp(dd$est), 0.995))
    ) +
    facet_wrap(~year) +
    ggtitle("Prediction (fixed effects + all random effects)")
})

plot_map(stock_preds, "est_rf") +
  scale_fill_gradient2(
  ) +
  ggtitle("Prediction (spatial random effects only)") +
  facet_wrap(~species)



# index from summer
summer_years <- dat_in %>% 
  filter(season_f == "su",
         #remove 2021 survey years (too few tows)
         !year == "2021") %>%
  pull(year) %>% 
  unique()


# fix to HSS survey
ind_preds <- purrr::map(dat_tbl$st_mod, function (x) {
  predict(x, newdata = exp_grid %>% filter(survey_f == "hss"), 
          return_tmb_object = TRUE)
})
index_list <- purrr::map(ind_preds, get_index, bias_correct = TRUE)
index_df <- purrr::map2(index_list, dat_tbl$species, function (x, sp) {
  x$species <- sp
  return(x)
}) %>% 
  bind_rows()

index_plot <- ggplot(index_df %>% filter(year %in% summer_years), 
                     aes(year, est)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  labs(x = "Year", y = "Count") +
  ggsidekick::theme_sleek() +
  facet_wrap(~species, scales = "free_y")


# separate survey effects
ind_preds2 <- purrr::map(dat_tbl$st_mod, function (x) {
  predict(x, newdata = exp_grid %>% filter(fake_survey == "0"), 
          return_tmb_object = TRUE)
})
index_list2 <- purrr::map(ind_preds2, get_index, bias_correct = TRUE)
index_df2 <- purrr::map2(index_list2, dat_tbl$species, function (x, sp) {
  x$species <- sp
  return(x)
}) %>% 
  bind_rows()


index_list <- list(survey_eff = index_df, no_survey_eff = index_df2)

index1 <- index_list[[1]] %>% mutate(survey = "hss")
index_combined <- index_list[[2]] %>% 
  mutate(survey = "ipes") %>% 
  filter(year > 2016) %>% 
  rbind(., index1) %>% 
  filter(year %in% summer_years)

comb_index_plot <- ggplot(index_combined, aes(year, est)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr, fill = survey), shape = 21) +
  labs(x = "Year", y = "Count") +
  ggsidekick::theme_sleek() +
  facet_wrap(~species, scales = "free_y")


pdf(here::here("figs", "diagnostics", "st_index_surv_all_sp.pdf"), height = 7,
    width = 9)
index_plot
comb_index_plot
dev.off()