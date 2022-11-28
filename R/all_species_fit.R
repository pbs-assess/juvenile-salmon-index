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
      n_juv ~ dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        s(target_depth, bs = "tp", k = 4) +
        day_night +
        survey_f,
      offset = x$effort,
      data = x,
      mesh = spde_in,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      anisotropy = FALSE,
      priors = sdmTMBpriors(
        matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
      ),
      knots = list(
        week = c(0, 52)
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


### FIT SATURATED --------------------------------------------------------------

st_mod <- furrr::future_pmap(
  list(dat_tbl$data, dat_tbl$bspde, dat_tbl$species),
  function (x, spde_in, sp_in) {
    sdmTMB(
      n_juv ~ 1 +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
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
      knots = list(
        week = c(0, 52)
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

purrr::map(dat_tbl$sims, ~ mean(.x[,1]))


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
 week = seq(0, 52, length.out = 100),
 dist_to_coast_km = median(dat_in$dist_to_coast_km),
 day_night = "DAY",
 target_depth = 0,
 survey_f = "hss",
 year = 2012L
)
dist_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = seq(min(dat_in$dist_to_coast_km), 
                         max(dat_in$dist_to_coast_km), length.out = 100),
  day_night = "DAY",
  target_depth = 0,
  survey_f = "hss",
  year = 2012L
)
day_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = c("DAY", "NIGHT"),
  target_depth = 0,
  survey_f = "hss",
  year = 2012L
)
target_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = "DAY",
  target_depth = seq(min(dat_in$target_depth), 
                         max(dat_in$target_depth), length.out = 100),
  survey_f = "hss",
  year = 2012L
)
survey_dat <- data.frame(
  week = median(dat_in$week),
  dist_to_coast_km = median(dat_in$dist_to_coast_km),
  day_night = "DAY",
  target_depth = 0,
  survey_f = c("hss", "ipes"),
  year = 2012L
)

# pred_tbl <- tibble(
#   var = c("week",
#           "dist_to_coast_km",
#           "day_night",
#           "target_depth",
#           "survey_f"),
#   data = list(week_dat, dist_dat, day_dat, target_dat, survey_dat),
#   plot = c("line", "line", "dot", "line", "dot")
# )
# 
# # predict
# pred_list <- vector(length = nrow(pred_tbl), mode = "list")
# spatial_pred_list <- vector(length = nrow(pred_tbl), mode = "list")
# for (i in seq_along(pred_tbl$var)) {
#   pred_list[[i]] <- furrr::future_map2(
#     dat_tbl$st_mod, dat_tbl$species, function(x , sp) {
#       predict(x, newdata = pred_tbl$data[[i]], se_fit = T, re_form = NA) %>%
#         mutate(species = sp,
#                up = est + 1.96 * est_se,
#                lo = est - 1.96 * est_se,
#                exp_est = exp(est),
#                exp_up = exp(up),
#                exp_lo = exp(lo))
#     }
#   ) %>%
#     bind_rows()
#   spatial_pred_list[[i]] <- furrr::future_map2(
#     dat_tbl$spatial_mod, dat_tbl$species, function(x , sp) {
#       predict(x, newdata = pred_tbl$data[[i]], se_fit = T, re_form = NA) %>%
#         mutate(species = sp,
#                up = est + 1.96 * est_se,
#                lo = est - 1.96 * est_se,
#                exp_est = exp(est),
#                exp_up = exp(up),
#                exp_lo = exp(lo))
#     }
#   ) %>%
#     bind_rows()
# }
# 
# pred_tbl$pred_dat <- pred_list
# pred_tbl$spatial_pred_dat <- spatial_pred_list

# takes a long time so save
# saveRDS(pred_tbl, here::here("data", "preds", "fe_preds.rds"))
pred_tbl <- readRDS(here::here("data", "preds", "fe_preds.rds"))


# make all plots for quick visualization
plot_eff <- function (dat, x_var, type = c("line", "dot")) {
  p <- dat %>% 
    # mutate(
    #   exp_est = exp(est),
    #   up = exp(est + 1.96 * est_se),
    #   lo = exp(est - 1.96 * est_se)
    # ) %>% 
    ggplot(
    ., 
    aes_string(x_var, "exp_est", ymin = "lo", ymax = "up")
  ) +
    facet_wrap(~species, scales = "free_y") +
    labs(title = x_var) +
    ggsidekick::theme_sleek()
  if (type == "line") {
    p2 <- p +
      geom_line() +
      geom_ribbon(alpha = 0.3)
  }
  if (type == "dot") {
    p2 <- p + 
      geom_pointrange() 
  }
  return(p2)
} 

pdf(here::here("figs", "fixed_effects.pdf"), width = 8, height = 5)
purrr::pmap(list(pred_tbl$pred_dat, pred_tbl$var, pred_tbl$plot), plot_eff)
dev.off()

pdf(here::here("figs", "spatial_fixed_effects.pdf"), width = 8, height = 5)
purrr::pmap(list(pred_tbl$spatial_pred_dat, pred_tbl$var, pred_tbl$plot), 
            plot_eff)
dev.off()


# make individual plots for ms
pred_tbl$adj_pred_dat <- purrr::map(
  pred_tbl$pred_dat,
  ~ .x %>% 
    mutate(
      up = est + 1.96 * est_se,
      lo = est - 1.96 * est_se,
      exp_est = exp(est),
      exp_up = exp(up),
      exp_lo = exp(lo)
    )
) 
p_dat <- pred_tbl %>% 
  select(var, adj_pred_dat) %>% 
  unnest(cols = adj_pred_dat) %>% 
  mutate(species = tolower(species),
         day_night = tolower(day_night),
         survey_f = toupper(survey_f))

col_pal <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0')
names(col_pal) <- c('chinook','pink','chum','coho','sockeye')

# line plot function
pp_foo <- function(dat, x_var, y_var = "exp_est", ymin_var = "exp_lo", 
                     ymax_var = "exp_up", fill_var = "species") {
  ggplot(
    dat,
    aes_string(x_var, y_var, ymin = ymin_var, ymax = ymax_var,
               fill = fill_var)
  ) +
    ggsidekick::theme_sleek() +
    scale_fill_manual(values = col_pal) +
    ylab("Abundance Index")
}

png(here::here("figs", "ms_figs", "week_preds.png"), height = 4, width = 5,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "week"),
  x_var = "week"
) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Week") 
dev.off()


png(here::here("figs", "ms_figs", "depth_preds.png"), height = 4, width = 5,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "target_depth"),
  x_var = "target_depth"
) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Target Headrope Depth (m)") 
dev.off()


png(here::here("figs", "ms_figs", "dist_preds.png"), height = 4, width = 5,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "dist_to_coast_km"),
  x_var = "dist_to_coast_km"
) +
  geom_line() +
  geom_ribbon(alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Distance to Nearest Coastline (km)") 
dev.off()


png(here::here("figs", "ms_figs", "dn_preds.png"), height = 3, width = 8,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "day_night"),
  x_var = "day_night"
) +
  geom_pointrange(shape = 21) +
  facet_wrap(~species, nrow = 1) +
  xlab("Diel Effects") +
  theme(
    legend.position = "none"
  )
dev.off()


png(here::here("figs", "ms_figs", "surv_preds.png"), height = 3, width = 8,
    units = "in", res = 200)
pp_foo(
  dat = p_dat %>% filter(var == "survey_f"),
  x_var = "survey_f"
) +
  geom_pointrange(shape = 21) +
  facet_wrap(~species, nrow = 1) +
  xlab("Survey Effects") +
  theme(
    legend.position = "none"
  )
dev.off()


## MAKE SPATIAL PREDICTIONS ----------------------------------------------------

# spatial distribution of residuals
grid <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  mutate(utm_x_1000 = X / 1000,
         utm_y_1000 = Y / 1000,
         dist_to_coast_km = shore_dist / 1000)

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(dat_in$year),
  survey_f = unique(dat_in$survey_f),
  week = c(25, 42)
) %>%
  mutate(id = row_number()) %>%
  split(., .$id) %>%
  purrr::map(., function (x) {
    grid %>% 
      mutate(
        year = x$year,
        survey_f = x$survey_f,
        target_depth = 0,
        week = x$week,
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


dat_tbl$spatial_preds <- purrr::map(dat_tbl$st_mod, function (x) {
  predict(x, newdata = exp_grid, se_fit = FALSE, re_form = NULL)
})

stock_preds <- dat_tbl %>%
  select(species, spatial_preds) %>%
  unnest(cols = spatial_preds) %>%
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

dd <- stock_preds %>% 
  filter(year == "2012") %>% 
  mutate(
    exp_est = exp(est)
  )
plot_map(dd,
         "est") +
  scale_fill_viridis_c(
    # trans = "sqrt",
    # limits = c(0, quantile(exp(dd$est), 0.995))
  ) +
  facet_wrap(~species)


# index from summer
summer_years <- dat_in %>% 
  filter(season_f == "su",
         #remove 2021 survey years (too few tows)
         !year == "2021") %>%
  pull(year) %>% 
  unique()
fall_years <- dat_in %>% 
  filter(season_f == "wi",
         #remove 2021 survey years (too few tows)
         !year == "2021") %>%
  pull(year) %>% 
  unique()


# fix to HSS survey (can't combine because predictions shouldn't be passed 
# duplicates, but require tmb_object stored in preds)
# ind_preds_sum <- purrr::map(dat_tbl$st_mod, function (x) {
#   predict(x, newdata = exp_grid %>% filter(survey_f == "hss", week == "25"), 
#           return_tmb_object = TRUE)
# })
# index_list_sum <- purrr::map(ind_preds_sum, get_index, bias_correct = TRUE)
# 
# index_df_summ <- purrr::map2(index_list_sum, dat_tbl$species, function (x, sp) {
#   x$species <- sp
#   return(x)
# }) %>% 
#   bind_rows()
# 
# 
# ind_preds_fall <- purrr::map(dat_tbl$st_mod, function (x) {
#   predict(x, newdata = exp_grid %>% filter(survey_f == "hss", week == "42"), 
#           return_tmb_object = TRUE)
# })
# index_list_fall <- purrr::map(ind_preds_fall, get_index, bias_correct = TRUE)
# 

# index_lists <- c(index_list_sum,
#                     index_list_fall)
# saveRDS(index_lists, here::here("data", "fits", "index_list.rds"))
index_lists <- readRDS(here::here("data", "fits", "index_list.rds")) 


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
    group_mean = mean(est),
    scaled_est = scale(est, center = TRUE, scale = TRUE) %>% 
      as.numeric(),
    anomaly = case_when(
      lwr > group_mean ~ "pos",
      upr < group_mean ~ "neg",
      TRUE ~ "avg"
    ) %>% 
      as.factor()
    ) %>% 
  ungroup()


index_plot <- ggplot(index_dat, 
                     aes(year, est)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), shape = 21, fill = "white") +
  labs(x = "Year", y = "Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") 


col_pal <- c("#f5f5f5", "#d8b365", "#5ab4ac")
names(col_pal) <- levels(index_dat$anomaly)

index_scaled <- ggplot(index_dat) +
  geom_point(aes(year, scaled_est, fill = anomaly), shape = 21, size = 2) +
  labs(x = "Year", y = "Scaled Abundance Index") +
  ggsidekick::theme_sleek() +
  facet_grid(species~season, scales = "free_y") +
  scale_fill_manual(values = col_pal)


png(here::here("figs", "ms_figs", "hss_index.png"))
index_plot
dev.off()

png(here::here("figs", "ms_figs", "hss_index_scaled.png"))
index_scaled
dev.off()



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