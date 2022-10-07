### Juvenile Chinook model fit -- testing
## Original model fitting with negative binomial had convergence issues 
## when fitting negbin with volume swept offset; switch to CPUE 
## See previous script for testing of different data inputs and anisotrop vs.
## barrier mesh comparison
## 1) Evaluate link functions (tweedie vs. delta LN vs. delta gamma)
## 2) Fit and generate annual index
## May 26, 2022


library(tidyverse)
library(sdmTMB)
library(ggplot2)


# downscale data and predictive grid
dat <- readRDS(here::here("data", "chin_catch_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    volume_m3 = volume_km3 * 1000,
    # effort = log(volume_m3),
    week = lubridate::week(date),
    vessel = as.factor(vessel),
    ck_juv_cpue = ck_juv / volume_m3
  ) %>% 
  droplevels()

dat_trim <- dat %>%
  filter(# sampling coverage very sparse early in time series
         !year < 1998,
         year < 2018
         ) %>% 
  sample_n(., 250)

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -137, ymin = 47, xmax = -121.25, ymax = 57) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))



## MAKE MESH -------------------------------------------------------------------

# have to account for landmasses so use predictive grid as baseline 
# NOTE: switch to anisotropy given preliminary analysis of improved model fit

spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
                   cutoff = 10, type = "kmeans")

# add barrier mesh
bspde <- add_barrier_mesh(
  spde, coast_utm, range_fraction = 0.1,
  # scaling = 1000 since UTMs were rescaled above
  proj_scaling = 1000, plot = TRUE
)

mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
mesh_df_land <- bspde$mesh_sf[bspde$barrier_triangles, ]
ggplot(coast_utm) +
  geom_sf() +
  geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
  geom_sf(data = mesh_df_land, size = 1, colour = "green") +
  geom_point(data = dat_trim, aes(x = utm_x, y = utm_y),
             colour = "red", alpha = 0.4)


## FIT SIMPLE MODEL ------------------------------------------------------------

library(furrr)
plan(multisession, workers = 6)

# Fit simplified spatial model to identify ideal link function
fam_list <- list(
  tweedie(link = "log"),
  delta_lognormal(),
  delta_gamma(link1 = "logit", link2 = "log")
)
fit_list <- future_map(fam_list, function (x) {
  sdmTMB(
    ck_juv_cpue ~ 1 +  
    s(dist_to_coast_km, bs = "tp", k = 4) +
    s(month, bs = "cc", k = 4) +
    survey_f
    ,
    family = x,
    data = dat_trim,
    mesh = bspde,
    spatial = "on",
    anisotropy = FALSE)
})



### FIT SATURATED --------------------------------------------------------------

# keep early years in dataset for now but exclude 2022 and consider dropping 
# years with v. limited summer data
dat_in <- dat %>% filter(!year == "2022")
spde <- make_mesh(dat_in, c("utm_x_1000", "utm_y_1000"), 
                  cutoff = 10, type = "kmeans")
bspde <- add_barrier_mesh(
  spde, coast_utm, range_fraction = 0.1,
  # scaling = 1000 since UTMs were rescaled above
  proj_scaling = 1000, plot = TRUE
)


fit <- sdmTMB(
  ck_juv ~ 1 +  
    s(bath_depth_mean_m, bs = "tp", k = 4) +
    s(dist_to_coast_km, bs = "tp", k = 4) + 
    s(month, bs = "cc", k = 4) +
    survey_f,
  offset = dat_in$effort,
  data = dat_in,
  mesh = bspde,
  time = "year",
  # infill 1996
  extra_time = c(1996L),
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
  )
)
saveRDS(fit, here::here("data", "fits", "fit_st_full.rds"))
#AIC = 16874.48

fit_v <- sdmTMB(
  ck_juv ~ 1 +  
    s(bath_depth_mean_m, bs = "tp", k = 4) +
    s(dist_to_coast_km, bs = "tp", k = 4) + 
    s(month, bs = "cc", k = 4) +
    survey_f  +
    (1 | vessel),
  offset = dat_in$effort,
  data = dat_in,
  mesh = spde,
  time = "year",
  # infill 1996
  extra_time = c(1996L),
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = FALSE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 10, sigma_lt = 80)
  )
)


## simulate and check residuals
sims <- simulate(fit, nsim = 250)
dharma_residuals(sims, fit)

samps <- sample.int(250, size = 9)
breaks_vec <- c(seq(0, 200, by = 10))
obs_dat <- dat_in %>% filter(ck_juv < 200) %>% pull(ck_juv)
par(mfrow = c(3, 3))
for (i in seq_along(samps)) {
  dum <- sims[, i ]
  hist(obs_dat, col="green", pch=20, cex=4, breaks=breaks_vec)
  hist(dum[dum < 200], pch=20, cex=4, breaks=breaks_vec, col=rgb(1,0,0,0.5), add=TRUE)
}

for (i in seq_along(samps)) {
  dum <- sims[, i ]
  plot(log(dum[-1]) ~ log(dat_in$ck_juv))
  abline(0, 1, col = "red")
}




visreg::visreg(fit2, xvar = "month", scale = "response")


pdf(here::here("figs", "diagnostics", "fits_st.pdf"))
simulate(fit_st, nsim = 500) %>% 
  dharma_residuals(fit_st)
visreg::visreg(fit_st, xvar = "season_f", scale = "response")
visreg::visreg(fit_st, xvar = "dist_to_coast_km", by = "season_f", 
               xlim = c(0, 200))
visreg::visreg(fit_st, xvar = "dist_to_coast_km", by = "season_f", 
               scale = "response", xlim = c(0, 200), nn = 200)
dev.off()


# spatial distribution of residuals
dat_trim$resids <- residuals(fit_st)

pdf(here::here("figs", "diagnostics", "spatial_resids_st.pdf"))
ggplot(dat_trim, aes(utm_x_1000, utm_y_1000, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~year) +
  coord_fixed()
dev.off()

ggplot(dat_trim, aes(x = year, y = resids)) +
  geom_point() +
  ggsidekick::theme_sleek()

grid <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  mutate(utm_x_1000 = X / 1000,
         utm_y_1000 = Y / 1000,
         dist_to_coast_km = shore_dist / 1000)

# add unique years and seasons
exp_grid <- expand.grid(
  year = unique(dat_trim$year),
  season_f = unique(dat_trim$season_f),
  survey_f = unique(dat_trim$survey_f)
) %>%
  mutate(id = row_number()) %>%
  split(., .$id) %>% 
  purrr::map(., function (x) {
    grid %>% 
      mutate(
        year = x$year,
        season_f = x$season_f,
        survey_f = x$survey_f
      )
  }) %>%
  bind_rows() %>% 
  # left_join(., 
  #           dat_trim %>% dplyr::select(year, survey_f, season_f) %>% distinct(),
  #           by = c("year", "season_f")) %>% 
  filter(!year %in% c("1995", "2022")) %>% 
  mutate(
    fake_survey = case_when(
      season_f == "su" & year > 2016 & survey_f == "hss" ~ "1",
      season_f == "su" & year < 2017 & survey_f == "ipes" ~ "1",
      TRUE ~ "0")
  )

preds <- predict(fit_st_survey, newdata = exp_grid)
preds_trim <- preds %>% 
  filter(season_f == "su") 

plot_map <- function(dat, column) {
  ggplot(dat, aes_string("utm_x_1000", "utm_y_1000", fill = column)) +
    geom_raster() +
    coord_fixed() +
    ggsidekick::theme_sleek()
}

plot_map(preds_trim, "exp(est)") +
  scale_fill_viridis_c(
    trans = "sqrt",
    # trim extreme high values to make spatial variation more visible
    na.value = "yellow", limits = c(0, quantile(exp(preds_trim$est), 0.995))
  ) +
  facet_wrap(~year) +
  ggtitle("Prediction (fixed effects + all random effects)")

plot_map(preds_trim, "exp(est_non_rf)") +
  scale_fill_viridis_c(
    trans = "sqrt"
    ) +
  ggtitle("Prediction (fixed effects only)")  

plot_map(preds_trim, "est_rf") +
  scale_fill_gradient2(
  ) +
  ggtitle("Prediction (spatial random effects only)")  

plot_map(preds_trim, "epsilon_st") +
  scale_fill_gradient2() +
  facet_wrap(~year) +
  ggtitle("Prediction (spatiotemporal random effects)")


# index from summer
p_st <- predict(fit_st_survey, 
                newdata = exp_grid %>% filter(season_f == "su",
                                              survey_f == "hss"), 
                return_tmb_object = TRUE)
index <- get_index(p_st)

index_plot <- ggplot(index, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  labs(x = "Year", y = "Count") +
  ggsidekick::theme_sleek()

pdf(here::here("figs", "diagnostics", "st_index_surv.pdf"))
index_plot
dev.off()


# as above but not accounting for survey effects
p_st2 <- predict(fit_st_survey, 
                newdata = exp_grid %>% filter(season_f == "su",
                                              fake_survey == "0"), 
                return_tmb_object = TRUE)
index2 <- get_index(p_st2)
index_list <- list(survey_eff = index, no_survey_eff = index2)


## export to Rmd
saveRDS(exp_grid, here::here("data", "spatial", "exp_pred_ipes_grid_utm.rds"))
saveRDS(preds, here::here("data", "preds", "st_survey_spatial_preds.rds"))
saveRDS(index_list, here::here("data", "preds", "st_survey_index_list.rds"))


## ANISOTROPY VS BARRIER MESH --------------------------------------------------

# use summer dataset to evaluate relative support for anisotropy vs. barrier 
# mesh

fit_summer_barrier <- sdmTMB(ck_juv ~ 1 +  
                       s(dist_to_coast_km, bs = "tp", k = 4) + 
                       # switch to thin plate 
                       s(month, bs = "tp", k = 4) +
                       survey_f,
                     offset = dat_summer$effort,
                     data = dat_summer %>% filter(year < 2018),
                     mesh = bspde_summer,
                     time = "year",
                     extra_time = c(2016L, 2018L, 2019L),
                     family = sdmTMB::nbinom2(),
                     spatial = "on",
                     spatiotemporal = "ar1"
)
fit_summer_anisotropy <- sdmTMB(ck_juv ~ 1 +  
                               s(dist_to_coast_km, bs = "tp", k = 4) + 
                               # switch to thin plate 
                               s(month, bs = "tp", k = 4) +
                               survey_f,
                             offset = dat_summer$effort,
                             data = dat_summer %>% filter(year < 2018),
                             anisotropy = TRUE,
                             mesh = spde_summer,
                             time = "year",
                             extra_time = c(2016L, 2018L, 2019L),
                             family = sdmTMB::nbinom2(),
                             spatial = "on",
                             spatiotemporal = "ar1"
)
AIC(fit_summer_anisotropy, fit_summer_barrier)

summer_fit_tbl <- tibble(
  name = c("barrier", "anisotropy"
  ),
  mod = list(fit_summer_barrier,
             fit_summer_anisotropy)) %>% 
  mutate(
    preds = purrr::map(mod, function (x) {
      predict(x, newdata = dat_test)
    }),
    test_preds = purrr::map(preds, function (x) {
      x %>% 
        filter(year %in% c("2018", "2019"),
               ipes_grid == TRUE,
               season_f == "su")
    })
  )

# calculate RMSE and plot predictions
purrr::map(summer_fit_tbl$test_preds, function (x) {
  Metrics::rmse(x$ck_juv, exp(x$est))
}) %>% 
  unlist()

plot_list <- purrr::map(summer_fit_tbl$preds, function (x) {
  plot(ck_juv ~ exp(est), data = x %>% filter(ck_juv < 100))
  abline(a = 0, b = 1, col = "red")
})

hist_list <- purrr::map(summer_fit_tbl$preds, function (x) {
  dum <- x %>% filter(ck_juv < 50)
  p1 <- hist(exp(dum$est))                     
  p2 <- hist(dum$ck_juv)                     
  plot( p1, col=rgb(0,0,1,1/4)#, xlim=c(0,20)
        )  # first histogram
  plot( p2, col=rgb(1,0,0,1/4)#, xlim=c(0,20)
        , add=T)  # second
})


# compare spatial predictions
grid <- readRDS(here::here("data", "spatial", "pred_ipes_grid.RDS")) %>% 
  mutate(utm_x_1000 = X / 1000,
         utm_y_1000 = Y / 1000,
         dist_to_coast_km = shore_dist / 1000)

# add unique years and seasons
summer_exp_grid <- expand.grid(
  year = seq(1998, 2019, by = 1)
) %>%
  mutate(id = row_number()) %>%
  split(., .$id) %>% 
  purrr::map(., function (x) {
    grid %>% 
      mutate(
        year = x$year,
        month = 7,
        survey_f = "hss"
      )
  }) %>%
  bind_rows() 

summer_fit_tbl$spatial_preds <- purrr::map2(
  summer_fit_tbl$name,
  summer_fit_tbl$mod, 
  function (x, y) {
    predict(y, summer_exp_grid) %>% 
      mutate(name = x)
  }
)


preds_2015 <- summer_fit_tbl$spatial_preds %>% 
  bind_rows() %>% 
  filter(year == "2015")

plot_map <- function(dat, column) {
  ggplot(dat, aes_string("utm_x_1000", "utm_y_1000", fill = column)) +
    geom_raster() +
    coord_fixed() +
    ggsidekick::theme_sleek()
}

out_map <- plot_map(preds_2015, "exp(est)") +
  scale_fill_viridis_c(
    trans = "sqrt",
    na.value = "yellow", limits = c(0, quantile(exp(preds_2015$est), 0.995))
  ) +
  facet_wrap(~name) +
  ggtitle("Prediction (fixed effects + all random effects)")


png(here::here("figs", "anisotropy_barrier_comparison_map.png"), 
    height = 5, width = 8, res = 250, 
    units = "in")
out_map
dev.off()