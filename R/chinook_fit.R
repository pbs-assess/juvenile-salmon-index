### Juvenile Chinook model fit
## Fit sdmTMB to 
## 1) Develop standardized index of abundance
## 2) Estimate vessel effects (requires developing interannual smooth or 
## AR-term, otherwise confounded with year effects)
## May 26, 2022


library(tidyverse)
library(sdmTMB)
library(ggplot2)


# downscale data and predictive grid
dat_trim <- readRDS(here::here("data", "chin_catch_sbc.rds")) %>% 
  mutate(utm_x_1000 = utm_x / 1000,
         utm_y_1000 = utm_y / 1000,
         effort = log(distance_travelled),
         week = lubridate::week(date)
         ) %>%
  filter(#!is.na(depth_mean_m),
         #!is.na(dist_to_coast_km),
         !effort < 0,
         !is.na(effort),
         # drop 1995 so that '96 doesn't have to be interpolated
         !year %in% c("1995", "1997"),
         year < 2018)

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -137, ymin = 47, xmax = -121.25, ymax = 57) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))



## MAKE MESH -------------------------------------------------------------------

# have to account for landmasses so use predictive grid as baseline 
# NOTE: switch to anisotropy given preliminary analysis of improved model fit

# spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
#                            n_knots = 250, type = "kmeans")
spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
                   cutoff = 10, type = "kmeans")
plot(spde)

# spde2 <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
#                   cutoff = 20, type = "kmeans")
# plot(spde2)

# add barrier mesh
# bspde <- add_barrier_mesh(
#   spde, coast_utm, range_fraction = 0.1,
#   # scaling = 1000 since UTMs were rescaled above
#   proj_scaling = 1000, plot = TRUE
# )
# 
# mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
# mesh_df_land <- bspde$mesh_sf[bspde$barrier_triangles, ]
# ggplot(coast_utm) +
#   geom_sf() +
#   geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
#   geom_sf(data = mesh_df_land, size = 1, colour = "green") +
#   geom_point(data = dat_trim, aes(x = utm_x, y = utm_y), 
#              colour = "red", alpha = 0.4)


#compare 
# plot(dat_trim_spde)
# plot(jchin1_spde_nch)
# plot(jchin1_spde$mesh, main = NA, edge.color = "grey60", asp = 1)
# plot(jchin1_spde_nch$mesh, main = NA, edge.color = "grey60", asp = 1)


## alternative meshes that use data subsets
dat_ipes <- dat_trim %>% filter(ipes_grid == TRUE)
spde_ipes <-  make_mesh(dat_ipes,
                        c("utm_x_1000", "utm_y_1000"), 
                        cutoff = 10, type = "kmeans")
# bspde_ipes <- add_barrier_mesh(
#   spde_ipes, coast_utm, range_fraction = 0.1,
#   proj_scaling = 1000, plot = TRUE
# )

dat_summer <- dat_trim %>% filter(season_f == "su")
spde_summer <-  make_mesh(dat_summer,
                          c("utm_x_1000", "utm_y_1000"), 
                          cutoff = 10, type = "kmeans")
# bspde_summer <- add_barrier_mesh(
#   spde_summer, coast_utm, range_fraction = 0.1,
#   proj_scaling = 1000, plot = TRUE
# )

dat_summer_ipes <- dat_trim %>% filter(season_f == "su",
                                       ipes_grid == TRUE)
spde_summer_ipes <-  make_mesh(dat_summer_ipes,
                               c("utm_x_1000", "utm_y_1000"), 
                               cutoff = 10, type = "kmeans")
# bspde_summer_ipes <- add_barrier_mesh(
#   spde_summer_ipes, coast_utm, range_fraction = 0.1,
#   proj_scaling = 1000, plot = TRUE
# )

par(mfrow = c(2, 2))
plot(spde)
plot(spde_ipes)
plot(spde_summer)
plot(spde_summer_ipes)



## FIT SIMPLE MODEL ------------------------------------------------------------

# correlation between depth and distance to coast may be significant; probably
# should fit one or the other
dum <- dat_trim[, c("depth_mean_m", "dist_to_coast_km")]
cor(dum)

# Fit simplified spatial model to most constrained dataset to estimate relative 
# benefits of NB vs zero-inflated NB 
# (exclude spatiotemporal to minimize fit time)
# fit_nb <- sdmTMB(ck_juv ~ s(dist_to_coast_km, k = 4) +  survey_f, 
#                  offset = dat_summer$effort,
#                  data = dat_summer,
#                  mesh = bspde_summer,
#                  family = sdmTMB::nbinom2(),
#                  spatial = "on"
# )
# 
# # fails to converge
# fit_nb0 <- sdmTMB(ck_juv ~ s(dist_to_coast_km, k = 4) + survey_f, 
#                  offset = dat_summer$effort,
#                  data = dat_summer,
#                  mesh = bspde_summer,
#                  family = sdmTMB::delta_truncated_nbinom2(),
#                  spatial = "on",
#                  control = sdmTMBcontrol(
#                    nlminb_loops = 2,
#                    newton_loops = 2
#                  )
# )
# 
# # inferior AIC
# fit_nb1 <- sdmTMB(ck_juv ~ s(dist_to_coast_km, k = 4) +  survey_f, 
#                   offset = dat_summer$effort,
#                   data = dat_summer,
#                   mesh = bspde_summer,
#                   family = sdmTMB::nbinom1(),
#                   spatial = "on"
# )
# # inferior AIC (even with ST version)
# fit_ln0 <- sdmTMB(ck_juv ~ s(dist_to_coast_km, k = 4) + survey_f,
#                   offset = dat_summer$effort,
#                   data = dat_summer,
#                   mesh = bspde_summer,
#                   family = sdmTMB::delta_lognormal(),
#                   spatial = "on",
#                   time = "year",
#                   extra_time = c(2016L),
#                   spatiotemporal = "ar1",
#                   control = sdmTMBcontrol(
#                     nlminb_loops = 2,
#                     newton_loops = 1
#                   )
# )


## COMPARE FIXED EFFECT STRUCTURES ---------------------------------------------

# week vs month (spatial only)
fit_month <- sdmTMB(ck_juv ~ 1 +  
                 s(dist_to_coast_km, bs = "tp", k = 4) + 
                 s(month, bs = "cc", k = 4) +
                 survey_f,
               offset = dat_trim$effort,
               data = dat_trim,
               mesh = bspde,
               family = sdmTMB::nbinom2(),
               spatial = "on"
               )
fit_week <- sdmTMB(ck_juv ~ 1 +  
                 s(dist_to_coast_km, bs = "tp", k = 4) + 
                 s(week, bs = "cc", k = 4) +
                 survey_f,
               offset = dat_trim$effort,
               data = dat_trim,
               mesh = bspde,
               family = sdmTMB::nbinom2(),
               spatial = "on"
)
AIC(fit_month, fit_week)
# month outperforms week

# bathy vs. distance
fit_dist <- fit_month
fit_bathy <- sdmTMB(ck_juv ~ 1 +  
                      s(depth_mean_m, bs = "tp", k = 4) + 
                      s(month, bs = "cc", k = 4) +
                      survey_f,
                    offset = dat_trim$effort,
                    data = dat_trim,
                    mesh = bspde,
                    family = sdmTMB::nbinom2(),
                    spatial = "on"
) 
AIC(fit_bathy, fit_dist)
# distance outperforms bathymetry

# seasonal structure 
fit_a <- fit_dist
fit_b <- sdmTMB(ck_juv ~ 1 +
                s(dist_to_coast_km, bs = "tp", k = 4, by = season_f) +
                s(month, bs = "cc", k = 4) +
                 survey_f,
                offset = dat_trim$effort,
                data = dat_trim,
                mesh = bspde,
                family = sdmTMB::nbinom2(),
                spatial = "on",
                control = sdmTMBcontrol(
                  nlminb_loops = 2#,
                  # newton_loops = 2
                )
)
fit_c <- sdmTMB(ck_juv ~ 1 +  
                 s(dist_to_coast_km, bs = "tp", k = 4, by = season_f) + 
                 season_f +
                 survey_f,
                offset = dat_trim$effort,
                data = dat_trim,
                mesh = bspde,
                family = sdmTMB::nbinom2(),
                spatial = "on",
                control = sdmTMBcontrol(
                  nlminb_loops = 2
                )
)
AIC(fit_a, fit_b, fit_c)
# seasonal interactions not supported with current compelx structure


## COMPARE DATA INPTUS ---------------------------------------------------------

# compare predictions from four different input data structures, then compare to
# observed
# a) coastwide data, all seasons
# b) coastwide data, summer only
# c) IPES data, all seasons
# d) IPES data, summer only


# fit models with different fixed effects due to different time scales
fit_full <- sdmTMB(
  ck_juv ~ 1 +  
    s(dist_to_coast_km, bs = "tp", k = 4) + 
    s(month, bs = "cc", k = 4) +
    survey_f,
  offset = dat_trim$effort,
  data = dat_trim %>% filter(year < 2018),
  mesh = spde,
  time = "year",
  extra_time = c(2018L, 2019L),
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 2, sigma_lt = 5)
  )
)
fit_summer <- sdmTMB(
  ck_juv ~ 1 +  
    s(dist_to_coast_km, bs = "tp", k = 4) + 
    # switch to thin plate 
    s(month, bs = "tp", k = 4) +
    survey_f,
  offset = dat_summer$effort,
  data = dat_summer %>% filter(year < 2018),
  mesh = spde_summer,
  time = "year",
  extra_time = c(2016L, 2018L, 2019L),
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 2, sigma_lt = 5)
  )
)
fit_ipes <- sdmTMB(
  ck_juv ~ 1 +  
    s(dist_to_coast_km, bs = "tp", k = 4) + 
    s(month, bs = "cc", k = 4) +
    survey_f,
  offset = dat_ipes$effort,
  data = dat_ipes %>% filter(year < 2018),
  mesh = spde_ipes,
  time = "year",
  extra_time = c(2018L, 2019L),
  family = sdmTMB::nbinom2(),
  spatial = "on",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 2, sigma_lt = 5)
  )
)
# fit_summer_ipes <- sdmTMB(
#   ck_juv ~ 1 +  
#     s(dist_to_coast_km, bs = "tp", k = 4) + 
#     s(month, bs = "tp", k = 4) +
#     survey_f,
#   offset = dat_summer_ipes$effort,
#   data = dat_summer_ipes %>% filter(year < 2018),
#   mesh = spde_summer_ipes,
#   time = "year",
#   extra_time = c(2016L, 2018L, 2019L),
#   family = sdmTMB::nbinom2(),
#   spatial = "on",
#   spatiotemporal = "ar1",
#   anisotropy = TRUE,
#   priors = sdmTMBpriors(
#     matern_s = pc_matern(range_gt = 2, sigma_lt = 5)
#   )
# )

fit_list <- list(fit_full,
                 fit_summer,
                 fit_ipes#,
                 # fit_summer_ipes
                 )

saveRDS(fit_list, here::here("data", "fit_list.rds"))

## testing data are observed catches in 2018/19 in summer in the ipes_grid
# due to sdmTMB constraints, all time elements have to be passed when generating
# predictions
dat_test <- readRDS(here::here("data", "chin_catch_sbc.rds")) %>% 
  mutate(utm_x_1000 = utm_x / 1000,
         utm_y_1000 = utm_y / 1000,
         effort = log(distance_travelled)
  ) %>%
  filter(
    !effort < 0,
    !is.na(effort),
    !year %in% c("1995", "1997"),
    year < 2020)


# generate predictions for each model/dataset
fit_tbl <- tibble(
  name = c("full", "sum", "ipes"#, "sum_ipes"
  ),
  mod = fit_list) %>% 
  mutate(
    preds = purrr::map(mod, function (x) {
      test_preds <- predict(x,
                            newdata = dat_test)
      
      test_preds %>% 
        filter(year %in% c("2018", "2019"),
               ipes_grid == TRUE,
               season_f == "su")
    })
  )

# calculate RMSE and plot predictions
# fit_tbl$rmse <- purrr::map(fit_tbl$preds, function (x) {
#   Metrics::rmse(x$ck_juv, exp(x$est))
# }) %>% 
#   unlist()

# calculate log-lik
fit_tbl$log_lik <- purrr::map(fit_tbl$preds, function (x) {

})



plot_list <- purrr::map(fit_tbl$preds, function (x) {
  plot(ck_juv ~ exp(est), data = x)
})





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