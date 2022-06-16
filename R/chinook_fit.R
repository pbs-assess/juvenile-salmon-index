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
         survey_f = ifelse(
           year > 2016 & season_f == "su", "ipes", "hss") %>% as.factor()
         ) %>% 
  filter(!is.na(depth_mean_m),
         !is.na(dist_to_coast_km),
         !is.na(effort),
         # drop 1995 so that '96 doesn't have to be interpolated
         !year %in% c("1995", "2022"
                      ))

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -137, ymin = 47, xmax = -121.25, ymax = 57) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))



## MAKE MESH -------------------------------------------------------------------

# have to account for landmasses so use predictive grid as baseline 

# spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
#                            n_knots = 250, type = "kmeans")
spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
                   cutoff = 10, type = "kmeans")
plot(spde)


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


#compare 
plot(dat_trim_spde)
plot(jchin1_spde_nch)
plot(jchin1_spde$mesh, main = NA, edge.color = "grey60", asp = 1)
plot(jchin1_spde_nch$mesh, main = NA, edge.color = "grey60", asp = 1)


# as above but for summer survey only
spde_summer <-  make_mesh(dat_trim %>% filter(season_f == "su"),
                          c("utm_x_1000", "utm_y_1000"), 
                          cutoff = 10, type = "kmeans")
bspde_summer <- add_barrier_mesh(
  spde_summer, coast_utm, range_fraction = 0.1,
  # scaling = 1000 since UTMs were rescaled above
  proj_scaling = 1000, plot = TRUE
)


## FIT SIMPLE MODEL ------------------------------------------------------------

# correlation between depth and distance to coast may be significant; probably
# should fit one or the other
dum <- dat_trim[, c("depth_mean_m", "dist_to_coast_km")]
cor(dum)

# Fit offset only
# TODO: ask Sean why offsets don't work...
fit0 <- sdmTMB(ck_juv ~ depth_mean_m,
              offset = dat_trim$effort,
              data = dat_trim,
              mesh = bspde,
              family = nbinom2(link = "log"),
              spatial = "off")


# Fit spatial only model with environmental covariates and no effort offset
# (month, bathymetry and distance to coast)
# Explored monthly smooth (instead of season), but issues with convergence
fit <- sdmTMB(ck_juv ~ s(depth_mean_m, by = season_f) + 
                s(dist_to_coast_km, by = season_f) + 
                season_f,
              # offset = effort,
              data = dat_trim,
              mesh = bspde,
              # silent = FALSE,
              family = nbinom2(link = "log"),
              spatial = "on")

fit
sanity(fit)


# Fit spatiotemporal model with only one environmental covariate
fit_st <- sdmTMB(ck_juv ~ s(dist_to_coast_km, by = season_f, k = 3) + 
                season_f,
              # offset = dat_trim$effort,
              data = dat_trim,
              mesh = bspde,
              time = "year",
              family = nbinom2(link = "log"),
              spatial = "on",
              spatiotemporal = "ar1")

sanity(fit_st)


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



# compare fit_st to models that include a) fixed effect for survey and b) 
# constrained to only summer months (also requires infilling)
fit_st_survey <- sdmTMB(ck_juv ~ s(dist_to_coast_km, by = season_f, k = 3) + 
                          season_f + survey_f,
                        # offset = dat_trim$effort,
                        data = dat_trim,
                        mesh = bspde,
                        time = "year",
                        family = nbinom2(link = "log"),
                        spatial = "on",
                        spatiotemporal = "ar1")
sanity(fit_st_survey)
#export do explore in Rmd
# saveRDS(fit_st_survey, here::here("data", "fits", "fit_st_survey.RDS"))
fit_st_survey <- readRDS(here::here("data", "fits", "fit_st_survey.RDS"))


# ggpredict throws variable type errors
ggeffects::ggpredict(fit_st_survey,
                     terms = "dist_to_coast_km") %>%
  plot()

season_p <- visreg::visreg(fit_st_survey, xvar = "season_f", 
                           scale = "response")
dist_p <- visreg::visreg(fit_st_survey, xvar = "dist_to_coast_km", 
                         by = "season_f", xlim = c(0, 200), scale = "response")
plot_list <- list(season_p, dist_p)
saveRDS(plot_list, here::here("figs", "st_survey_counterfacs_list.RDS"))



## SPATIAL PREDS ---------------------------------------------------------------

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


