### Explore different mesh configurations for high seas juvenile data
## Assume range = ~120; use pink data since most problematic
## Compare: 
## 1) INLA mesh based on max_edge (<1/5 max)
## 2) sdmTMB n_knots = ~250
## 3) sdmTMB cutoff = ~30
## April 30, 2020
## Updated Nov 28, 2022

library(tidyverse)
library(sdmTMB)
library(ggplot2)

# downscale data and predictive grid
dat <- readRDS(here::here("data", "catch_survey_sbc.rds")) %>% 
  mutate(
    utm_x_1000 = utm_x / 1000,
    utm_y_1000 = utm_y / 1000,
    effort = log(volume_m3),
    week = lubridate::week(date)
  ) %>% 
  droplevels() 

dat_in <- dat %>% 
  filter(!year == "2022",
         !bath_depth_mean_m < 0)

pink_dat <- dat_in %>% 
  filter(grepl("PINK", species))


## Build different meshes ------------------------------------------------------

# INLA meshes
max_edge = 30 # based on 0.2 * 180 
bound_outer = 150


inla_mesh_raw <- INLA::inla.mesh.2d(
  loc = cbind(pink_dat$utm_x_1000, pink_dat$utm_y_1000),
  max.edge = c(1, 5) * max_edge,
  # - use 5 times max.edge in the outer extension/offset/boundary
  cutoff = max_edge / 5,
  offset = c(max_edge, bound_outer)
) 
inla_mesh <- make_mesh(pink_dat, 
                       c("utm_x_1000", "utm_y_1000"),
                       mesh = inla_mesh_raw)
# even though this is barely within the 20% of range estimate they suggest,
# still extremely high resolution

inla_mesh2_raw <- INLA::inla.mesh.2d(
  loc = cbind(pink_dat$utm_x_1000, pink_dat$utm_y_1000),
  max.edge = c(1, 5) * 1000,
  # - use 5 times max.edge in the outer extension/offset/boundary
  cutoff = max_edge / 5,
  offset = c(max_edge, bound_outer)
) 
inla_mesh2 <- make_mesh(pink_dat, 
                       c("utm_x_1000", "utm_y_1000"),
                       mesh = inla_mesh2_raw)
# even with much larger edge size than recommended n > 1000


# sdmTMB meshes
sdm_mesh <- make_mesh(pink_dat,
                  c("utm_x_1000", "utm_y_1000"),
                  type = "kmeans",
                  n_knots = 250)

sdm_mesh2 <- make_mesh(pink_dat,  
                        c("utm_x_1000", "utm_y_1000"),
                        type = "cutoff",
                        cutoff = 30)

mesh_list <- list(inla_mesh, inla_mesh2, sdm_mesh, sdm_mesh2)
names(mesh_list) <- c("inla_vfine", "inla_fine", "sdm_coarse", "sdm_vcoarse")

purrr::map(mesh_list, ~ .x$mesh$n)


## Fit both models -------------------------------------------------------------

## Evaluate mesh impacts using spatial model initial
fit_list <- purrr::map(
  mesh_list, ~ 
    sdmTMB(
      n_juv ~ 1 +
        as.factor(year) +
        dist_to_coast_km +
        s(week, bs = "cc", k = 5) +
        target_depth +
        day_night +
        survey_f,
      offset = pink_dat$effort,
      data = pink_dat,
      mesh = .x,
      family = sdmTMB::nbinom2(),
      spatial = "on",
      spatiotemporal = "iid",
      time = "year",
      anisotropy = TRUE,
      knots = list(
        week = c(0, 52)
      ),
      control = sdmTMBcontrol(
        newton_loops = 1
      )
    )
)

# some differences in intercept estimates 
purrr::map(
  fit_list, ~ .x$sd_report$par.fixed[3]
)

# qqplot of residuals
purrr::map(
  fit_list, 
  ~ {
    resid <- residuals(.x)
    qqnorm(resid)
    abline(a = 0, b = 1, col = "red")
  }
)

# some differences in intercept estimates 
purrr::map(
  fit_list, sanity
)


## SPATIAL PREDICTIONS ---------------------------------------------------------

# grid of summer survey area
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

hss_pred <- exp_grid %>% 
  filter(
    # use 2012 as arbitrary reference year
    year %in% c("2010", "2011", "2012"),
    week == "42",
    survey_f == "hss"
  )

# generate predictions for each mesh type fixed to reference values for cov
pred_list <- purrr::map2(
  fit_list, names(mesh_list),
  ~ predict(.x, newdata = hss_pred, se_fit = FALSE, re_form = NULL) %>% 
    mutate(
      mesh = .y
    )
)
pred_dat <- pred_list %>% 
  bind_rows() %>% 
  mutate(
    exp_est = exp(est)
  )


plot_map <- function(dat, column) {
  ggplot(dat, aes_string("utm_x_1000", "utm_y_1000", fill = column)) +
    geom_raster() +
    coord_fixed() +
    ggsidekick::theme_sleek()
}

mean_pred <- plot_map(pred_dat, "exp(est)") +
  scale_fill_viridis_c(
    trans = "sqrt",
    limits = c(0, quantile(exp(pred_dat$est), 0.995))
  ) +
  facet_wrap(year~mesh) +
  ggtitle("Prediction (fixed effects + all random effects)")
omega_pred <- plot_map(pred_dat %>% filter(year == "2012"), "omega_s") +
  scale_fill_gradient2() +
  facet_wrap(~mesh) +
  ggtitle("Prediction (spatial random effects)")
epsilon_pred <- plot_map(pred_dat, "epsilon_st") +
  scale_fill_gradient2() +
  facet_wrap(year~mesh) +
  ggtitle("Prediction (spatiotemporal random effects)")

pdf(here::here("figs", "mesh_comparison", "mesh_spatial_preds.pdf"))
mean_pred
omega_pred
epsilon_pred
dev.off()


## INDEX ESTIMATES -------------------------------------------------------------

hss_pred2 <- exp_grid %>% 
  filter(
    week == "42",
    survey_f == "hss"
  )

index_preds <- purrr::map(
  fit_list, 
  ~ predict(.x, newdata = hss_pred2, return_tmb_object = TRUE)
)

index_preds2 <- purrr::map2(
  index_preds, names(mesh_list),
  ~ get_index(.x, bias_correct = TRUE) %>% 
    mutate(
      mesh = .y
    )
)
inds <- index_preds2 %>% 
  bind_rows() %>% 
  # remove early 90s estimates
  filter(!year < 1998)


pdf(here::here("figs", "mesh_comparison", "mesh_indices.pdf"))
ggplot(inds, aes(x = as.factor(year), y = est , fill = mesh)) + 
  geom_pointrange(aes(ymin = lwr , ymax = upr , fill = mesh),
                  position = position_dodge(width = 0.75), shape = 21) +
  ggsidekick::theme_sleek()
dev.off()


## COEFFICIENT ESTIMATES -------------------------------------------------------

fix_pars <- purrr::map2(
  fit_list, names(mesh_list),
  ~ tidy(.x, conf.int = T) %>% 
    mutate(
      mesh = .y
    )
) %>% 
  bind_rows() %>% 
  filter(!grepl("year", term))

ran_pars <- purrr::map2(
  fit_list, names(mesh_list),
  ~ tidy(.x, "ran_pars", conf.int = T) %>% 
    mutate(
      mesh = .y
    )
) %>% 
  bind_rows()


pdf(here::here("figs", "mesh_comparison", "mesh_par_ests.pdf"))
ggplot(fix_pars, aes(x = mesh, y = estimate, fill = mesh)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), shape = 21) +
  ggsidekick::theme_sleek() +
  facet_wrap(~term, scales = "free_y")
ggplot(ran_pars, aes(x = mesh, y = estimate, fill = mesh)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), shape = 21) +
  ggsidekick::theme_sleek() +
  facet_wrap(~term, scales = "free_y")
dev.off()



## COMPARE ANISOTROPY V. BARRIER MESH ------------------------------------------

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -137, ymin = 47, xmax = -121.25, ymax = 57) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))


pink_dat_trim <- pink_dat %>% 
  filter(!year < 1998)

# first mesh is same as finer res sdmTMB mesh above
spde <- make_mesh(pink_dat_trim,
                  c("utm_x_1000", "utm_y_1000"),
                  type = "kmeans",
                  n_knots = 250)
bspde <- add_barrier_mesh(spde, coast_utm, range_fraction = 0.1,
                          # scaling = 1000 since UTMs were rescaled above
                          proj_scaling = 1000)



## Compare anisotropy with finer resolution sdmTMB mesh
fit_list_ani <- purrr::map2(
  list(spde, bspde), list(TRUE, FALSE), 
  ~ sdmTMB(
    n_juv ~ 1 +
      as.factor(year) +
      dist_to_coast_km +
      s(week, bs = "cc", k = 5) +
      target_depth +
      day_night +
      survey_f,
    offset = pink_dat_trim$effort,
    data = pink_dat_trim,
    mesh = .x,
    family = sdmTMB::nbinom2(),
    spatial = "on",
    spatiotemporal = "iid",
    time = "year",
    anisotropy = .y,
    knots = list(
      week = c(0, 52)
    ),
    control = sdmTMBcontrol(
      newton_loops = 1
    )
  )
)

purrr::map(fit_list_ani, sanity)


# spatial predictions for anisotropy
pred_list_ani <- purrr::map2(
  fit_list_ani, c("anisotropy", "barrier"),
  ~ predict(.x, newdata = hss_pred, se_fit = FALSE, re_form = NULL) %>% 
    mutate(
      mesh = .y
    )
)
pred_dat_ani <- pred_list_ani %>% 
  bind_rows() %>% 
  mutate(
    exp_est = exp(est)
  )

pdf(here::here("figs", "mesh_comparison", "anisotropy_spatial_preds.pdf"))
plot_map(pred_dat_ani, "exp(est)") +
  scale_fill_viridis_c(
    trans = "sqrt",
    limits = c(0, quantile(exp(pred_dat_ani$est), 0.995))
  ) +
  facet_grid(year~mesh) +
  ggtitle("Prediction (fixed effects + all random effects)")
plot_map(pred_dat_ani %>% filter(year == "2012"), "omega_s") +
  scale_fill_gradient2() +
  facet_wrap(~mesh) +
  ggtitle("Prediction (spatial random effects)")
plot_map(pred_dat_ani, "epsilon_st") +
  scale_fill_gradient2() +
  facet_grid(year~mesh) +
  ggtitle("Prediction (spatiotemporal random effects)")
dev.off()


# index predictions
index_preds_ani <- purrr::map(
  fit_list_ani, 
  ~ predict(.x, 
            newdata = dum <- hss_pred2 %>% filter(year %in% pink_dat_trim$year),
            return_tmb_object = TRUE)
)

index_preds2_ani <- purrr::map2(
  index_preds_ani,  c("anisotropy", "barrier"),
  ~ get_index(.x, bias_correct = TRUE) %>% 
    mutate(
      mesh = .y
    )
)
inds_ani <- index_preds2_ani %>% 
  bind_rows() 


pdf(here::here("figs", "mesh_comparison", "anisotropy_indices.pdf"))
ggplot(inds_ani, aes(x = as.factor(year), y = est , fill = mesh)) + 
  geom_pointrange(aes(ymin = lwr , ymax = upr , fill = mesh),
                  position = position_dodge(width = 0.75), shape = 21) +
  ggsidekick::theme_sleek()
dev.off()
