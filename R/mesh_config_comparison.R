### Explore different mesh configurations for high seas juvenile data
## April 30, 2020

library(tidyverse)
library(sdmTMB)
library(ggplot2)


bridge <- readRDS(here::here("data", "ipes_hs_merged_bridge.rds")) 

# use bridge data alone for now since ignoring stock composition
jchin <- bridge %>%
  #scale UTM coords
  mutate(xUTM_start = xUTM_start / 10000,
         yUTM_start = yUTM_start / 10000,
         season = as.factor(
           case_when(
             month %in% c("2", "3") ~ "winter",
             month %in% c("5", "6", "7", "8") ~ "summer",
             month %in% c("9", "10", "11" , "12") ~ "fall")
         ),
         yday_z = as.vector(scale(yday)[,1]),
         yday_z2 = yday_z^2,
         time_f = as.factor(time_f),
         juv_cpue = ck_juv / dur,
         offset = dur
  ) %>% 
  #remove extra vars and stock ppn data
  dplyr::select(station_id, stable_station:date, time_f, yday, yday_z, yday_z2,
                offset, 
                month:dur, head_depth,
                ck_juv, juv_cpue)

#trim dataset
jchin1 <- jchin %>% 
  filter(stable_station == "1",
         month %in% c(6, 7),
         head_depth < 21,
         !juv_cpue > 500 #remove large values that may be skewing results
  )


## Build different meshes ------------------------------------------------------

#alternative mesh using non-convex boundary based on predictive grid (excludes
#land masses)
source(here::here("R", "prepPredictionGrid.R"))
grid <- jchin1 %>% 
  mutate(xUTM_start = xUTM_start * 10000, #scale up to match reality
         yUTM_start = yUTM_start * 10000) %>% 
  make_pred_grid()
grid <- grid / 10000

bnd <- INLA::inla.nonconvex.hull(as.matrix(pred_grid), convex = -0.03)
plot(bnd$loc)
mesh.loc <- SpatialPoints(as.matrix(cbind(jchin1$xUTM_start, 
                                          jchin1$yUTM_start)))
mesh <- INLA::inla.mesh.2d(loc=mesh.loc,
                     boundary=list(
                       bnd,
                       NULL),
                     max.edge=c(2.5, 4.5),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=1.51, ## Filter away adjacent points.
                     offset=c(0.1, 0.3)) ## Offset for extra boundaries, if needed.
plot(mesh)
jchin1_spde_nch <- make_spde(jchin1$xUTM_start, jchin1$yUTM_start, mesh = mesh)

#default mesh w// approximately same number of nodes
jchin1_spde <- make_spde(jchin1$xUTM_start, jchin1$yUTM_start, n_knots = 250)

#compare 
plot_spde(jchin1_spde)
plot_spde(jchin1_spde_nch)
plot(jchin1_spde$mesh, main = NA, edge.color = "grey60", asp = 1)
plot(jchin1_spde_nch$mesh, main = NA, edge.color = "grey60", asp = 1)


## Fit both models -------------------------------------------------------------

m <- sdmTMB(ck_juv ~ 0 + as.factor(year) + yday_z + yday_z2 + offset,
                data = jchin1,
                time = "year",
                spde = jchin1_spde,
                silent = FALSE,
                anisotropy = TRUE,
                include_spatial = TRUE,
                ar1_fields = FALSE,
                family = nbinom2(link = "log"))
m_nch <- sdmTMB(ck_juv ~ 0 + as.factor(year) + yday_z + yday_z2 + offset,
                data = jchin1,
                time = "year",
                spde = jchin1_spde_nch,
                silent = FALSE,
                anisotropy = TRUE,
                include_spatial = TRUE,
                ar1_fields = FALSE,
                family = nbinom2(link = "log"))

# examine residuals
check_res <- function(mod, dat) {
  dat$resid <- residuals(mod)
  #histogram
  hist(dat$resid)
  # qqplot
  qqnorm(dat$resid)
  abline(a = 0, b = 1, col = "red")
  # map of resids
  ggplot(dat, aes(xUTM_start, yUTM_start, col = resid)) + 
    scale_colour_gradient2() +
    geom_point() + facet_wrap(~year) + coord_fixed()
}

check_res(m, jchin1)
check_res(m_nch, jchin1)

# generate predictions using grid from above
pred_grid <- grid %>% 
  expand(nesting(X, Y), #make a predictive dataframe
         year = unique(jchin1$year)) %>%
  mutate(yday_z = 0,
         yday_z2 = 0,
         offset = mean(jchin1$dur))

# look at annual index predictions from each model
pred <- predict(m, newdata = pred_grid, return_tmb_object = TRUE)
pred_nch <- predict(m_nch, newdata = pred_grid, return_tmb_object = TRUE)

# Compare predictions from the two datasets
ind1 <- get_index(pred, bias_correct = FALSE) %>% 
  mutate(dataset = "def")
ind2 <- get_index(pred_nch, bias_correct = FALSE) %>% 
  mutate(dataset = "nch")  
inds <- rbind(ind1, ind2)

scale <- 2 * 2  # 2 x 2 km grid 
ggplot(inds, aes(year, est*scale)) + 
  geom_line(aes(group = 1)) +
  geom_ribbon(aes(ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  xlab('Year') + ylab('Abundance Estimate') +
  facet_wrap(~dataset)

# Look at output maps
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

# Predictions incorporating all fixed and random effects
# Generally the impacts of different meshes considered here were negligible
plot_map(pred$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")
plot_map(pred_nch$data, "exp(est)") +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

# Spatial random effects (temporally stable factors driving changes in abundance)
plot_map(pred$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()
plot_map(pred_nch$data, "omega_s") +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# Spatiotemporal random effects (dynamic drivers)
plot_map(pred$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()
plot_map(pred_nch$data, "epsilon_st") +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()
