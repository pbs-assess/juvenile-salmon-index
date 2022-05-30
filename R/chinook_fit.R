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
         year_f = as.factor(year)) %>% 
  filter(!is.na(depth_mean_m),
         !is.na(dist_to_coast_km))

grid <- readRDS(here::here("data", "spatial", "pred_bathy_grid.RDS")) %>% 
  mutate(utm_x_1000 = X / 1000,
         utm_y_1000 = Y / 1000)

coast_utm <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -132, ymin = 47, xmax = -121.25, ymax = 52.5) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))



## MAKE MESH -------------------------------------------------------------------

# have to account for landmasses so use predictive grid as baseline 
# TODO: generate grid without SOG based on BSCI shape file

# spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
#                            n_knots = 250, type = "kmeans")
spde <- make_mesh(dat_trim, c("utm_x_1000", "utm_y_1000"), 
                   cutoff = 5, type = "kmeans")
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
  geom_sf(data = mesh_df_land, size = 1, colour = "green")


#compare 
plot(dat_trim_spde)
plot(jchin1_spde_nch)
plot(jchin1_spde$mesh, main = NA, edge.color = "grey60", asp = 1)
plot(jchin1_spde_nch$mesh, main = NA, edge.color = "grey60", asp = 1)


## FIT SIMPLE MODEL ------------------------------------------------------------

# Include covariates for monthly smooth, year factor, distance travelled, 
# bathymetry and distance to coast
fit <- sdmTMB(ck_juv ~ s(depth_mean_m) + s(dist_to_coast_km) + 
                s(month),
              data = dat_trim,
              spde = bspde,
              # silent = FALSE,
              # anisotropy = TRUE,
              # ar1_fields = FALSE,
              family = nbinom2(link = "log"))


