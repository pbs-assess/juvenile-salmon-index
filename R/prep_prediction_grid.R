## Make predictive grid for sdmTMB models
# 1) Prep bathymetric data based on chinTagging/R/prep_bathymetry 
# (ignore American data for now); use UTM to ensure equal spacing in grid
# 2) Back-convert to lat/lon to use coastdistance function
# 3) Combine and export
# Updated May 17, 2022


library(sf)
library(raster)
library(tidyverse)
library(ncdf4)
library(maptools)
library(rmapshaper)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)


dat_trim <- readRDS(here::here("data", "chin_catch_sbc.rds"))


# parallelize based on operating system (should speed up some of the spatial
# processing calculations)
library("parallel")
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library("doParallel")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  doMC::registerDoMC(ncores)
}


## GENERATE BATHYMETRY RASTER _-------------------------------------------------

## import depth data from netcdf file for BC and US PNW
big_bathy_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/BC NetCDF"
ncin <- nc_open(
  paste(big_bathy_path, "british_columbia_3_msl_2013.nc", sep = "/"))
ncin_us <- nc_open(
  paste(big_bathy_path, "usgsCeCrm8_703b_754f_94dc.nc", sep = "/"))

#specify lat/long
dep_list_bc <- list(lon = ncvar_get(ncin, "lon"),
                    lat = ncvar_get(ncin, "lat"),
                    dep = ncvar_get(ncin, "Band1"))
dep_list_us <- list(lon = ncvar_get(ncin_us, "longitude"),
                    lat = ncvar_get(ncin_us, "latitude"),
                    dep = ncvar_get(ncin_us, "topo"))

# function create dataframe
dep_dat_f <- function(x) {
  expand.grid(lon = x$lon, lat = x$lat) %>%
    as.matrix(.) %>%
    cbind(., depth = as.vector(x$dep)) %>%
    data.frame()
}

dep_dat_bc <- dep_dat_f(dep_list_bc)
#trim BC data
dep_dat_bc2 <- dep_dat_bc %>%
  filter(lat <= ceiling(max(dat_trim$mean_lat)),
         lon >= floor(min(dat_trim$mean_lon))#,
         # !depth < -1000
         ) %>% 
  mutate(depth = -1 * depth)
# dep_dat_us <- dep_dat_f(dep_list_us)
# #trim US data
# dep_dat_us2 <- dep_dat_us %>%
#   filter(lat <= 48, #min(dep_dat_bc$lat),
#          # lon >= floor(min(dat_trim$mean_lon)),
#          !is.na(depth),
#          !depth > 0.00#,
#          # !depth < -1000
#          ) %>% 
#   mutate(depth = -1 * depth)

# convert each to raster, downscale, and add terrain data to both
# then convert back to dataframes
bc_raster <- rasterFromXYZ(dep_dat_bc2, 
                           crs = sp::CRS("+proj=longlat +datum=WGS84"))
bc_raster_utm <- projectRaster(bc_raster,
                               crs = sp::CRS("+proj=utm +zone=9 +units=m"),
                               # convert to 1000 m resolution
                               res = 1000)

# # merge and add aspect/slope
# coast_raster <- merge(bc_raster, us_raster)
bc_raster_slope <- terrain(bc_raster_utm, opt = 'slope', unit = 'degrees',
                              neighbors = 8)
bc_raster_aspect <- terrain(bc_raster_utm, opt = 'aspect', unit = 'degrees',
                               neighbors = 8)
bc_raster_list <- list(depth = bc_raster_utm,
                          slope = bc_raster_slope,
                          aspect = bc_raster_aspect)
# downscaled version
low_bc_raster_list <- purrr::map(bc_raster_list,
                                 aggregate, 
                                 fac = 20, fun = mean)


# check plots
par(mfrow = c(1, 3))
purrr::map(bc_raster_list, plot)
purrr::map(low_bc_raster_list, plot)
purrr::map(low_bc_raster_list2, plot)


saveRDS(bc_raster_list,
        here::here("data", "spatial", "bathy_raster_utm_1000m.RDS"))
saveRDS(low_bc_raster_list,
        here::here("data", "spatial", "bathy_lowres_bc_rasters.RDS"))


## GENERATE GRID ---------------------------------------------------------------

# lower resolution raster list generated above
# low_bc_raster_list <- readRDS(
#   here::here("data", "spatial", "bathy_lowres_bc_rasters.RDS"))
low_bc_raster_list <- readRDS(
  here::here("data", "spatial", "bathy_raster_utm_1000m.RDS"))


# boundary box for receiver locations
min_lat <- min(dat_trim$mean_lat - 0.1)
max_lat <- max(dat_trim$mean_lat + 0.1)
min_lon <- min(dat_trim$mean_lon - 0.1)
max_lon <- max(dat_trim$mean_lon + 0.1)


bathy_low_sf_list <- purrr::map2(
  low_bc_raster_list,
  names(low_bc_raster_list),
  function (x, y) {
    as(x, 'SpatialPixelsDataFrame') %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("x", "y"),
               crs = sp::CRS("+proj=utm +zone=9 +units=m"))
    # rename(lon = x, lat = y) %>%
    #   st_as_sf(., coords = c("lon", "lat"),
    #            crs = sp::CRS("+proj=longlat +datum=WGS84"))
  }
)


# join depth and slope data
bathy_low_sf <- st_join(bathy_low_sf_list$depth, bathy_low_sf_list$slope) 

## convert to lat/lon for coast distance function 
bathy_low_sf_deg <- bathy_low_sf %>%
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

# coast sf for plotting and calculating distance to coastline 
# (has to be lat/lon for dist2Line)
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = min_lon, ymin = 48, xmax = max_lon, ymax = max_lat) 
  
# calculate distance to coastline
coast_dist <- geosphere::dist2Line(p = sf::st_coordinates(bathy_low_sf_deg),
                                   line = as(coast, 'Spatial'))


coast_utm <- coast %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))



# combine all data
bathy_dat <- data.frame(
  st_coordinates(bathy_low_sf[ , 1]),
  depth = bathy_low_sf$depth,
  slope = bathy_low_sf$slope,
  shore_dist = coast_dist[, "distance"]
)

ggplot() + 
  geom_sf(data = coast_utm) +
  geom_raster(data = bathy_dat, aes(x = X, y = Y, fill = shore_dist)) +
  scale_fill_viridis_c() +
  geom_point(data = dat_trim, aes(x = utm_x, y = utm_y), 
             fill = "white", shape = 21) +
  ggsidekick::theme_sleek()


# export grid
saveRDS(bathy_dat, here::here("data", "spatial", "pred_bathy_grid.RDS"))


