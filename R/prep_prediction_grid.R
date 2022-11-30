## Make predictive grid for sdmTMB models
# 1) Prep bathymetric data based on chinTagging/R/prep_bathymetry 
# (ignore American data for now); use UTM to ensure equal spacing in grid
# 2) Back-convert to lat/lon to use coastdistance function
# 3) Combine and export
# Updated May 17, 2022


library(sf)
library(raster)
library(rgdal)
library(tidyverse)
library(ncdf4)
# library(maptools)
library(rmapshaper)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)


# chinook dataset
dat_trim <- readRDS(here::here("data", "chin_catch_sbc.rds"))

# shapefiles for IPES survey grid and combined WCVI/IPES grid
ipes_grid_raw <- readOGR(
  here::here("data", "spatial", "ipes_shapefiles", "IPES_Grid_UTM9.shp"))
ipes_wcvi_grid_raw <- readOGR(
  here::here("data", "spatial", "wcvi_ipes_shapefiles", 
             "IPES_WCVI_boundary_UTM9.shp"))


# parallelize based on operating system (should speed up some of the spatial
# processing calculations)
library("parallel")
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library("doParallel")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  big_bathy_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/BC NetCDF"
  creel_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/creel_areas/"
} else {
  doMC::registerDoMC(ncores)
  big_bathy_path <- "/Users/cam/Google Drive/spatial/BC NetCDF"
  creel_path <- "/Users/cam/Google Drive/spatial/creel_areas/"
}



## GENERATE BATHYMETRY RASTER _-------------------------------------------------

## import depth data from netcdf file for BC and US PNW
## UPDATE: switch to single integrated lower res file
# ncin <- nc_open(
#   paste(big_bathy_path, "british_columbia_3_msl_2013.nc", sep = "/"))
# ncin_us <- nc_open(
#   paste(big_bathy_path, "usgsCeCrm8_703b_754f_94dc.nc", sep = "/"))
ncin <- nc_open(
  paste(
    big_bathy_path, 
    "gebco_2022_n58.9197_s47.033_w-137.9622_e-121.9747_ipes.nc",
    sep = "/"
  )
)

#specify lat/long
dep_list <- list(
  lon = ncvar_get(ncin, "lon"),
  lat = ncvar_get(ncin, "lat"),
  dep = ncvar_get(ncin, "elevation")
)

# function create dataframe
dep_dat_f <- function(x) {
  expand.grid(lon = x$lon, lat = x$lat) %>%
    as.matrix(.) %>%
    cbind(., depth = as.vector(x$dep)) %>%
    data.frame()
}

dep_dat_full <- dep_dat_f(dep_list) %>% 
  mutate(depth = -1 * depth) %>% 
  # remove land data
  filter(depth > 0)


# convert each to raster, downscale, and add terrain data to both
# then convert back to dataframes
bc_raster <- rasterFromXYZ(
  dep_dat_full, 
  crs = sp::CRS("+proj=longlat +datum=WGS84")
)
bc_raster_utm <- projectRaster(
  bc_raster,
  crs = sp::CRS("+proj=utm +zone=9 +units=m"),
  # convert to 1000 m resolution
  res = 1000
)

# save RDS for manuscript figs
saveRDS(bc_raster, 
        here::here("data", "spatial", "full_coast_raster_latlon_1000m.RDS"))


# crop to survey grid
dum <- crop(bc_raster_utm, extent(ipes_grid_raw))
ipes_raster_utm <- mask(dum, ipes_grid_raw)
wcvi_ipes_raster_utm <- mask(dum, ipes_wcvi_grid_raw)


# # merge and add aspect/slope
wcvi_ipes_raster_slope <- terrain(
  wcvi_ipes_raster_utm, opt = 'slope', unit = 'degrees',
  neighbors = 8
)
wcvi_ipes_raster_aspect <- terrain(
  wcvi_ipes_raster_utm, opt = 'aspect', unit = 'degrees',
  neighbors = 8
)
wcvi_ipes_raster_list <- list(
  depth = wcvi_ipes_raster_utm,
  slope = wcvi_ipes_raster_slope,
  aspect = wcvi_ipes_raster_aspect
)
ipes_raster_slope <- terrain(
  ipes_raster_utm, opt = 'slope', unit = 'degrees',
  neighbors = 8
)
ipes_raster_aspect <- terrain(
  ipes_raster_utm, opt = 'aspect', unit = 'degrees',
  neighbors = 8
)
ipes_raster_list <- list(
  depth = ipes_raster_utm,
  slope = ipes_raster_slope,
  aspect = ipes_raster_aspect
)

# combine ipes and wcvi lists for subsequent cleaning
raster_list <- c(wcvi_ipes_raster_list, ipes_raster_list)
names(raster_list) <- c(
  paste("wcvi", names(wcvi_ipes_raster_list), sep = "_"),
  paste("ipes", names(ipes_raster_list), sep = "_")
)

saveRDS(bc_raster_utm,
        here::here("data", "spatial", "coast_raster_utm_1000m.RDS"))
saveRDS(raster_list,
        here::here("data", "spatial", "wcvi_ipes_raster_utm_1000m.RDS"))



## GENERATE GRID ---------------------------------------------------------------

# raster list generated above
raster_list <- readRDS(
        here::here("data", "spatial", "wcvi_ipes_raster_utm_1000m.RDS"))

sf_list <- purrr::map2(
  raster_list,
  names(raster_list),
  function (x, y) {
    as(x, 'SpatialPixelsDataFrame') %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("x", "y"),
               crs = sp::CRS("+proj=utm +zone=9 +units=m"))
  }
)


# join depth and slope data
wcvi_sf <- st_join(sf_list$wcvi_depth, sf_list$wcvi_slope) 
ipes_sf <- st_join(sf_list$ipes_depth, sf_list$ipes_slope) 


# coast sf for plotting and calculating distance to coastline 
# (has to be lat/lon for dist2Line)
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = min_lon, ymin = 48, xmax = max_lon, ymax = max_lat) 


## convert to lat/lon for coast distance function 
wcvi_sf_deg <- wcvi_sf %>%
  st_transform(., crs = st_crs(coast))

# calculate distance to coastline for WCVI SF object
coast_dist <- geosphere::dist2Line(
  p = sf::st_coordinates(wcvi_sf_deg),
  line = as(coast, 'Spatial')
)

# combine all data
wcvi_grid <- data.frame(
  st_coordinates(wcvi_sf[ , 1]),
  depth = wcvi_sf$depth,
  slope = wcvi_sf$slope,
  shore_dist = coast_dist[, "distance"]
)
# use coast_dist of wcvi_grid for ipes
ipes_grid <- data.frame(
  st_coordinates(ipes_sf[ , 1]),
  depth = ipes_sf$depth,
  slope = ipes_sf$slope
) %>% 
  left_join(., 
            wcvi_grid %>% select(X, Y, shore_dist),
            by = c("X", "Y")) 

# interpolate missing data in both 
wcvi_grid_interp <- VIM::kNN(wcvi_grid, k = 5)
ipes_grid_interp <- VIM::kNN(ipes_grid, k = 5)


coast_utm <- coast %>% 
  sf::st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

ggplot() + 
  geom_sf(data = coast_utm) +
  geom_sf(data = ipes_sf, aes(fill = depth)) +
  # geom_raster(data = wcvi_grid_interp, aes(x = X, y = Y, fill = depth)) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek()


# export grids as list
grid_list <- list(
  wcvi_grid = wcvi_grid_interp %>% 
    select(-ends_with("imp")),
  ipes_grid = ipes_grid_interp %>% 
    select(-ends_with("imp"))
)

saveRDS(grid_list,
        here::here("data", "spatial", "pred_ipes_grid.RDS"))
saveRDS(coast_utm, here::here("data", "spatial", "coast_trim_utm.RDS"))


