## Make predictive grid for sdmTMB models
# 1) Prep bathymetric data based on chinTagging/R/prep_bathymetry 
# (ignore American data for now)
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
bc_raster_utm <- projectRaster()

# # merge and add aspect/slope
# coast_raster <- merge(bc_raster, us_raster)
bc_raster_slope <- terrain(bc_raster, opt = 'slope', unit = 'degrees',
                              neighbors = 8)
bc_raster_aspect <- terrain(bc_raster, opt = 'aspect', unit = 'degrees',
                               neighbors = 8)
bc_raster_list <- list(depth = bc_raster,
                          slope = bc_raster_slope,
                          aspect = bc_raster_aspect)
# downscaled version
low_bc_raster_list <- purrr::map(bc_raster_list,
                                 aggregate, 
                                 fac = 20, fun = mean)
low_bc_raster_list2 <- purrr::map(bc_raster_list,
                                 aggregate, 
                                 fac = 40, fun = mean)

# check plots
par(mfrow = c(1, 3))
purrr::map(bc_raster_list, plot)
purrr::map(low_bc_raster_list, plot)
purrr::map(low_bc_raster_list2, plot)


# saveRDS(low_coast_raster_list,
#         here::here("data", "spatial", "bathy_lowres_rasters.RDS"))
saveRDS(low_bc_raster_list,
        here::here("data", "spatial", "bathy_lowres_bc_rasters.RDS"))


## GENERATE GRID ---------------------------------------------------------------

# lower resolution raster list generated above
low_bc_raster_list <- readRDS(
  here::here("data", "spatial", "bathy_lowres_bc_rasters.RDS"))


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
      rename(lon = x, lat = y) %>%
      st_as_sf(., coords = c("lon", "lat"),
               crs = sp::CRS("+proj=longlat +datum=WGS84")) #%>%
    # convert to utm
    # st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m")) 
  }
)

# join depth and slope data
bathy_low_sf <- st_join(bathy_low_sf_list$depth, bathy_low_sf_list$slope) 

# coast sf for plotting and calculating distance to coastline 
# (has to be lat/lon for dist2Line)
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = min_lon, ymin = 48, xmax = max_lon, ymax = max_lat)
# saveRDS(coast, 
#         here::here("data", "spatial_data",
#                    "crop_coast_sf.RDS"))

coast_dist <- geosphere::dist2Line(p = sf::st_coordinates(bathy_low_sf), 
                                   line = as(coast, 'Spatial'))

bathy_dat <- data.frame(
  st_coordinates(bathy_low_sf[ , 1]),
  depth = bathy_low_sf$depth,
  slope = bathy_low_sf$slope,
  shore_dist = coast_dist[, "distance"]
)


ggplot() + 
  geom_sf(data = coast) +
  geom_raster(data = bathy_dat, aes(x = X, y = Y, fill = shore_dist)) +
  scale_fill_viridis_c() +
  geom_point(data = dat_trim, aes(x = mean_lon, y = mean_lat), 
             fill = "white", shape = 21) +
  ggsidekick::theme_sleek()


saveRDS(bathy_dat, here::here("data", "spatial", "pred_bathy_grid.RDS"))










## as above but with UTM 
bathy_low_sf_utm <- bathy_low_sf %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m")) 

# coast sf for plotting and calculating distance to coastline
bathy_dat_utm <- data.frame(
  st_coordinates(bathy_low_sf_utm[ , 1]),
  depth = bathy_low_sf_utm$layer,
  slope = bathy_low_sf_utm$slope,
  shore_dist = coast_dist[, "distance"]
) 

coast_utm <- coast %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=10 +units=m")) 

ggplot() + 
  geom_sf(data = coast_utm) +
  geom_point(data = bathy_dat_utm, aes(x = X, y = Y), size = 0.05) +
  # geom_tile(data = bathy_dat_utm, aes(x = X, y = Y, fill = shore_dist), 
  #           width = 1500, height = 1500) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek()





# wCan <- map_data("world", region = "canada") %>%
#   filter(long < -110)



make_pred_grid <- function(dat, file_name_out = "trimmedSurveyGrid.rds") {
  ## Select boundary box based on survey area
  minLat2 <- min(floor(dat$utm_y))
  maxLat2 <- max(floor(dat$utm_y))
  minLong2 <- min(floor(dat$utm_x))
  maxLong2 <- max(floor(dat$utm_x))
  
  # First generate grid for entire area that just excludes landmass
  projCRS <- "+proj=utm +zone=10 +datum=WGS84"
  coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                           returnclass = "sf"), 
                 rnaturalearth::ne_states( "Canada", returnclass = "sf"))
  coastUTM <- st_transform(coast, crs = projCRS)
  cropR <- raster(extent(minLong2, maxLong2, minLat2, maxLat2),
                  crs = projCRS, res = 10000)
  g <- fasterize(coastUTM, cropR)
  
  ## fast conversion pixel to polygons
  p <- spex::polygonize(!is.na(g))
  ## layer is whether we are on land or not
  plot(subset(p, !layer)$geometry)
  plot(coastUTM$geometry, add = TRUE)
  
  
  # Now subset grid based on whether the cells have contained a set in recent years
  # fullGrid <- subset(p, !layer)$geometry %>% 
  #   st_sf(ID = seq(1, length(.), by = 1))
  # sets <- dat %>% 
  #   #remove early years that include a large number of offshore and northern sets
  #   dplyr::select(year, utm_x, utm_y) %>% 
  #   st_as_sf(., coords = c("utm_x", "utm_y"), crs = projCRS)
  # 
  # # Subset spatially
  # subGrid <- fullGrid %>%
  #   st_join(., sets, join = st_intersects) %>% 
  #   filter(!is.na(year)) %>% 
  #   distinct()
  # plot(st_geometry(subGrid))
  # 
  # Extract just the coordinates and export
  # subGridCoords <- subGrid %>%
  #   st_coordinates(.)
  
  grid_coords <- st_coordinates(p)
  
  gridOut <- data.frame(X = grid_coords[, "X"],
                        Y = grid_coords[, "Y"])
  return(gridOut)
}

