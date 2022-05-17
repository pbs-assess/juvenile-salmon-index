## Make predictive grid for sdmTMB models
# Goal is to constrain the grid to cells with a lot of data in a relatively 
# constrained space. Begin focusing on WCVI, QCS and JS, excluding grid cells
# that are a) on land or b) infrequently surveyed 
# Nov. 6 2019

library(tidyverse)
library(sf)
library(fasterize)
library(raster)
library(spData)

wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)

dat_trim <- readRDS(here::here("data", "chin_catch_sbc.rds"))
coast_trim <- readRDS(here::here("data", "sbc_sf_utm.rds"))


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

