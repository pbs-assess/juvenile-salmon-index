library(ggplot2)
library(ggmap)
library(sp)
library(sf)
library(raster)
library(dplyr)
library(spData)

library(tidyverse)
library(sf)
library(fasterize)
library(raster)
library(spData)

wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)
jchin <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  filter(stableStation == "Y")

## Select boundary box based on survey area
minLat2 <- min(floor(jchin$yUTM_start))
maxLat2 <- max(floor(jchin$yUTM_start))
minLong2 <- min(floor(jchin$xUTM_start))
maxLong2 <- max(floor(jchin$xUTM_start))


# coast <- rnaturalearth::ne_coastline(110, returnclass = "sf")
projCRS <- "+proj=utm +zone=9 +datum=WGS84"
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf"))
coastUTM <- st_transform(coast, crs = projCRS)
cropR <- raster(extent(minLong2, maxLong2, minLat2, maxLat2),
                crs = projCRS, res = 20000)
g <- fasterize(coastUTM, cropR)

## fast conversion pixel to polygons
p <- spex::polygonize(!is.na(g))
## layer is whether we are on land or not
plot(subset(p, !layer)$geometry)
plot(coastUTM$geometry, add = TRUE)
plot(jchin)

fullGrid <- subset(p, !layer)$geometry %>% 
  st_sf(ID = seq(1, length(.), by = 1))
sets <- jchin %>% 
  #remove early years that include a relatively large number of offshore and N
  # sets
  filter(!year < 2013) %>%
  dplyr::select(year, xUTM_start, yUTM_start) %>% 
  st_as_sf(., coords = c("xUTM_start", "yUTM_start"), crs = projCRS)

subGrid <- fullGrid %>%
  st_join(., sets, join = st_intersects) %>% 
  filter(!is.na(year)) 
plot(st_geometry(subGrid))

subGridCoords <- subGrid %>%
  st_coordinates(.)
gridOut <- data.frame(X = subGridCoords[, "X"],
                      Y = subGridCoords[, "Y"])
saveRDS(gridOut, here::here("data", "trimmedSurveyGrid.rds"))


#alternative method doesn't filter by year and drops cells based on whether 
#x number of years occur in bounds, but results in gaps 
#subGrid <- fullGrid %>%
#   st_join(., sets, join = st_intersects) %>% 
#   filter(!is.na(year)) %>% 
#   group_by(ID, year) %>% 
#   tally() %>% 
#   filter(!n < 6)
