library(ggplot2)
library(ggmap)
library(sp)
library(sf)
library(raster)
library(dplyr)
library(spData)

wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)
jchin <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  filter(stableStation == "Y")

# minLat <- min(floor(jchin$start_lat))
# maxLat <- max(floor(jchin$start_lat))
# minLong <- min(floor(jchin$start_long))
# maxLong <- max(floor(jchin$start_long))
## UTMs
minLat2 <- min(floor(jchin$yUTM_start))
maxLat2 <- max(floor(jchin$yUTM_start))
minLong2 <- min(floor(jchin$xUTM_start))
maxLong2 <- max(floor(jchin$xUTM_start))


library(sf)
library(fasterize)
library(raster)
# coast <- rnaturalearth::ne_coastline(110, returnclass = "sf")
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf"))
coastUTM <- st_transform(coast, crs = "+proj=utm +zone=9 +datum=WGS84")
cropR <- raster(extent(minLong2, maxLong2, minLat2, maxLat2),
                crs = "+proj=utm +zone=9 +datum=WGS84", res = 20000)
g <- fasterize(coastUTM, cropR)

## fast conversion pixel to polygons
p <- spex::polygonize(!is.na(g))
## layer is whether we are on land or not
plot(subset(p, !layer)$geometry)
plot(coastUTM$geometry, add = TRUE)

dum <- subset(p, !layer)$geometry %>% 
  st_transform(., crs = "+proj=utm +zone=9 +datum=WGS84") %>% 
  st_coordinates(.)
gridOut <- data.frame(X = dum[, "X"],
                      Y = dum[, "Y"])

saveRDS(gridOut, here::here("data", "spatialData", "trimmedSurveyGrid.rds"))


  
### TRASHY CODE ###

# coords <- wCan %>% 
#   select(yUTM = lat, xUTM = long)
# sp::coordinates(coords) <- c("xUTM", "yUTM")
# sp::proj4string(coords) <- sp::CRS("+proj=longlat +datum=WGS84")
# #SE Van Island extends into zone 10; not sure if that's an issue or not
# coords2 <- sp::spTransform(coords, sp::CRS("+proj=utm +zone=9 ellps=WGS84")) %>%
#   as(., "SpatialPoints")
# 
# wCanOut <- wCan %>% 
#   cbind(., coords2@coords) %>% 
#   mutate(xUTM = xUTM / 10000,
#          yUTM = yUTM / 10000)

xlims <- c(-131, -124)
ylims <- c(48, 53)

p0 <- ggplot(wCan) +
  geom_path(aes(x = long, y = lat, group = group), color = "black", 
            size = 0.25) +
  coord_map(projection = "mercator") + 
  scale_x_continuous(limits = xlims, expand = c(0, 0)) + 
  scale_y_continuous(limits = ylims, expand = c(0, 0)) + 
  labs(list(title = "", x = "Longitude", y = "Latitude"))

# Pull eez data 
# very slow
eezPath <- here::here("data", "spatialData", "eezShapeFile")
eez <- rgdal::readOGR(dsn = eezPath,
                      layer = tools::file_path_sans_ext("eez_v10.shp"))
eezCan <- eez[eez@data$Territory1 == "Canada", ]
eezCan@data$id <- rownames(eezCan@data)
eezCanWest <- fortify(eezCan) %>%
  inner_join(., eezCan@data, by = "id") %>%
  filter(id == "228",
         long > -132 & long < -120,
         lat < 55) %>%
  dplyr::select(long, lat, group, hole, piece)

coords <- eezCanWest %>%
  dplyr::select(yUTM = lat, xUTM = long)
sp::coordinates(coords) <- c("xUTM", "yUTM")
sp::proj4string(coords) <- sp::CRS("+proj=longlat +datum=WGS84")
#zone is only relevant for NW VI
coords2 <- sp::spTransform(coords, sp::CRS("+proj=utm +zone=9 ellps=WGS84")) %>%
  as(., "SpatialPoints")
eezCanWestOut <- cbind(eezCanWest, coords2@coords) %>%
  mutate(xUTM = xUTM / 10000,
         yUTM = yUTM / 10000)
  
grid <- raster(extent(eezCanWestOut))

# save cropped shape file
# saveRDS(eezCanWestOut ,
#        here::here("data", "spatialData", "canadianEEZ.rds"))
# eezCanWestOut <- readRDS(here::here("data", "spatialData", "canadianEEZ.rds"))

p0 + 
  geom_path(data = eezCanWestOut,
            aes(x = long, y = lat, group = group),
            colour = "blue", size = 0.75)

jchin <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  filter(stableStation == "Y") %>%
  mutate(xUTM_start = xUTM_start / 10000,
         yUTM_start = yUTM_start / 10000)




trimGrid <- sp::point.in.polygon(surv_grid,
                                 eez_verts)

sp::over(surv_grid,
        eez_verts)

surv_grid_trim <- surv_grid %>% 
  mutate(trimGrid = trimGrid) %>% 
  filter(trimGrid == 1) %>% 
  dplyr::select(-trimGrid)
# saveRDS(surv_grid_trim, 
#         here::here("data", "spatialData", "trimmedSurveyGrid.rds"))
# surv_grid_trim <- readRDS(here::here("data", "spatialData", 
#                    "trimmedSurveyGrid.rds"))

ggplot(surv_grid_trim) +
  geom_tile(aes(x = X, y = Y)) +
  coord_fixed()

ggplot(eezCanWestOut) +
  geom_path(aes(x = xUTM, y = yUTM, group = group)) +
  # geom_point(data = surv_grid, aes(x = X, y = Y)) +
  coord_fixed()


geom_path(data = eezCanWestOut,
          aes(x = long, y = lat, group = group),
          colour = "blue", size = 0.75)

# open polygon:
point.in.polygon(1:10,1:10,c(3,5,5,3),c(3,3,5,5))
# closed polygon:
point.in.polygon(1:10,rep(4,10),c(3,5,5,3,3),c(3,3,5,5,3))