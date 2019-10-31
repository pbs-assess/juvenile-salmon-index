library(ggplot2)
library(ggmap)
library(sp)
library(sf)
library(raster)
library(dplyr)
library(spData)

wCan <- map_data("world", region = "canada") %>%
  filter(long < -110)

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
eezPath <- here::here("data", "spatialData", "eezShapeFile")
# very slow
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

surv_grid <- expand.grid(
  X = seq(from = min(jchin$xUTM_start),
          to = max(jchin$xUTM_start),
          length.out = 50),
  Y = seq(from = min(jchin$yUTM_start),
          to = max(jchin$yUTM_start),
          length.out = 50))

trimGrid <- sp::point.in.polygon(surv_grid$X, surv_grid$Y,
                                 eezCanWestOut$xUTM, eezCanWestOut$yUTM)

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
  geom_point(data = surv_grid, aes(x = X, y = Y)) +
  coord_fixed()


geom_path(data = eezCanWestOut,
          aes(x = long, y = lat, group = group),
          colour = "blue", size = 0.75)