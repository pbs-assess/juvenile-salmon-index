## Clean BCSI Data
# Uses .csv files produced by Amy Tabata May 2022 (updated Oct 2022)

library(tidyverse)
library(sf)
library(sp)

# original depths data generated separately but has been merged with bridge data
# depths <- read.csv(
#   here::here("data", "event_depths_from_bathymetry_hires.csv")
# ) %>% 
#   janitor::clean_names()
chin <- read.csv(
  here::here("data", "BCSI_Juv_CHINOOK_Counts_True20220504.csv")
) %>% 
  janitor::clean_names()
bridge_raw <- read.csv(here::here("data", "BCSI_Bridge_Info_20220929.csv")) %>% 
  janitor::clean_names() 
synoptic_stations <- read.csv(here::here("data", "synoptic_stations.csv")) %>% 
  janitor::clean_names()


## INTERSECTION WITH IPES GRID -------------------------------------------------

## import and transform IPES grid
# ipes_sf_deg <- readRDS(here::here("data", "spatial", "ipes_sf_list_deg.RDS"))
ipes_grid_raw <- rgdal::readOGR(
  here::here("data", "spatial", "ipes_shapefiles", "IPES_Grid_UTM9.shp"))

# convert shapefile to sf, then points, then combine into single polygon
ipes_sf <- as(ipes_grid_raw, 'SpatialPolygonsDataFrame')
ipes_sf_pts <- st_as_sf(ipes_grid_raw, as_points = TRUE, merge = TRUE)
ipes_sf_poly <- st_union(ipes_sf_pts, by_feature = FALSE)


## add UTM to bridge and convert to sf
get_utm <- function(x, y, zone, loc){
  points = SpatialPoints(cbind(x, y),
                         proj4string = CRS("+proj=longlat +datum=WGS84"))
  points_utm = spTransform(
    points, CRS(paste0("+proj=utm +zone=",zone[1]," +units=m"))
  )
  if (loc == "x") {
    return(coordinates(points_utm)[,1])
  } else if (loc == "y") {
    return(coordinates(points_utm)[,2])
  }
}

bridge <- bridge_raw %>% 
  mutate(mean_lat = ifelse(is.na(end_latitude),
                           start_latitude,
                           (start_latitude + end_latitude) / 2),
         mean_lon = ifelse(is.na(end_longitude),
                           start_longitude,
                           (start_longitude + end_longitude) / 2),
         utm_x = get_utm(mean_lon, mean_lat, zone = "9", loc = "x"),
         utm_y = get_utm(mean_lon, mean_lat, zone = "9", loc = "y")) 

bridge_sf <- bridge %>% 
  select(unique_event, mean_lat, mean_lon) %>% 
  st_as_sf(., coords = c("mean_lon", "mean_lat"), 
           crs = sp::CRS("+proj=longlat +datum=WGS84")) 

bridge_sf2 <- bridge_sf %>% 
  st_transform(., crs = st_crs(ipes_sf_poly))
    
ggplot() +
  geom_sf(data = st_intersection(bridge_sf2, ipes_sf_poly))

# extract events within ipes survey grid
ipes_grid_events <- st_intersection(bridge_sf2, ipes_sf_poly) %>% 
  pull(., unique_event)


## CLEAN AND MERGE -------------------------------------------------------------

# impute missing bridge data
imp_dat <- bridge %>% 
  select(utm_x, utm_y, trip_id, tow_length_hour, mouth_height_use:volume_km3) 
imp_dat2 <- VIM::kNN(imp_dat, k = 5)
imp_dat2 <- imp_dat2 %>% 
  select(-ends_with("imp")) %>% 
  mutate(unique_event = bridge$unique_event)
  
dat <- bridge %>% 
  # remove values and replace with imputed data above
  select(-c(utm_x, utm_y, trip_id, tow_length_hour, 
           mouth_height_use:volume_km3)) %>% 
  left_join(., imp_dat2, by = "unique_event") %>% 
  left_join(., chin, by = "unique_event") %>%
  mutate(
    date = as.POSIXct(date,
                           format = "%Y-%m-%d",
                           tz = "America/Los_Angeles"),
    month = lubridate::month(date),
    week = lubridate::week(date),
    season = case_when(
      month %in% c("2", "3", "4") ~ "sp",
      month %in% c("5", "6", "7", "8") ~ "su",
      month %in% c("9", "10", "11", "12") ~ "wi"
    ),
    season_f = as.factor(season),
    # define core area as southern BC and SEAK excluding SoG
    synoptic_station = ifelse(
      mean_lat > 47 & mean_lat < 56 & !grepl("GS", station_name) & 
        mean_lon > -135,
      TRUE,
      FALSE
      ),
    # define IPES based on intersections with grid
    ipes_grid = ifelse(
      unique_event %in% ipes_grid_events,
      TRUE,
      FALSE
    ),
    survey_f = ifelse(
      year > 2016 & season_f == "su", "ipes", "hss") %>% as.factor(),
    year_f = as.factor(year),
    vessel_name = tolower(vessel_name),
    vessel = case_when(
      grepl("crest", vessel_name) ~ "sea crest",
      grepl("franklin", vessel_name) ~ "franklin",
      grepl("ricker", vessel_name) ~ "ricker",
      grepl("viking", vessel_name) ~ "viking storm",
      grepl("nordic pearl", vessel_name) ~ "nordic pearl",
      grepl("anita", vessel_name) ~ "anita j",
      grepl("frosti", vessel_name) ~ "frosti",
      grepl("columbia", vessel_name) ~ "columbia",
      grepl("ocean selector", vessel_name) ~ "ocean selector"
    )
  ) %>% 
  dplyr::select(unique_event, date, year, month, week, day, day_night, season_f,
         stratum:station_name, pfma = dfo_stat_area_code,
         synoptic_station, ipes_grid, survey_f,
         mean_lat, mean_lon, utm_x, utm_y,
         vessel, distance_travelled, vessel_speed, 
         bath_depth_mean_m, dist_to_coast_km, 
         mouth_height = mouth_height_use, mouth_width = mouth_width_use, 
         volume_km3, ck_juv = n_juv, ck_ad = n_ad) 


# subset to core area and remove rows with missing data
dat_trim <- dat %>%
  filter(synoptic_station == TRUE,
         # exclude SoG and Puget Sound PFMAs
         !pfma %in% c("13", "14", "15", "16", "17", "18", "19", "PS")) %>%
  droplevels()


min_lat <- min(floor(dat_trim$mean_lat) - 0.1)
max_lat <- max(dat_trim$mean_lat) + 0.1
min_lon <- min(floor(dat_trim$mean_lon) - 0.1)
max_lon <- max(dat_trim$mean_lon) + 0.1

coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))
  

set_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = dat_trim,
             aes(x = utm_x, y = utm_y, fill = ipes_grid), 
             shape = 21, alpha = 0.4) +
  facet_wrap(~season_f) +
  ggsidekick::theme_sleek() +
  theme(#legend.position = "top",
        axis.title = element_blank())

png(here::here("figs", "set_map.png"), height = 4, width = 8, res = 250, 
    units = "in")
set_map
dev.off()


# export subsetted version to use for initial fitting
saveRDS(dat_trim, here::here("data", "chin_catch_sbc.rds"))
saveRDS(coast, here::here("data", "spatial", "sbc_sf_utm.rds"))


## DATA EXPLORE ----------------------------------------------------------------

# some potential issues with limited data for certain vessels (don't assume
# models will perform well)
vessel_cov <- dat_trim %>% 
  group_by(year, month, vessel) %>% 
  summarize(n_tows = length(unique_event)) %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(x = month, y = year, size = n_tows),
              alpha = 0.3, width = 0.25) +
  facet_wrap(~vessel) +
  ggsidekick::theme_sleek()


# check seasonal coverage
temp_cov <- dat_trim %>% 
  group_by(year, month, survey_f) %>% 
  summarize(n_tows = length(unique_event)) %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(x = month, y = year, size = n_tows, fill = survey_f),
              alpha = 0.3, width = 0.25, shape = 21) +
  ggsidekick::theme_sleek()


dat_trim$log_dist_coast <- log(dat_trim$dist_to_coast_km)
dat_trim$log_depth <- log(dat_trim$depth_mean_m)

plot_foo <- function(x_in) {
  ggplot(data = dat_trim) +
    geom_point(aes_string(x = x_in, y = "ck_juv")) +
    ggsidekick::theme_sleek() +
    facet_wrap(~season_f, scales = "free")
}

catch_eff <- plot_foo("distance_travelled")
catch_dist <- plot_foo("dist_to_coast_km")
catch_bathy <- plot_foo("depth_mean_m")


pdf(here::here("figs", "exp_figs.pdf"))
vessel_cov
temp_cov
catch_eff
catch_dist
catch_bathy
dev.off()


png(here::here("figs", "temporal_coverage.png"), height = 4, width = 4, 
    res = 250, 
    units = "in")
temp_cov
dev.off()


# VISUALIZE DIFFERENT DATA COVERAGE -------------------------------------------- 

base_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  facet_wrap(~season_f, drop = FALSE) +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text = element_blank())

map1 <- base_map +
  geom_point(data = dat_trim,
             aes(x = utm_x, y = utm_y, fill = ipes_grid), 
             shape = 21, alpha = 0.4)
map2 <- base_map +
  geom_point(data = dat_trim %>% filter(ipes_grid == "TRUE"),
             aes(x = utm_x, y = utm_y, fill = ipes_grid), 
             shape = 21, alpha = 0.4)
map3 <- base_map +
  geom_point(data = dat_trim %>% filter(season_f == "su"),
             aes(x = utm_x, y = utm_y, fill = ipes_grid), 
             shape = 21, alpha = 0.4)
map4 <- base_map +
  geom_point(data = dat_trim %>% filter(season_f == "su",
                                        ipes_grid == "TRUE"),
             aes(x = utm_x, y = utm_y, fill = ipes_grid), 
             shape = 21, alpha = 0.4)

png(here::here("figs", "data_input_map.png"), height = 5, width = 8, res = 250, 
    units = "in")
cowplot::plot_grid(map1, map2, map3, map4, nrow = 2)
dev.off()