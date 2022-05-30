## Clean BCSI Data
# Uses .csv files produced by Amy Tabata May 2022


library(tidyverse)
library(sf)
library(sp)

depths <- read.csv(
  here::here("data", "event_depths_from_bathymetry_hires.csv")
) %>% 
  janitor::clean_names()
chin <- read.csv(
  here::here("data", "BCSI_Juv_CHINOOK_Counts_True20220504.csv")
) %>% 
  janitor::clean_names()
bridge <- read.csv(here::here("data", "BCSI_Bridge_Info_20220505.csv")) %>% 
  janitor::clean_names() 
synoptic_stations <- read.csv(here::here("data", "synoptic_stations.csv")) %>% 
  janitor::clean_names()

dat <- bridge %>% 
  left_join(., 
            depths %>% dplyr::select(unique_event, depth_mean_m, dist_to_coast_km,
                              bridge_dist_km),
            by = "unique_event") %>% 
  left_join(., chin, by = "unique_event") %>% 
  mutate(
    date = as.POSIXct(date,
                           format = "%Y-%m-%d",
                           tz = "America/Los_Angeles"),
    month = lubridate::month(date),
    mean_lat = ifelse(is.na(end_latitude),
                      start_latitude,
                      (start_latitude + end_latitude) / 2),
    mean_lon = ifelse(is.na(end_longitude),
                      start_longitude,
                      (start_longitude + end_longitude) / 2),
    # incomplete flagging use lat for now
    # synoptic_station = ifelse(unique_event %in% synoptic_stations$station_id,
    #                           TRUE,
    #                           FALSE),
    synoptic_station = ifelse(
      mean_lat > 48 & mean_lat < 52 & !grepl("GS", station_name) & 
        mean_lon > -131,
      TRUE,
      FALSE
      ),
    season = case_when(
      month %in% c("4", "5", "6") ~ "sp",
      month %in% c("7", "8", "9") ~ "su",
      month %in% c("10", "11", "12") ~ "fa",
      month %in% c("1", "2", "3") ~ "wi"
    ),
    # missing bathymetry for a few cases
    depth_mean_m = ifelse(!is.na(depth_mean_m),
                          as.numeric(depth_mean_m),
                          start_bottom_depth),
    # bridge_dist_km = as.numeric(bridge_dist_km),
    dist_to_coast_km = as.numeric(dist_to_coast_km),
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
  dplyr::select(unique_event, date, year, month, day, day_night, season,
         stratum:station_name, synoptic_station, mean_lat, mean_lon, 
         vessel, distance_travelled, vessel_speed, 
         depth_mean_m, dist_to_coast_km, 
         mouth_height = mouth_height_use, mouth_width = mouth_width_use,
         ck_juv = n_juv, ck_ad = n_ad)


## PREP SPATIAL ----------------------------------------------------------------

# NOTE full dataset spans many UTM zones so only apply to southern regions
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
 
dat_trim <- dat %>%
  filter(synoptic_station == TRUE) %>% 
  mutate(utm_x = get_utm(mean_lon, mean_lat, "9", loc = "x"),
         utm_y = get_utm(mean_lon, mean_lat, "9", loc = "y"))


min_lat <- min(floor(dat_trim$mean_lat) - 0.1)
max_lat <- max(dat_trim$mean_lat) + 0.1
min_lon <- min(floor(dat_trim$mean_lon) - 0.1)
max_lon <- max(dat_trim$mean_lon) + 0.1

coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))
  

ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = dat_trim,
             aes(x = utm_x, y = utm_y, fill = year), 
             shape = 21, alpha = 0.4)


# some potential issues with limited data for certain vessels (don't assume
# models will perform well)
dat_trim %>% 
  group_by(year, month, vessel) %>% 
  summarize(n_tows = length(unique_event)) %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(x = month, y = year, size = n_tows),
              alpha = 0.3, width = 0.25) +
  facet_wrap(~vessel) +
  ggsidekick::theme_sleek()



# export subsetted version to use for initial fitting
saveRDS(dat_trim, here::here("data", "chin_catch_sbc.rds"))
saveRDS(coast, here::here("data", "spatial", "sbc_sf_utm.rds"))
