## Clean BCSI Data
# Uses .csv files produced by Amy Tabata May 2022 (updated Feb 2023)

library(tidyverse)
library(sf)
library(sp)


bridge_raw <- read.csv(here::here("data", "BCSI_TowInfo_2023215.csv")) %>% 
  #read.csv(here::here("data", "BCSI_Bridge_Info_20220929.csv")) %>% 
  janitor::clean_names()  %>%
  # exclude some absurdly deep tows
  filter(!target_depth > 75) %>% 
  mutate(
    # missing day/night classifier for 2022 fall survey 
    day_night = ifelse(trip_id == "2022-011", "DAY", day_night),
    # calculate mean depth 
    mean_bottom_depth = ifelse(
      !is.na(end_bottom_depth),
      (start_bottom_depth + end_bottom_depth) / 2,
      start_bottom_depth
    )
  )

synoptic_stations <- read.csv(here::here("data", "synoptic_stations.csv")) %>% 
  janitor::clean_names()

# target headrope depths generated separately (now included in updated tow info)
# target_depth <- read.csv(
#   here::here("data", "bcsi_bridge_target_depth_20221009.csv")) %>% 
#   janitor::clean_names()




## INTERSECTION WITH IPES GRID -------------------------------------------------

## import and transform IPES grid
# ipes_sf_deg <- readRDS(here::here("data", "spatial", "ipes_sf_list_deg.RDS"))
ipes_sf_poly <- st_read(
  here::here("data", "spatial", "ipes_shapefiles", "IPES_Grid_UTM9.shp")) %>% 
  #convert into single polygon
  st_union(., by_feature = FALSE)


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
  mutate(
    utm_x = get_utm(start_longitude, start_latitude, zone = "9", loc = "x"),
    utm_y = get_utm(start_longitude, start_latitude, zone = "9", loc = "y")
  ) 

bridge_sf <- bridge %>% 
  select(unique_event, utm_y, utm_x) %>% 
  st_as_sf(., coords = c("utm_x", "utm_y"), 
           crs = st_crs(ipes_sf_poly))

# check
ggplot() +
  geom_sf(data = st_intersection(bridge_sf, ipes_sf_poly))

# extract events within ipes survey grid
ipes_grid_events <- st_intersection(bridge_sf, ipes_sf_poly) %>% 
  pull(., unique_event)


## MISSING BATHY DATA ----------------------------------------------------------

missing_bathy <- bridge_raw %>% filter(
  is.na(mean_bottom_depth)
)

bathy_grid <- marmap::getNOAA.bathy(
  lon1 = min(missing_bathy$start_longitude),
  lon2 = max(missing_bathy$start_longitude),
  lat1 = min(missing_bathy$start_latitude),
  lat2 = max(missing_bathy$start_latitude),
  resolution = 0.5
)

bathy_df <- data.frame(
  lon = rep(rownames(bathy_grid), 
            each = length(unique(colnames(bathy_grid)))),
  lat = rep(colnames(bathy_grid), 
            times = length(unique(rownames(bathy_grid)))),
  depth = NA
) 
for (i in 1:nrow(bathy_grid)) {
  bathy_df$depth[i] <- bathy_grid[bathy_df$lon[i], bathy_df$lat[i]]
}

# subset to remove missing depths that are on land or inshore of 50 m isobath
bathy_df_trim <- bathy_df %>% 
  filter(!is.na(depth),
         !depth > -50)

# identify closest depth data
tt <- hutilscpp::match_nrst_haversine(
  lat = missing_bathy$start_latitude,
  lon = missing_bathy$start_longitude,
  addresses_lat = as.numeric(bathy_df_trim$lat),
  addresses_lon = as.numeric(bathy_df_trim$lon)
)

missing_bathy$est_mean_bottom_depth <- -1 * bathy_df_trim$depth[tt$pos] 

ggplot(missing_bathy) +
  geom_point(aes(x = start_longitude, y = start_latitude, fill = est_mean_bottom_depth),
             shape = 21)
ggplot() +
  geom_point(data = bridge_raw, 
             aes(x = start_longitude, y = start_latitude, fill = mean_bottom_depth),
             shape = 21) +
  geom_point(data = missing_bathy,
             aes(x = start_longitude, y = start_latitude, fill = est_mean_bottom_depth),
             shape = 23, colour = "red")

## CLEAN AND MERGE -------------------------------------------------------------

# impute missing bridge data (generally because distance travelled not tracked);
# based on vessel name, year, location, target_depth, net_desc, mouth_km2
imp_dat <- bridge %>% 
  mutate(
    # replace 0s
    volume_km3 = ifelse(volume_km3 == "0", NaN, volume_km3)
  ) %>% 
  select(utm_x, utm_y, trip_id, target_depth, tow_length_minute, net_desc, 
         mouth_km2, 
         bath_depth_mean_m, dist_to_coast_km, distance_km:height_km) 
imp_dat2 <- VIM::kNN(imp_dat, k = 5)
imp_dat2 <- imp_dat2 %>% 
  select(-ends_with("imp")) %>% 
  mutate(unique_event = bridge$unique_event)
  
dat <- bridge %>% 
  # remove values and replace with imputed data above
  select(-c(utm_x, utm_y, trip_id, tow_length_hour, target_depth, vessel_speed, 
            bath_depth_mean_m, dist_to_coast_km, distance_km:height_km)) %>% 
  left_join(., imp_dat2, by = "unique_event") %>% 
  # left_join(., chin, by = "unique_event") %>%
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
      year > 2016 & season_f == "su" & ipes_grid == TRUE,
      "ipes", 
      "hss"
    ) %>% 
      as.factor(),
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
    ),
    volume_m3 = (distance_km * 1000) * (height_km * 1000) * (width_km * 1000)
  ) %>% 
  dplyr::select(unique_event, date, year, month, week, day, day_night, season_f,
         stratum:station_name, pfma = dfo_stat_area_code,
         synoptic_station, ipes_grid, survey_f,
         mean_lat, mean_lon, utm_x, utm_y,
         vessel, distance_km, vessel_speed, target_depth, target_depth_bin, 
         bath_depth_mean_m, dist_to_coast_km, height_km, width_km, 
         volume_m3
         # , ck_juv = n_juv, ck_ad = n_ad
         ) 


# subset to core area and remove rows with missing data
dat_trim <- dat %>%
  filter(synoptic_station == TRUE,
         # exclude SoG and Puget Sound PFMAs
         !pfma %in% c("13", "14", "15", "16", "17", "18", "19", "PS")) %>%
  droplevels()


saveRDS(dat_trim, here::here("data", "survey_sbc.rds"))
saveRDS(coast, here::here("data", "spatial", "sbc_sf_utm.rds"))


## SURVEY FIGS -----------------------------------------------------------------


dat_trim <- readRDS(here::here("data", "survey_sbc.rds"))
coast <- readRDS(here::here("data", "spatial", "sbc_sf_utm.rds"))

min_lat <- min(floor(dat_trim$mean_lat) - 0.1)
max_lat <- max(dat_trim$mean_lat) + 0.1
min_lon <- min(floor(dat_trim$mean_lon) - 0.1)
max_lon <- max(dat_trim$mean_lon) + 0.1

coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

# 
# ggplot() +
#   geom_sf(data = coast, color = "black", fill = "white") +
#   geom_point(data = dat_trim %>% filter(!bath_depth_mean_m > 0),
#              aes(x = utm_x, y = utm_y, fill = survey_f),
#              shape = 21) +
#   scale_fill_discrete(name = "Survey") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title = element_blank()) 
  


# map of set locations
set_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = dat_trim,
             aes(x = utm_x, y = utm_y, fill = survey_f), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Survey") +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
  
png(here::here("figs", "ms_figs", "set_map.png"), height = 4, width = 4, 
    units = "in", res = 250)
set_map
dev.off()

  
# set locations by month
set_map_month <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = dat_trim,
             aes(x = utm_x, y = utm_y, fill = survey_f), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Survey") +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  facet_wrap(~month)

seasonal_effort <- dat_trim %>% 
  group_by(week) %>% 
  tally() %>% 
  ggplot(.) +
  geom_point(aes(x = week, y = n)) +
  labs(y = "Number of Tows in Dataset")

pdf(here::here("figs", "set_coverage.pdf"), height = 4, width = 8)
set_map_month
seasonal_effort
dev.off()



# map of spatial covariates (exclude for now)
# bc_raster <- readRDS(
#   here::here("data", "spatial", "full_coast_raster_latlon_1000m.RDS"))
# 
# # # merge and add aspect/slope
# ipes_raster_slope <- terrain(ipes_raster_utm, opt = 'slope', unit = 'degrees',
#                              neighbors = 8)
# ipes_raster_aspect <- terrain(ipes_raster_utm, opt = 'aspect', unit = 'degrees',
#                               neighbors = 8)
# ipes_raster_list <- list(depth = ipes_raster_utm,
#                          slope = ipes_raster_slope,
#                          aspect = ipes_raster_aspect)


# stacked bar plots for headrope depth and ppn day/night
stacked_headrope_depth <- dat_trim %>%
  group_by(target_depth_bin, year) %>%
  summarize(n = length(unique_event), .groups = "drop") %>% 
  ggplot(., aes(x = as.factor(year), y = n, fill = target_depth_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  scale_fill_brewer(type = "seq", palette = 3, 
                    name = "Target\nHeadrope\nDepth",
                    direction = -1) +
  labs(y = "Number of Tows") +
  theme(
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(breaks = seq(1995, 2020, by = 5))


png(here::here("figs", "ms_figs", "headrope.png"), height = 3.5, width = 7, 
    units = "in", res = 250)
stacked_headrope_depth
dev.off()


stacked_daynight <- dat_trim %>%
  group_by(day_night, year) %>%
  summarize(n = length(unique_event), .groups = "drop") %>% 
  ggplot(., aes(x = as.factor(year), y = n, fill = day_night)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  scale_fill_brewer(type = "qual", palette = 2, 
                    name = "") +
  labs(y = "Number of Tows") +
  theme(
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(breaks = seq(1995, 2020, by = 5))

png(here::here("figs", "ms_figs", "dn_tows.png"), height = 3.5, width = 7, 
    units = "in", res = 250)
stacked_daynight
dev.off()


# bubble plots of temporal coverage
bubble_temp_coverage <- dat_trim %>% 
  group_by(year, week, survey_f) %>% 
  summarize(n_tows = length(unique_event), .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(.) +
  geom_jitter(aes(x = week, y = year, size = n_tows, 
                  fill = survey_f),
              alpha = 0.3, width = 0.25, shape = 21) +
  ggsidekick::theme_sleek() +
  scale_size(name = "Number\nof\nTows") +
  scale_fill_discrete(name = "Survey") +
  theme(
    axis.title.y = element_blank()
  ) +
  labs(x = "Week")

png(here::here("figs", "ms_figs", "temp_cov.png"), height = 5.5, width = 5.5, 
    units = "in", res = 250)
bubble_temp_coverage
dev.off()


# histogram of swept volume
hist_vol_swept <- ggplot(dat_trim) +
  geom_histogram(aes(x = volume_m3), bins = 50, alpha = 0.6) +
  ggsidekick::theme_sleek() +
  labs(x = "Volume Swept (m^3)", y = "Number of Sets")

png(here::here("figs", "ms_figs", "vol_swept.png"), height = 3.5, width = 5.5, 
    units = "in", res = 250)
hist_vol_swept
dev.off()


## ADD CATCH DATA --------------------------------------------------------------

BCSI_SalmonCounts_2023215.csv

# import all species
sp_files <- list.files(path = here::here("data"), pattern = "BCSI_Juv_")
sp_catch <- purrr::map(
  sp_files, function (x) {
    read.csv(paste(here::here("data"), x, sep = "/"),
             fileEncoding="UTF-8-BOM")
  }
) %>% 
  bind_rows() %>% 
  janitor::clean_names()

sp_catch2 <- read.csv(here::here("data", "BCSI_SalmonCounts_2023215.csv"),
                      fileEncoding="UTF-8-BOM")

catch_dat <- left_join(dat_trim, sp_catch, by = "unique_event") 

saveRDS(catch_dat, here::here("data", "catch_survey_sbc.rds"))


## DATA EXPLORE ----------------------------------------------------------------

catch_dat_plot <- catch_dat %>% 
  mutate(
    log_n_juv = n_juv + 0.0001,
    non_zero = ifelse(n_juv > 0, "yes", "no")
  )
shape_pal <- c(21, 3)
names(shape_pal) <- c("yes", "no")

sp_dist <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = catch_dat_plot,
             aes(x = utm_x, y = utm_y, 
                 # fill = survey_f, 
                 size = log_n_juv,
                 shape = non_zero), 
             # shape = 21,
             fill = "red",
             alpha = 0.4) +
  # scale_fill_discrete(name = "Survey") +
  scale_shape_manual(values = shape_pal, name = "Salmon Caught") +
  scale_size_continuous(name = "Individuals Caught", trans = "log") +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  facet_grid(species~season_f)
sp_dist2 <- sp_dist +
  scale_size_continuous(name = "Individuals Caught")
sp_dist3 <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_point(data = catch_dat_plot %>% filter(non_zero == "yes"),
             aes(x = utm_x, y = utm_y, 
                 size = log_n_juv),
             shape = 21,
             fill = "red",
             alpha = 0.4) +
  scale_size_continuous(name = "Individuals Caught") +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  facet_grid(species~season_f)


pdf(here::here("figs", "ms_figs", "raw_catch.pdf"), height = 12, width = 8)
sp_dist
sp_dist2
sp_dist3
dev.off()  


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


png(here::here("figs", "data_input_map_year.png"), 
    height = 5, width = 16, res = 250, 
    units = "in")
base_map +
  geom_point(data = dat %>% filter(!year == "2022"),
             aes(x = utm_x, y = utm_y, fill = survey_f), 
             shape = 21, alpha = 0.4) +
  facet_grid(season_f ~ year)
dev.off()