## Clean BCSI Data
# Uses .csv files produced by Amy Tabata May 2022 (updated Feb 2023)

library(tidyverse)
library(sf)
library(sp)


bridge_raw <- read.csv(here::here("data", "BCSI_TowInfo_2023215.csv")) %>% 
  #read.csv(here::here("data", "BCSI_Bridge_Info_20220929.csv")) %>% 
  janitor::clean_names()  %>%
  # exclude some absurdly deep tows
  filter(!target_depth > 75,
         !usable == "N") 


# distance to coast and mean depth (still some gaps)
depth_bathy_dat <- read.csv(
  here::here("data", "full_depths_from_bathy_2023215.csv")) %>% 
  janitor::clean_names() 


# event dates (missing from updated bridge log)
event_dates <- read.csv(here::here("data", "BCSI_event_dates.csv")) %>% 
  janitor::clean_names()


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

bridge <- left_join(bridge_raw, event_dates, by = "unique_event") %>% 
  mutate(
    utm_x = get_utm(start_longitude, start_latitude, zone = "9", loc = "x"),
    utm_y = get_utm(start_longitude, start_latitude, zone = "9", loc = "y"),
    # missing day/night classifier for 2022 fall survey 
    day_night = ifelse(trip_id == "2022-011", "DAY", day_night),
    date = as.POSIXct(date, format="%m/%d/%Y"),
    month = lubridate::month(date),
    week = lubridate::week(date)
  ) 

# still have some missing dates
missing_dates <- bridge %>% 
  filter(is.na(date))

bridge_sf <- bridge %>% 
  select(unique_event, utm_y, utm_x) %>% 
  st_as_sf(., coords = c("utm_x", "utm_y"), 
           crs = st_crs(ipes_sf_poly))

# check
# ggplot() +
#   geom_sf(data = st_intersection(bridge_sf, ipes_sf_poly))

# extract events within ipes survey grid
ipes_grid_events <- st_intersection(bridge_sf, ipes_sf_poly) %>% 
  pull(., unique_event)


## MISSING BATHY DATA ----------------------------------------------------------

missing_bathy <- bridge %>%
  left_join(., depth_bathy_dat, by = "unique_event") %>% 
  filter(
    is.na(depth_mean_m) | depth_mean_m < 0)

bathy_grid <- marmap::getNOAA.bathy(
  lon1 = min(missing_bathy$start_longitude),
  lon2 = max(missing_bathy$start_longitude),
  lat1 = min(missing_bathy$start_latitude),
  lat2 = max(missing_bathy$start_latitude),
  resolution = 1
)

bathy_matrix <- matrix(data = NA, 
                       nrow = nrow(bathy_grid), ncol = ncol(bathy_grid))
for (i in 1:nrow(bathy_matrix)) {
  bathy_matrix[i, ] <- bathy_grid[i, ] %>% as.numeric()
}
dimnames(bathy_matrix) <- dimnames(bathy_grid)
bathy_df <- data.frame(
  lon = rownames(bathy_matrix),
  bathy_matrix) %>% 
  pivot_longer(., cols = -(lon), names_to = "lat",
               names_prefix = "X", values_to = "depth") %>%
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat)) %>% 
  glimpse()

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
missing_bathy$pred_depth_mean_m <- -1 * bathy_df_trim$depth[tt$pos]


## CLEAN AND MERGE -------------------------------------------------------------

# impute missing bridge data (generally because distance travelled not tracked);
# based on vessel name, year, location, target_depth, net_desc, mouth_km2
imp_eff <- left_join(bridge, 
                     depth_bathy_dat %>% select(unique_event, dist_to_coast_km), 
                     by = "unique_event") %>% 
  mutate(
    # replace 0s
    volume_km3 = ifelse(volume_km3 == "0", NaN, volume_km3),
    mouth_km2 = ifelse(mouth_km2 == "0", NaN, mouth_km2),
    tow_length_minute = ifelse(tow_length_minute == "0", NaN, tow_length_minute),
    distance_km = ifelse(distance_km == "0", NaN, distance_km),
    dist_to_coast_km = ifelse(dist_to_coast_km == "0", NaN, dist_to_coast_km),
  ) %>% 
  select(unique_event, year, week,
         utm_x, utm_y, target_depth, vessel_name, net_desc, mouth_km2, 
         tow_length_minute, distance_km, volume_km3, dist_to_coast_km) %>% 
  VIM::kNN(.) 
imp_eff$unique_event <- bridge$unique_event

sum(imp_eff$mouth_km2_imp) / nrow(imp_eff)
sum(imp_eff$tow_length_minute_imp) / nrow(imp_eff)
sum(imp_eff$distance_km_imp) / nrow(imp_eff)
sum(imp_eff$volume_km3_imp) / nrow(imp_eff)
sum(imp_eff$dist_to_coast_km_imp) / nrow(imp_eff)
# impute <1% to 3%


bridge2 <- bridge %>% 
  #remove data imputed above
  select(-c(utm_x, utm_y, target_depth, vessel_name, net_desc, mouth_km2, 
            tow_length_minute, distance_km, volume_km3, year, week)) %>% 
  # add bathy data
  left_join(
    ., 
    depth_bathy_dat %>% select(unique_event, depth_mean_m),
    by = "unique_event"
  ) %>% 
  # add missing bathy data
  left_join(
    ., 
    missing_bathy %>% select(unique_event, pred_depth_mean_m), 
    by = "unique_event"
  ) %>% 
  # add imputed effort data
  left_join(
    ., 
    imp_eff %>% select(-ends_with("imp")),
    by = "unique_event") %>%
  mutate(
    depth_mean_m = ifelse(
      is.na(pred_depth_mean_m), depth_mean_m, pred_depth_mean_m
    ),
    season_f = case_when(
      month %in% c("2", "3", "4") ~ "sp",
      month %in% c("5", "6", "7", "8") ~ "su",
      month %in% c("9", "10", "11", "12") ~ "wi"
    ) %>% 
      as.factor(.),
    # define core area as southern BC and CBC excluding SoG
    synoptic_station = ifelse(
      start_latitude > 47 & start_latitude < 53 & !grepl("GS", unique_event) & 
        start_longitude > -135,
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
    #used for plots even though not fit in model
    target_depth_bin = cut(
      target_depth, 
      breaks = c(-Inf, 14, 29, 44, 59, Inf), 
      labels=c("0", "15", "30", "45", "60")
    ) %>% 
      as.factor(),
    volume_m3 = volume_km3 * 1e09
  ) %>% 
  dplyr::select(
    unique_event, date, year, year_f, month, week, day, day_night, season_f,
    pfma = dfo_stat_area_code, synoptic_station, ipes_grid, survey_f,
    lat = start_latitude, lon = start_longitude, utm_x, utm_y,
    target_depth, target_depth_bin, dist_to_coast_km, volume_m3, volume_km3, 
    depth_mean_m
  )
  
  
# subset to core area and remove rows with missing data
dat_trim <- bridge2 %>%
  filter(synoptic_station == TRUE,
         # exclude low sample size years
         year >= 1998,
         # exclude SoG and Puget Sound PFMAs
         !pfma %in% c("13", "14", "15", "16", "17", "18", "19", "PS"),
         # remove one station clearly on land
         !unique_event  == "HS201466-JF02") %>%
  droplevels()


saveRDS(dat_trim, here::here("data", "survey_sbc.rds"))


## SURVEY FIGS -----------------------------------------------------------------


dat_trim <- readRDS(here::here("data", "survey_sbc.rds")) 
  
min_lat <- min(floor(dat_trim$lat) - 0.25)
max_lat <- max(dat_trim$lat) + 0.25
min_lon <- min(floor(dat_trim$lon) - 0.25)
max_lon <- max(dat_trim$lon) + 0.25

coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  st_crop(., xmin = min_lon, ymin = min_lat, xmax = max_lon, ymax = max_lat) %>%
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))
  # st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m"))

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
  geom_sf(data = coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = dat_trim,
             aes(x = lon, y = lat, fill = survey_f), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Survey") +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15), 
           xlim = c(min_lon + 0.15, max_lon - 0.15))
  
png(here::here("figs", "ms_figs_season", "set_map.png"), height = 5, width = 5, 
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

pdf(here::here("figs", "diagnostics",  "set_coverage.pdf"), height = 4, width = 8)
set_map_month
seasonal_effort
dev.off()



# set locations by year
set_map_year <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white", size = 1.25) +
  geom_point(data = dat_trim,
             aes(x = lon, y = lat, fill = season_f), 
             shape = 21, alpha = 0.4) +
  scale_fill_discrete(name = "Survey") +
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~year) +
  # hacky way to ensure borders are correct
  coord_sf(ylim = c(min_lat + 0.15, max_lat - 0.15), 
           xlim = c(min_lon + 0.15, max_lon - 0.15))


pdf(here::here("figs", "diagnostics",  "set_coverage_year.pdf"))
set_map_year
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
  group_by(target_depth_bin, year_f) %>%
  summarize(n = length(unique_event), .groups = "drop") %>% 
  ggplot(., aes(x = year_f, y = n, fill = target_depth_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  scale_fill_brewer(type = "seq", palette = 3, 
                    name = "Target\nHeadrope\nDepth",
                    direction = -1) +
  labs(y = "Number of Tows") +
  theme(
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(breaks = seq(1998, 2020, by = 5))


png(here::here("figs", "ms_figs_season", "headrope.png"), height = 3.5, width = 7, 
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

png(here::here("figs", "ms_figs_season", "dn_tows.png"), height = 3.5, width = 7, 
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

png(here::here("figs", "ms_figs_season", "temp_cov.png"), height = 5.5, width = 5.5, 
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

sp_dat <- read.csv(here::here("data", "BCSI_SalmonCounts_2023215.csv")) %>% 
  janitor::clean_names() %>% 
  mutate(
    species = case_when(
      scientific_name == "ONCORHYNCHUS GORBUSCHA" ~ "pink",
      scientific_name == "ONCORHYNCHUS KETA" ~ "chum",
      scientific_name == "ONCORHYNCHUS KISUTCH" ~ "coho",
      scientific_name == "ONCORHYNCHUS MYKISS" ~ "steelhead",
      scientific_name == "ONCORHYNCHUS NERKA" ~ "sockeye",
      scientific_name == "ONCORHYNCHUS TSHAWYTSCHA" ~ "chinook"
    )
  ) %>% 
  select(
    -scientific_name
  )

# # import all species
# sp_files <- list.files(path = here::here("data"), pattern = "BCSI_Juv_")
# sp_catch <- purrr::map(
#   sp_files, function (x) {
#     read.csv(paste(here::here("data"), x, sep = "/"),
#              fileEncoding="UTF-8-BOM")
#   }
# ) %>% 
#   bind_rows() %>% 
#   janitor::clean_names()


catch_dat <- left_join(dat_trim, sp_dat, by = "unique_event") 

# catch data excludes stations w/ 0 catches for all species; correct and check
missing_catches <- catch_dat %>% 
  filter(is.na(n_total))
sp_vec <- unique(sp_dat$species)
catch_list <- vector(mode = "list", length = length(sp_vec))
for (i in seq_along(catch_list)) {
  catch_list[[i]] <- missing_catches %>% 
    mutate(
      n_juv = 0,
      n_ad = 0,
      n_total = 0,
      species = sp_vec[i]
    )
} 

catch_out <- catch_dat %>% 
  filter(!is.na(n_total)) %>% 
  rbind(., bind_rows(catch_list)) 


saveRDS(catch_out, here::here("data", "catch_survey_sbc.rds"))


## DATA EXPLORE ----------------------------------------------------------------

catch_dat <- readRDS(here::here("data", "catch_survey_sbc.rds"))


# spatial distribution of catches
ggplot() +
  geom_sf(data = coast) +
  geom_point(data = catch_dat %>% filter(n_juv > 0, !species == "steelhead"),
             aes(x = lon, y = lat, size = n_juv)) +
  facet_grid(season_f ~ species)


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