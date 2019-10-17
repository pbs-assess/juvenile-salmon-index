## Clean and merge data from IPES and High Seas Access databases
# Note: 32-bit R needs to be used to import data from Access

library(RODBC); library(tidyverse); library(ggplot2)

# Access databases saved to regional drive to import bridge (i.e. total
# catches) and genetics data from historical (high seas) and contemporary (IPES)
# databases
pathHighSeas <- "R:/SCIENCE/IPES/DB/chinookIndex/HSSALMON.accdb" 
conHS <- odbcConnectAccess2007(pathHighSeas)
pathIPES <- "R:/SCIENCE/IPES/DB/IPES_TrawlDB_v19.07f_2017_18_19.mdb"
conIPES <- odbcConnectAccess2007(pathIPES)

# High seas bridge data - includes catch totals
bridgeQry <- "SELECT BRIDGE.STATION_ID, BRIDGE.EVENT, STATION_INFO.REGION, 
  BRIDGE.STATION, BRIDGE.YEAR, BRIDGE.MONTH, BRIDGE.DAY, BRIDGE.DATE, 
  BRIDGE.JULIAN_DATE, BRIDGE.START_TIME, BRIDGE.START_LAT, BRIDGE.START_LONG, 
  BRIDGE.END_LAT, BRIDGE.END_LONG, BRIDGE.DISTANCE, BRIDGE.DUR, 
  BRIDGE.[SOG-KTS], BRIDGE.HEADING, BRIDGE.START_BOT_DEPTH, 
  BRIDGE.END_BOT_DEPTH, BRIDGE.NET_OPENING_WIDTH, BRIDGE.NET_OPENING_HEIGHT, 
  BRIDGE.HEAD_DEPTH, BRIDGE.CK_JUV, BRIDGE.CK_ADULT, BRIDGE.GEAR_TYPE, 
  BRIDGE.COMMENTS
FROM STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = 
  BRIDGE.STATION_ID;"
bridgeHS <- sqlQuery(conHS, bridgeQry) 

# IPES bridge data (no catch); separate from Chin catches to keep sets with 0s
bridgeIPESQuery <- "SELECT BRIDGE_LOG_ID,
  BRIDGE_LOG_FIELD_ID, TRIP_ID, 
  EVENT_DATE, EVENT_NUMBER, TOW_NUMBER, 
  EVENT_TYPE, EVENT_SUB_TYPE, GEAR_CODE, 
  BEGIN_DEPLOYMENT_TIME, END_DEPLOYMENT_TIME, 
  BEGIN_RETRIEVAL_TIME, END_RETRIEVAL_TIME, 
  TOW_DURATION, 
  START_LATITUDE, START_LONGITUDE, 
  END_LATITUDE, END_LONGITUDE, 
  START_BOTTOM_DEPTH, END_BOTTOM_DEPTH, 
  AVG_GEAR_DEPTH, AVG_TOW_DEPTH, 
  DISTANCE_TRAVELLED, USABILITY_CODE, COMMENTS
FROM BRIDGE_LOG
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow'));
"
sqlQuery(conIPES, bridgeIPESQuery) %>% 
  glimpse()
bridgeIPES <- sqlQuery(conIPES, bridgeIPESQuery)

# IPES Chinook catches (only non-zero tows)
chinIPESQuery <- "SELECT CATCH_ID, CATCH_FIELD_ID, 
  BRIDGE_LOG_ID, SPECIES_CODE, CATCH_COUNT, 
  JUVENILE_CATCH_COUNT, CATCH_WEIGHT, JUVENILE_CATCH_WEIGHT, 
  COMMENTS
FROM CATCH
WHERE (((CATCH.SPECIES_CODE)='124'));
"
sqlQuery(conIPES, chinIPESQuery) %>% 
  glimpse()
trimChinIPES <- sqlQuery(conIPES, chinIPESQuery) %>% 
  select(BRIDGE_LOG_ID, ADULT_CATCH = CATCH_COUNT, 
         JUV_CATCH = JUVENILE_CATCH_COUNT)


## Clean IPES data 
# add catches to bridge; trim to match high seas
trimBridgeIPES <- bridgeIPES %>% 
  left_join(., trimChinIPES, by = "BRIDGE_LOG_ID") %>% 
  replace_na(list(ADULT_CATCH = 0, JUV_CATCH = 0)) %>% 
  mutate(STATION_ID = paste("IPES", TRIP_ID, BRIDGE_LOG_ID, sep = "-"),
         AVG_BOTTOM_DEPTH = (START_BOTTOM_DEPTH + END_BOTTOM_DEPTH) / 2,
         DATASET = "IPES") %>% 
  select(STATION_ID, DATASET, DATE = EVENT_DATE, 
         START_TIME = BEGIN_DEPLOYMENT_TIME,
         START_LAT = START_LATITUDE, START_LONG = START_LONGITUDE,
         # END_LAT = END_LATITUDE, END_LONG = END_LONGITUDE, 
         DUR = TOW_DURATION,
         AVG_BOTTOM_DEPTH, HEAD_DEPTH = AVG_GEAR_DEPTH, CK_JUV = JUV_CATCH,
         CK_ADULT = ADULT_CATCH)


## Plot overlap between IPES and high seas
library(ggmap)
nAm <- map_data("world") %>% 
  filter(region %in% c("Canada", "USA"))

ggplot(bridgeHS) +
  geom_point(aes(x = START_LONG, y = START_LAT, color = REGION)) +
  geom_point(data = trimBridgeIPES,
             aes(x = START_LONG, y = START_LAT)) +
  geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
           color = "black", fill = "gray80") + 
  coord_fixed(xlim = c(-129.5, -123), ylim = c(48, 52), ratio = 1.3)


## Clean High Seas data 
# retain only regions that overlap
keepRegions <- c("INSIDE VANCOUVER ISLAND", "JOHNSTONE STRAIT", "QUEEN 
                 CHARLOTTE SOUND", "QUEEN CHARLOTTE STRAIT", "VANCOUVER ISLAND")
# exclude tows based on information in comments of access database
excludeTows <- read.csv(here::here("data", "excludeTows.csv")) %>% 
  filter(EXCLUDE == "Y") 

trimBridgeHS <- bridgeHS %>% 
  filter(REGION %in% keepRegions,
         !STATION_ID %in% excludeTows$STATION_ID,
         !is.na(START_LAT)) %>% 
  mutate(AVG_BOTTOM_DEPTH = (START_BOT_DEPTH + END_BOT_DEPTH) / 2, 
         DATASET = "HighSeas") %>%
  replace_na(list(CK_ADULT = 0, CK_JUV = 0)) %>% 
  select(STATION_ID, DATASET, DATE, START_TIME, START_LAT, START_LONG, 
         # END_LAT, END_LONG,
         DUR, AVG_BOTTOM_DEPTH, HEAD_DEPTH, CK_JUV, CK_ADULT)

## Merge both bridge datasets
dum <- rbind(trimBridgeHS, trimBridgeIPES) %>% 
  mutate(MONTH = lubridate::month(DATE),
         YEAR = lubridate::year(DATE),
         jDay = lubridate::yday(DATE),
         START_TIME = strftime(START_TIME, format = "%H:%M:%S")) %>% 
  rename_all(., tolower)

## Convert coordinates to UTM
coords <- dum %>% 
  select(station_id, yUTM_start = start_lat, xUTM_start = start_long)
sp::coordinates(coords) <- c("xUTM_start", "yUTM_start")
sp::proj4string(coords) <- sp::CRS("+proj=longlat +datum=WGS84")
#SE Van Island extends into zone 10; not sure if that's an issue or not
coords2 <- sp::spTransform(coords, sp::CRS("+proj=utm +zone=9 ellps=WGS84")) %>%
  as(., "SpatialPoints")

bridgeOut <- dum %>% 
  cbind(., coords2@coords) %>% 
  select(station_id, dataset, date, jday, month, year, start_time, start_lat, 
       start_long, xUTM_start, yUTM_start, dur, avg_bottom_depth, head_depth,
       ck_juv, ck_adult)
saveRDS(bridgeOut, here::here("data", "mergedCatch.rds"))

## Clean genetics data ##

## NEED TO MATCH TO STATION NUMBER ##

# High seas Chinook GSI data
chinDNAQry <- "SELECT FISH_NUMBER, STATION_INFO.REGION, BRIDGE.Year, 
  BRIDGE.Month, BRIDGE.Day, BRIDGE.START_LAT, BRIDGE.START_LONG, 
  SPECIES, SHIP_FL, SHIP_TL, SHIP_WT, 
  [BATCH-DNA_NUMBER], STOCK_1, 
  STOCK_2, STOCK_3, 
  STOCK_4, STOCK_5, 
  REGION_1, REGION_2,
  REGION_3, REGION_4, 
  REGION_5, PROB_1, 
  PROB_2, PROB_3, 
  PROB_4, PROB_5
  FROM STATION_INFO 
INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON BRIDGE.STATION_ID = 
  BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN (BIOLOGICAL INNER JOIN 
  DNA_CHINOOK_STOCK_ID ON BIOLOGICAL.FISH_NUMBER = 
  DNA_CHINOOK_STOCK_ID.FISH_NUMBER) ON (BIOLOGICAL_JUNCTION.FISH_NUMBER = 
  DNA_CHINOOK_STOCK_ID.FISH_NUMBER) AND (BIOLOGICAL_JUNCTION.FISH_NUMBER = 
  BIOLOGICAL.FISH_NUMBER)) ON (STATION_INFO.STATION_ID = BRIDGE.STATION_ID) AND 
  (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((BIOLOGICAL.SPECIES)='chinook'));
"
chinHS <- sqlQuery(conHS, chinDNAQry)


# IPES Chinook GSI
dnaIPESQuery <- "SELECT DNA_STOCK_INDIVIDUAL_FISH_ID, BCSI_FISH_NUMBER, 
  BATCH_DNA_NUMBER, STOCK_1, STOCK_2, STOCK_3, STOCK_4, STOCK_5, REGION_1, 
  REGION_2, REGION_3, REGION_4, REGION_5, PROB_1, PROB_2, PROB_3, PROB_4, 
  PROB_5
FROM DNA_STOCK_INDIVIDUAL_FISH;
"
dnaIPES <- sqlQuery(conIPES, dnaIPESQuery)

#Samples from specimen table have differently formatted fish identifier
fishID <- dnaIPES$BCSI_FISH_NUMBER %>% 
  as.vector() %>% 
  strsplit(., split = "-") %>% 
  unlist() %>%
  matrix(., nrow = 5, ncol = length(dnaIPES$BCSI_FISH_NUMBER)) %>%
  t() %>%
  data.frame() %>%
  dplyr::rename("prog" = X1, "survey" = X2, "event" = X3, "species" = X4,
                "conFish" = X5) %>% 
  mutate(BCSI_FISH_NUMBER = dnaIPES$BCSI_FISH_NUMBER) %>% 
  filter(species %in% c("124", "124J", "124A")) 

dnaIPESTrim <- dnaIPES %>% 
  right_join(., fishID, by = "BCSI_FISH_NUMBER") 
# %>% 
#   mutate(
#     year = as.numeric(
#       case_when(
#         survey == "201873" ~ "2018",
#         survey == "201742" ~ "2017",
#         TRUE ~ "NA")),
#     survey = case_when(
#       survey == "201873" ~ "73",
#       survey == "201742" ~ "42",
#       TRUE ~ "NA"),
#     age = case_when(
#       grepl("J", species) ~ "J",
#       grepl("A", species) ~ "A"),
#     # species = "124",
#     UNIVERSAL_FISH_LABEL = prog %>% #matches identifier from next query
#       paste(., year, sep = "") %>% 
#       paste(., survey, event, species, sep = "-") %>% 
#       paste(., paste("0", conFish, sep = ""), sep = "-")
#     )
  
# IPES sampling key (necessary to match to station_id)
chinIPESQuery <- "SELECT SPECIMEN.SPECIMEN_ID, BRIDGE_LOG.BRIDGE_LOG_ID, BRIDGE_LOG.EVENT_NUMBER, TRIP.TRIP_YEAR, TRIP.TRIP_START_DATE, SPECIMEN.CATCH_ID, SPECIMEN.FISH_ID, SPECIMEN.UNIVERSAL_FISH_LABEL, SPECIMEN.SPECIES_CODE, SPECIMEN.LENGTH, SPECIMEN.WEIGHT
FROM TRIP INNER JOIN ((BRIDGE_LOG INNER JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) INNER JOIN SPECIMEN ON CATCH.CATCH_ID = SPECIMEN.CATCH_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((SPECIMEN.SPECIES_CODE)='124'));"
chinIPES <- sqlQuery(conIPES, chinIPESQuery) 
  
# Join stock identification data, new survey ID and individual sampling data 
dum <- chinIPES$UNIVERSAL_FISH_LABEL %>% 
  as.vector() %>% 
  strsplit(., split = "-") %>% 
  unlist() %>%
  matrix(., nrow = 5, ncol = length(chinIPES$UNIVERSAL_FISH_LABEL)) %>%
  t() %>%
  data.frame() %>%
  dplyr::rename("prog" = X1, "survey" = X2, "event" = X3, "species" = X4,
                "conFish" = X5) %>% 
  mutate(UNIVERSAL_FISH_LABEL = chinIPES$UNIVERSAL_FISH_LABEL) %>% 
  filter(species %in% c("124", "124J", "124A")) %>% 
  mutate(
    year = as.numeric(
      case_when(
        grepl("2017", prog) ~ "2017",
        grepl("2018", prog) ~ "2018",
        grepl("2019", prog) ~ "2019")),
    prog = case_when(
      grepl("IPES", prog) ~ "IPES",
      TRUE ~ "NA"),
    conFish = formatC(conFish, width = 3, format = "d", flag = "0")
  )

dum2 <- chinIPES %>% 
  left_join(., dum, by = "UNIVERSAL_FISH_LABEL") %>% 
  mutate(
    age = case_when(
      grepl("J", species) ~ "J",
      grepl("A", species) ~ "A",
      LENGTH <= 300 ~ "J",
      LENGTH > 300 ~ "A",
      TRUE ~ "NA"
    ),
    species = paste("124", age, sep = ""),
    BCSI_FISH_NUMBER = paste(year, survey, sep = "") %>% 
      paste(prog, ., event, species, conFish, sep = "-")
  )

dnaIPESOut <- dum2 %>% 
  right_join(., dnaIPESTrim, by = "BCSI_FISH_NUMBER") 

