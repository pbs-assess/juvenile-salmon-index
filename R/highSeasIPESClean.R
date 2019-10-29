## Clean and merge data from IPES and High Seas Access databases
# Note: 32-bit R needs to be used to import data from Access

library(RODBC); library(tidyverse); library(ggplot2)

# Access databases saved to regional drive to import bridge (i.e. total
# catches) and genetics data from historical (high seas) and contemporary (IPES)
# databases
pathHighSeas <- "R:/SCIENCE/CFreshwater/chinookIndex/HSSALMON.accdb" 
conHS <- odbcConnectAccess2007(pathHighSeas)
# similar DB is saved in parent directory, but doesn't include GSI tables
pathIPES <- "R:/SCIENCE/CFreshwater/chinookIndex/IPES_TrawlDB_v19.07f_2017_18_19_wGSI.mdb"
conIPES <- odbcConnectAccess2007(pathIPES)


##### MERGE CATCH DATA #####

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

# IPES bridge data (no catch); separate from Chin catches to keep sets with 0s
bridgeIPESQuery <- "SELECT BRIDGE_LOG_ID,
  BRIDGE_LOG_FIELD_ID, BRIDGE_LOG.TRIP_ID, 
  EVENT_DATE, EVENT_NUMBER, TOW_NUMBER, 
  EVENT_TYPE, EVENT_SUB_TYPE, GEAR_CODE, 
  BEGIN_DEPLOYMENT_TIME, END_DEPLOYMENT_TIME, 
  BEGIN_RETRIEVAL_TIME, END_RETRIEVAL_TIME, 
  TOW_DURATION, 
  START_LATITUDE, START_LONGITUDE, 
  END_LATITUDE, END_LONGITUDE, 
  START_BOTTOM_DEPTH, END_BOTTOM_DEPTH, 
  AVG_GEAR_DEPTH, AVG_TOW_DEPTH, 
  DISTANCE_TRAVELLED, USABILITY_CODE
FROM BRIDGE_LOG
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow'));
"
sqlQuery(conIPES, bridgeIPESQuery) %>% 
  glimpse()
bridgeIPES <- sqlQuery(conIPES, bridgeIPESQuery) %>% 
  #add catches to bridge
  left_join(., trimChinIPES, by = "BRIDGE_LOG_ID") %>% 
  #remove unuseable tows
  filter(!USABILITY_CODE == "5") %>% 
  replace_na(list(ADULT_CATCH = 0, JUV_CATCH = 0)) %>% 
  mutate(STATION_ID = paste("IPES", TRIP_ID, BRIDGE_LOG_ID, sep = "-"),
         AVG_BOTTOM_DEPTH = (START_BOTTOM_DEPTH + END_BOTTOM_DEPTH) / 2,
         #convert to fractions of an hour to match HS database
         DUR = TOW_DURATION / 60, 
         DATASET = "IPES")

## Clean IPES data 
# Trime to match high seas
trimBridgeIPES <- bridgeIPES %>% 
  select(STATION_ID, BRIDGE_LOG_ID, DATASET, DATE = EVENT_DATE, 
         START_TIME = BEGIN_DEPLOYMENT_TIME,
         START_LAT = START_LATITUDE, START_LONG = START_LONGITUDE, DUR,
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

# ggplot(juvOut %>% filter(stableStation == "Y")) +
#   geom_point(aes(x = start_long, y = start_lat, color = dataset)) +
#   geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
#            color = "black", fill = "gray80") + 
#   coord_fixed(xlim = c(-129.5, -123), ylim = c(48, 52), ratio = 1.3)

## Clean High Seas data 
# exclude tows based on information in comments of access database (see bridge 
# table comments for details)
excludeTowsComments <- read.csv(here::here("data", "excludeTows.csv")) %>% 
  filter(EXCLUDE == "Y") 

trimBridgeHS <- bridgeHS %>% 
  filter(!STATION_ID %in% excludeTowsComments$STATION_ID,
         !is.na(START_LAT)) %>% 
  mutate(BRIDGE_LOG_ID = NA,
         AVG_BOTTOM_DEPTH = (START_BOT_DEPTH + END_BOT_DEPTH) / 2, 
         DATASET = "HighSeas") %>%
  replace_na(list(CK_ADULT = 0, CK_JUV = 0)) %>% 
  select(STATION_ID, BRIDGE_LOG_ID, DATASET, DATE, START_TIME, START_LAT, 
         START_LONG, DUR, AVG_BOTTOM_DEPTH, HEAD_DEPTH, CK_JUV, CK_ADULT)


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

excludeTowsRegion <- bridgeHS %>% 
  filter(!REGION %in% c("INSIDE VANCOUVER ISLAND", "JOHNSTONE STRAIT", 
                        "QUEEN CHARLOTTE SOUND", "QUEEN CHARLOTTE STRAIT", 
                        "VANCOUVER ISLAND")) %>% 
  select(STATION_ID, REGION)

bridgeOut <- dum %>% 
  cbind(., coords2@coords) %>% 
  mutate(
    #is station present in both IPES and HS datasets
    stableStation = case_when(
      station_id %in% excludeTowsRegion$STATION_ID ~ "N",
      TRUE ~ "Y"
    )) %>% 
  select(station_id, bridge_log_id, stableStation, dataset, date, jday, month, year, start_time, start_lat, 
       start_long, xUTM_start, yUTM_start, dur, avg_bottom_depth, head_depth,
       ck_juv, ck_adult)

### NOTE - CATCH ESTIMATES FROM BRIDGEOUT ARE SUSPECT 
## I.e. fish sampled from sets with catch recorded as 0. For now assume that 
# sampled fish is minimum catch estimate and do not use bridgeout independently.

##### MERGE GENETICS DATA #####

# High seas Chinook GSI data
chinDNAQry <- "SELECT BIOLOGICAL.FISH_NUMBER, STATION_INFO.STATION_ID,
  STATION_INFO.REGION, BRIDGE.Year, 
  BRIDGE.Month, BRIDGE.Day, BRIDGE.DATE, BRIDGE.START_LAT, BRIDGE.START_LONG, 
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
dnaHS <- sqlQuery(conHS, chinDNAQry) 
dnaHSOut <- dnaHS %>% 
  mutate(
    age = case_when(
      SHIP_FL <= 300 ~ "J",
      SHIP_FL > 300 ~ "A",
      TRUE ~ "NA"),
    jDay = lubridate::yday(DATE),
    month2 = lubridate::month(DATE),
    dataset = "HighSeas"
  ) %>% 
  rename_all(., tolower) %>% 
  select(fish_number, station_id, dataset, year, month = month2, jday, start_lat, 
         start_long, age, ship_fl, batch_dna_number = `batch-dna_number`, 
         stock_1:prob_5) 

# IPES Chinook GSI
dnaIPESQuery <- "SELECT DNA_STOCK_INDIVIDUAL_FISH_ID, BCSI_FISH_NUMBER, 
  BATCH_DNA_NUMBER, STOCK_1, STOCK_2, STOCK_3, STOCK_4, STOCK_5, REGION_1, 
  REGION_2, REGION_3, REGION_4, REGION_5, PROB_1, PROB_2, PROB_3, PROB_4, 
  PROB_5
FROM DNA_STOCK_INDIVIDUAL_FISH;
"
dnaIPES <- sqlQuery(conIPES, dnaIPESQuery)

#Samples from IPES specimen table have differently formatted fish identifier and 
#include non-Chinook species
dnaIPES_ck <- dnaIPES$BCSI_FISH_NUMBER %>% 
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
  filter(BCSI_FISH_NUMBER %in% dnaIPES_ck$BCSI_FISH_NUMBER) 
  
# IPES sampling key (necessary to match to station_id)
chinIPESQuery <- "SELECT SPECIMEN.SPECIMEN_ID, BRIDGE_LOG.BRIDGE_LOG_ID, 
  BRIDGE_LOG.EVENT_NUMBER, TRIP.TRIP_YEAR, TRIP.TRIP_START_DATE, 
  SPECIMEN.CATCH_ID, SPECIMEN.FISH_ID, SPECIMEN.UNIVERSAL_FISH_LABEL, 
  SPECIMEN.SPECIES_CODE, SPECIMEN.LENGTH, SPECIMEN.WEIGHT
FROM TRIP INNER JOIN ((BRIDGE_LOG INNER JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID 
  = CATCH.BRIDGE_LOG_ID) INNER JOIN SPECIMEN ON CATCH.CATCH_ID = 
  SPECIMEN.CATCH_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((SPECIMEN.SPECIES_CODE)='124'));"
chinIPES <- sqlQuery(conIPES, chinIPESQuery) 
  
# Create new identification number that can be used to join genetics data to 
# bridge data (bug in the database identification table)
newID <- chinIPES$UNIVERSAL_FISH_LABEL %>% 
  as.vector() %>% 
  strsplit(., split = "-") %>% 
  unlist() %>%
  matrix(., nrow = 5, ncol = length(chinIPES$UNIVERSAL_FISH_LABEL)) %>%
  t() %>%
  data.frame() %>%
  dplyr::rename("prog" = X1, "survey" = X2, "event" = X3, "species" = X4,
                "conFish" = X5) %>% 
  mutate(UNIVERSAL_FISH_LABEL = chinIPES$UNIVERSAL_FISH_LABEL,
         SPECIMEN_ID = chinIPES$SPECIMEN_ID) %>% 
  filter(species %in% c("124", "124J", "124A")) %>% 
  mutate(
    year = as.numeric(
      case_when(
        grepl("2017", prog) ~ "2017",
        grepl("2018", prog) ~ "2018",
        grepl("2019", prog) ~ "2019")),
    prog = case_when(
      grepl("IPES", prog) ~ "IPES",
      TRUE ~ "NA")
  ) %>%
  group_by(survey, event, year) %>%
  mutate(conFish = case_when(
    # 2017 data are not consecutive but 2018/2019 are; for 2017 make consecutive
    # for both pad with 0s
    as.character(year) == "2017" ~ row_number() %>%
           formatC(., width = 3, format = "d", flag = "0"),
    TRUE ~ formatC(conFish, width = 3, format = "d", flag = "0"))
  ) %>%
  ungroup()

# Join genetics data to new identification numbers, then to associated bridge
# data
dnaIPESOut <- chinIPES %>% 
  left_join(., newID, by = c("SPECIMEN_ID", "UNIVERSAL_FISH_LABEL")) %>% 
  mutate(
    #redefine age based on length (misentered by genetics lab)
    age = case_when(
      LENGTH <= 300 ~ "J",
      LENGTH > 300 ~ "A",
      TRUE ~ "NA"
    ),
    species = paste("124", age, sep = ""),
    BCSI_FISH_NUMBER = paste(year, survey, sep = "") %>% 
      paste(prog, ., event, species, conFish, sep = "-")
  ) %>% 
  left_join(dnaIPESTrim, ., by = "BCSI_FISH_NUMBER") %>% 
  #add lat longs and dates from ipes bridge data
  left_join(., 
            bridgeOut %>% 
              filter(dataset == "IPES") %>% 
              rename(BRIDGE_LOG_ID = bridge_log_id),
            by = c("BRIDGE_LOG_ID", "year")) %>% 
  select(fish_number = BCSI_FISH_NUMBER, station_id, dataset, year, month, jday, 
         start_lat, start_long, age, 
         ship_fl = LENGTH, BATCH_DNA_NUMBER,  
         STOCK_1:REGION_5, PROB_1:PROB_5) %>% 
  rename_all(., tolower) 

# Check whether there are fish with genetics data that are missing bridge data
# missed <- dnaIPESOut %>%
#   filter(is.na(station_id)) %>%
#   select(fish_number) %>%
#   as.vector()

# Merge high seas and genetics data
gsiBind <- dnaHSOut %>% 
  rbind(., dnaIPESOut) 

# Convert lat/long to utm
coords <- gsiBind %>% 
  select(fish_number, yUTM_start = start_lat, xUTM_start = start_long)
sp::coordinates(coords) <- c("xUTM_start", "yUTM_start")
sp::proj4string(coords) <- sp::CRS("+proj=longlat +datum=WGS84")
#SE Van Island extends into zone 10; not sure if that's an issue or not
coords2 <- sp::spTransform(coords, sp::CRS("+proj=utm +zone=9 ellps=WGS84")) %>%
  as(., "SpatialPoints")

# Convert all stocks to upper case
ups <- lapply(gsiBind %>% 
                select(stock_1:region_5),  
              toupper) %>% 
  data.frame() 

gsiOut <- gsiBind %>%
  select(-(stock_1:region_5)) %>% 
  cbind(., ups) %>% 
  cbind(., coords2@coords) %>% 
  select(fish_number:start_long, xUTM_start, yUTM_start, age:region_5)

# Add in regional roll ups that approximate what's used by CTC
gReg <- gsiOut %>% 
  select(fish_number, region_1:region_5) %>% 
  gather(key = "region_rank", value = "region", -fish_number) %>% 
  arrange(fish_number)
gStocks <- gsiOut %>% 
  select(fish_number, stock_1:stock_5) %>% 
  gather(key = "stock_rank", value = "stock", -fish_number) %>% 
  arrange(fish_number)
gsiLong <- gsiOut %>% 
  select(-(stock_1:region_5)) %>% 
  gather(key = "prob_rank", value = "prob", prob_1:prob_5) %>%
  arrange(fish_number) %>% 
  cbind(., gStocks[, -1], gReg[, -1]) %>% 
  # correct most ambiguous misspelled stocks
  mutate(
    stock = case_when(
      stock %in% c("BIG", "BIGQUL@LANG") ~ "BIG_QUALICUM",
      TRUE ~ as.character(stock)
  )) %>% 
  filter(!is.na(stock),
         #exclude individuals from bad sets based on comments
         !station_id %in% excludeTowsComments$STATION_ID)

# export distinct stocks to make key in makeStockKey.R
# change this to a sourced function?
# stockList <- gsiLong %>%
#   select(stock, region) %>%
#   distinct()
# saveRDS(stockList, here::here("data", "tempStockList.rds"))

# Import corrected stock list to calculate aggregate probabilities
cleanStockKey <- readRDS(here::here("data", "finalStockList.rds"))

gsiLongFull <- gsiLong %>% 
  select(-region_rank, - region) %>% 
  full_join(., cleanStockKey, by = "stock") 
saveRDS(gsiLongFull, here::here("data", "mergedGSI_highResLong.rds"))

# Use summed probabilities to focus on coarsest region 4 and assign individuals
# to single most likely region (as long as it exceeds 50%)
gsiLongAgg <- gsiLongFull %>% 
  group_by(fish_number, Region4Name) %>% 
  mutate(aggProb = sum(prob)) %>%
  group_by(fish_number) %>% 
  mutate(maxProb = max(aggProb)) %>% 
  select(-c(batch_dna_number:Region3Name)) %>%
  arrange(fish_number) %>% 
  ungroup() %>% 
  # arrange(year, month, jday, dataset, fish_number) %>% 
  distinct() %>%
  #remove probabilities that are not the max and individuals where max doesn't
  #exceed threshold (50%)
  filter(maxProb > 0.5,
         !aggProb < maxProb)


##### MERGE PROPORTIONS DATA #####

# Calculate catch proportions
stockComp <- gsiLongAgg %>% 
  select(-ship_fl, -c(aggProb:maxProb)) %>% 
  group_by(station_id, Region4Name, age) %>% 
  tally(name = "regCatch") %>% 
  group_by(station_id, age) %>% 
  mutate(samp_catch = sum(regCatch),
         regPpn = regCatch / samp_catch) %>% 
  pivot_wider(., names_from = Region4Name, values_from = regPpn) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>% 
  select(station_id, age, samp_catch:SEAK) 

# Add bridge data (which contains 0 catches) to juv and adult dataset separately
adultsOut <- stockComp %>% 
  filter(age == "A") %>%
  left_join(bridgeOut %>% 
              select(station_id, stableStation:head_depth, ck_adult),
            ., 
            by = "station_id") %>%
  mutate(
    #replace erroneous 0 catches (i.e. when bridge = 0 but fish were GSId)
    ck_adult = case_when(
      samp_catch > ck_adult ~ samp_catch,
      TRUE ~ ck_adult
    ),
    samp_ppn = case_when(
      samp_catch > "0" ~ samp_catch / ck_adult,
      samp_catch == "0" ~ 0)
  ) %>% 
  select(station_id:year, stableStation, start_lat:ck_adult, samp_catch, 
         samp_ppn, SalSea:SEAK)
juvOut <- stockComp %>% 
  filter(age == "J") %>% 
  left_join(bridgeOut %>% 
              select(station_id, stableStation:head_depth, ck_juv),
            ., 
            by = "station_id") %>%
  mutate(
    #replace erroneous 0 catches (i.e. when bridge = 0 but fish were GSId)
    ck_juv = case_when(
      samp_catch > ck_juv ~ samp_catch,
      TRUE ~ ck_juv
    ),
    samp_ppn = case_when(
      samp_catch > "0" ~ samp_catch / ck_juv,
      samp_catch == "0" ~ 0)
  ) %>% 
  select(station_id:year, stableStation, start_lat:ck_juv, samp_catch,
         samp_ppn, SalSea:SEAK)


saveRDS(adultsOut, here::here("data", "adultCatchGSI_reg4.rds"))
saveRDS(juvOut, here::here("data", "juvCatchGSI_reg4.rds"))
