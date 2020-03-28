## IPES and HS Merge
# Data import/clean portion of script written by E. Anderson to merge historical
# high seas and contemporary IPES data
# Original script: 
# https://github.com/ErikaDAnderson/SOPO/blob/master/RScripts/SOPO2019.R


#####################################
# load IPES data
#####################################

# CPUE data
# load as csv file since view built on views
# use for swept volume and join to catch 
# use BRIDGE_FIELD_ID becuase different database version

# March 2020 adjustments to view based on target depth averages instead of overall averages
# from IPES_TrawlBC_v20.02b database 
volume_ipes_orig <- read_csv("Input/2019/JB_VIEW_IPES_CPUE_BRIDGE_LOG_FIELD_ID.csv")

# estalish connection to IPES Access database
# used IPES Report Version since has extra GSI table that I included
db_ipes <- "C:/Users/andersoned/Documents/GitHub/IPES_Report/Input/2019/IPES_TrawlDB_v19.07f_2017_18_19.mdb"
myconn_ipes <- odbcConnectAccess2007(db_ipes)

tows_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, BRIDGE_LOG.EVENT_TYPE, 
BRIDGE_LOG.START_LATITUDE, BRIDGE_LOG.START_LONGITUDE, 
IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, 
BRIDGE_LOG.BLOCK_DESIGNATION, BRIDGE_LOG.STRATUM, BRIDGE_LOG.USABILITY_CODE
FROM TRIP LEFT JOIN BRIDGE_LOG ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND 
((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 
Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') 
AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((BRIDGE_LOG.USABILITY_CODE)<>5));
                          ")
# limited to daylight, usable tows within IPES area
cpue_ipes_orig <- sqlQuery(myconn_ipes, "SELECT BRIDGE_LOG.BRIDGE_LOG_ID, BRIDGE_LOG.BRIDGE_LOG_FIELD_ID, 
TRIP.TRIP_NAME, BRIDGE_LOG.EVENT_NUMBER, BRIDGE_LOG.EVENT_TYPE, BRIDGE_LOG.START_LATITUDE, 
BRIDGE_LOG.START_LONGITUDE, BRIDGE_LOG.BLOCK_DESIGNATION, BRIDGE_LOG.STRATUM, CATCH.SPECIES_CODE, 
CATCH.JUVENILE_CATCH_COUNT, IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, BRIDGE_LOG.USABILITY_CODE
FROM TRIP INNER JOIN (BRIDGE_LOG LEFT JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) 
ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND 
((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND 
((BRIDGE_LOG.USABILITY_CODE)<>5));
       ")

# load bio data for length weight 
lw_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, SPECIMEN.UNIVERSAL_FISH_LABEL, 
CATCH.SPECIES_CODE, SPECIMEN.LENGTH, SPECIMEN.WEIGHT
FROM TRIP LEFT JOIN ((BRIDGE_LOG LEFT JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) 
LEFT JOIN SPECIMEN ON CATCH.CATCH_ID = SPECIMEN.CATCH_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((CATCH.SPECIES_CODE)='108' Or (CATCH.SPECIES_CODE)='112' Or (CATCH.SPECIES_CODE)='115' Or 
(CATCH.SPECIES_CODE)='118' Or (CATCH.SPECIES_CODE)='124') AND ((BRIDGE_LOG.EVENT_TYPE)='midwater tow') 
AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND 
((BRIDGE_LOG.USABILITY_CODE)<>5));
                         ")

# GSI data from IPES for 2019
# IPES survey blocks only
# juvenile defined by length of coho, sockeye, chinook only
# daytime usable tows
gsi_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, SPECIMEN_COLLECTED.COLLECTED_ATTRIBUTE_CODE, 
SPECIMEN_COLLECTED.STORAGE_CONTAINER_SUB_ID, CATCH.SPECIES_CODE, BRIDGE_LOG.BRIDGE_LOG_ID, 
BRIDGE_LOG.BRIDGE_LOG_FIELD_ID, BRIDGE_LOG.TRIP_ID, BRIDGE_LOG.EVENT_DATE, BRIDGE_LOG.BLOCK_DESIGNATION, 
BRIDGE_LOG.USABILITY_CODE, IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, SPECIMEN.LENGTH
FROM TRIP LEFT JOIN (BRIDGE_LOG LEFT JOIN (CATCH LEFT JOIN (SPECIMEN LEFT JOIN SPECIMEN_COLLECTED ON 
SPECIMEN.SPECIMEN_ID = SPECIMEN_COLLECTED.SPECIMEN_ID) ON CATCH.CATCH_ID = SPECIMEN.CATCH_ID) ON 
BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((TRIP.TRIP_YEAR)=2019) AND ((SPECIMEN_COLLECTED.COLLECTED_ATTRIBUTE_CODE)=4) AND 
((CATCH.SPECIES_CODE)='115' Or (CATCH.SPECIES_CODE)='118' Or (CATCH.SPECIES_CODE)='124') AND 
((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((BRIDGE_LOG.USABILITY_CODE)<>5) AND 
((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND ((SPECIMEN.LENGTH)<350));
                          ")

# close database
close(myconn_ipes)

#####################################
# wrangle IPES data
#####################################

# find distinct tows for IPES data
tows_ipes <-  tows_ipes_orig %>%
  dplyr::select(TRIP_YEAR, BLOCK_DESIGNATION, START_LATITUDE, START_LONGITUDE)


# limit swept volume columns
volume_ipes <- volume_ipes_orig %>%
  filter(BLOCK_DESIGNATION > 0) %>%
  filter(DayNight == "Day") %>%
  filter(USABILITY_CODE == 1) %>%
  distinct(BRIDGE_LOG_ID, BRIDGE_LOG_FIELD_ID, TRIP_NAME, EVENT_NUMBER, EVENT_TYPE, BLOCK_DESIGNATION, STRATUM,
           OfficialVolumeSwept_km3)

# select salmon 
# calculate CPUE by swept volume
cpue_ipes_salmon <- cpue_ipes_orig %>%
  filter(SPECIES_CODE %in% c(108, 112, 115, 118, 124)) %>%
  left_join(., volume_ipes, by = c("BRIDGE_LOG_FIELD_ID", "TRIP_NAME", "EVENT_NUMBER", 
                                   "EVENT_TYPE", "BLOCK_DESIGNATION", "STRATUM")) %>%
  mutate(EVENT = str_c(TRIP_NAME, str_pad(EVENT_NUMBER, 3, side = "left", pad = 0), sep = "-"),
         CPUE = JUVENILE_CATCH_COUNT, 
         logCPUE1 = log(CPUE + 1),
         TRIP_YEAR = as.numeric(str_extract(TRIP_NAME, "[0-9]+"))) %>%
  dplyr::select(TRIP_YEAR, BLOCK_DESIGNATION, EVENT, START_LATITUDE, START_LONGITUDE, JUVENILE_CATCH_COUNT, 
                OfficialVolumeSwept_km3, SPECIES_CODE, CPUE, logCPUE1) %>%
  rename(SWEPT_VOLUME = OfficialVolumeSwept_km3)

#####################################
# load high seas data
#####################################
# estalish connection to high sea Access database
db_hs <- "C:/Users/andersoned/Documents/BCSI/High Seas Salmon Database/HSSALMON.accdb"
myconn_hs <- odbcConnectAccess2007(db_hs)

# get bridge data from high seas
# synoptic stations only
# June and July
# headrope depth <20 m
cpue_hs_orig <- sqlQuery(myconn_hs, "SELECT STATION_INFO.CRUISE, STATION_INFO.STATION_ID, 
STATION_INFO.REGION, STATION_INFO.REGION_CODE, STATION_INFO.SYNOPTIC_STATION, BRIDGE.Year, 
BRIDGE.Month, BRIDGE.START_LAT, BRIDGE.START_LONG, BRIDGE.DISTANCE, BRIDGE.START_BOT_DEPTH, 
BRIDGE.END_BOT_DEPTH, BRIDGE.NET_OPENING_WIDTH, BRIDGE.NET_OPENING_HEIGHT, BRIDGE.HEAD_DEPTH, 
BRIDGE.PK_JUV, BRIDGE.CM_JUV, BRIDGE.SE_JUV, BRIDGE.CO_JUV, BRIDGE.CK_JUV
FROM STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.Month)='JUN' Or (BRIDGE.Month)='JUL') 
AND ((BRIDGE.HEAD_DEPTH)<20));
                         ")

# get length weight data from high seas
# synoptic stations
# June and July
# headrope depth < 20 m
lw_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.Year, BIOLOGICAL_JUNCTION.FISH_NUMBER, BIOLOGICAL.SPECIES, 
BIOLOGICAL.SHIP_FL, BIOLOGICAL.SHIP_WT
FROM STATION_INFO INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON 
BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN BIOLOGICAL ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON 
(STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) AND (STATION_INFO.STATION_ID = BRIDGE.STATION_ID)
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.Month)='JUN' Or (BRIDGE.Month)='JUL') AND 
(BRIDGE.HEAD_DEPTH)<20);
                       ")

# get calorimetry data from high seas
cal_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.YEAR, BIOLOGICAL_JUNCTION.FISH_NUMBER, 
CALORIMETRY.HEAT_RELEASED_CAL, CALORIMETRY.HEAT_RELEASED_KJ, CALORIMETRY.DUPLICATE, 
CALORIMETRY.DATA_ISSUE
FROM STATION_INFO INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON 
BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN CALORIMETRY ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = CALORIMETRY.FISH_NUMBER) ON 
(STATION_INFO.STATION_ID = BRIDGE.STATION_ID) AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((CALORIMETRY.DATA_ISSUE)='N') AND ((STATION_INFO.SYNOPTIC_STATION)=True) AND 
((BRIDGE.MONTH)='JUN' Or (BRIDGE.MONTH)='JUL') AND (BRIDGE.HEAD_DEPTH)<20);
                        ")

# get calorimetry data for IPES from high seas table
# need to limit using bridge info to IPES blocks and usable tows
cal_ipes_orig <- sqlQuery(myconn_hs, "SELECT CALORIMETRY_IPES.FISH_NUMBER, CALORIMETRY_IPES.DUPLICATE, 
CALORIMETRY_IPES.HEAT_RELEASED_CAL, CALORIMETRY_IPES.HEAT_RELEASED_KJ, CALORIMETRY_IPES.DATA_ISSUE, 
CALORIMETRY_IPES.COMMENTS
FROM CALORIMETRY_IPES
WHERE (((CALORIMETRY_IPES.DATA_ISSUE)='N'));
                          ")

# get calorimetry data for hisoric data
# limited to synoptic stations
# June and July
# headrope depth <20 m
cal_historic_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.YEAR, BRIDGE.MONTH, BIOLOGICAL.SPECIES, 
BIOLOGICAL_JUNCTION.FISH_NUMBER, PROXIMATE_FISH.ENERGY_BOMB, PROXIMATE_FISH.ENERGY_BOMB_BLIND_DUPL
FROM STATION_INFO INNER JOIN (((BIOLOGICAL_JUNCTION INNER JOIN PROXIMATE_FISH ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = PROXIMATE_FISH.FISH_NUMBER) INNER JOIN BRIDGE ON 
BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE.STATION_ID) INNER JOIN BIOLOGICAL ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON (STATION_INFO.STATION_ID = BRIDGE.STATION_ID) 
AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((BRIDGE.MONTH)='JUN' Or (BRIDGE.MONTH)='JUL') AND ((STATION_INFO.SYNOPTIC_STATION)=True) AND 
(BRIDGE.HEAD_DEPTH)<20);
                          ")

# close database
close(myconn_hs)

#####################################
# wrangle high seas data
#####################################

# check for empty net dimensions and distance values
missing_hs <- cpue_hs_orig %>%
  filter(is.na(NET_OPENING_WIDTH)| is.na(NET_OPENING_HEIGHT) | is.na(DISTANCE))

# pull vector of cruises with missing info
missing_hs_vec <- missing_hs %>%
  pull(CRUISE)

#### all missing from cruise 201893 
# no missing distance values

# caluculate net averages for cruises at specific headrope depths
hs_net_avg <- cpue_hs_orig %>%
  group_by(CRUISE, HEAD_DEPTH) %>%
  summarize(AvgNET_OPENING_WIDTH = mean(NET_OPENING_WIDTH, na.rm = TRUE),
            AvgNET_OPENING_HEIGHT = mean(NET_OPENING_HEIGHT, na.rm = TRUE)) %>%
  filter(CRUISE %in% missing_hs_vec)

# use Sea Crest average net width and height from whole survey
# Jackie did this in past with net height = 14 and net width = 33 m
# gear comparison study was used initially 
# same vessel, Captain and net but the chain links and floats were different
# use average of 10.6 for headrope ~ 15 m


# replace missing net values
cpue_hs_net <- cpue_hs_orig %>%
  mutate(NET_OPENING_WIDTH = case_when(
    !(is.na(NET_OPENING_WIDTH)) ~ NET_OPENING_WIDTH,
    is.na(NET_OPENING_WIDTH) ~ 33), # from average for 201893 survey
    NET_OPENING_HEIGHT = case_when(
      !(is.na(NET_OPENING_HEIGHT)) ~ NET_OPENING_HEIGHT,
      is.na(NET_OPENING_HEIGHT) & HEAD_DEPTH <= 7~ 14, # from average for 201893 survey
      is.na(NET_OPENING_HEIGHT) & HEAD_DEPTH > 7 ~ 11), # average from same cruise rounded
    NET_AREA_KM = (NET_OPENING_WIDTH/1000) * (NET_OPENING_HEIGHT/1000),
    DISTANCE_KM = DISTANCE * 1.852,
    SWEPT_VOLUME = NET_AREA_KM*DISTANCE_KM)


# create function to rearrange format to match IPES for each salmon
wrangle_fn <- function(df, colName, speciesCode) {
  
  df %>%
    mutate(BLOCK_DESIGNATION = NA) %>%
    dplyr::select(Year, BLOCK_DESIGNATION, STATION_ID, START_LAT, START_LONG, colName, SWEPT_VOLUME) %>%
    mutate(SPECIES_CODE = speciesCode) %>%
    rename(TRIP_YEAR = Year,
           EVENT = STATION_ID,
           START_LATITUDE = START_LAT,
           START_LONGITUDE = START_LONG,
           JUVENILE_CATCH_COUNT = colName)
}

# apply to all salmon species  
cpue_hs_pk <- wrangle_fn(cpue_hs_net, "PK_JUV", 108)
cpue_hs_ck <- wrangle_fn(cpue_hs_net, "CK_JUV", 124)
cpue_hs_co <- wrangle_fn(cpue_hs_net, "CO_JUV", 115)
cpue_hs_cm <- wrangle_fn(cpue_hs_net, "CM_JUV", 112)
cpue_hs_se <- wrangle_fn(cpue_hs_net, "SE_JUV", 118)

# bind species back together
cpue_hs <- rbind(cpue_hs_ck, cpue_hs_cm, cpue_hs_co, cpue_hs_pk, cpue_hs_se)

# make CPUE by swept volume and log(cpue + 1)
cpue_hs <- cpue_hs %>%
  mutate(CPUE = JUVENILE_CATCH_COUNT / SWEPT_VOLUME,
         logCPUE1 = log(CPUE + 1)) # natural log by default
