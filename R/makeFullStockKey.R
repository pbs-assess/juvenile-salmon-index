## Short script to build a more complete stock key
# October 18, 2019

library(tidyverse)

scStockKey <- read.csv(here::here("data", "stockKeys", 
                                  "southCoastStockKey.csv")) %>%
    mutate(stock = toupper(Stock)) %>%
    select(stock, Region1Name, Region2Name, Region3Name)
# list of observed stocks from high seas
stockKeyHS <- readRDS(here::here("data", "stockKeys", "tempStockList.rds")) %>%
  left_join(., scStockKey, by = "stock") %>% 
  mutate(Region1Name = as.character(Region1Name),
         Region2Name = as.character(Region2Name),
         Region3Name = as.character(Region3Name)) %>% 
  distinct()
# list of observed stocks from WCVI troll
stockKey1 <- readRDS(here::here("data", "stockKeys", "wcviTrollStocks.rds")) %>% 
  full_join(., stockKeyHS, by = "stock") %>% 
  distinct()

# Associate misspelled and unknown stocks with higher level regions
stockKeyOut <- stockKey1 %>% 
    select(stock:Region1Name) %>%
    mutate(
      #add unknown stocks
      Region1Name = case_when(
        grepl("PIT", stock) ~ "LWFR-Su",
        stock == "SOLDUC_F" ~ "Washington_Coast",
        stock == "WALKER" ~ "UPFR",
        stock == "HOMATHKO" ~ "SOMN",
        stock == "SALOOMPT" ~ "NOMN",
        grepl("SNOHOMISH", stock) ~ "N_Puget_Sound",
        stock == "LYON'S_FERRY_F" ~ "Snake_R_fa",
        stock == "L_SHUS@U_ADAMS" ~ "SOTH",
        grepl("SKAGIT", stock) ~ "N_Puget_Sound",
        grepl("UMPQUA", stock) ~ "Mid_Oregon_Coast",
        grepl("MAMQUAM", stock) ~ "SOMN",
        grepl("SQUAMISH", stock) ~ "SOMN",
        grepl("SOOS_CR", stock) ~ "S_Puget_Sound",
        grepl("GREEN", stock) ~ "S_Puget_Sound",
        grepl("WILLA", stock) ~ "Washington_Coast",
        grepl("SKYKOMISH", stock) ~ "C_Puget_Sound",
        stock == "BUTE" ~ "SOMN",
        grepl("SILMIL", stock) ~ "U_Columbia_R_su/fa",
        grepl("HOH_RI", stock) ~ "Washington_Coast",
        grepl("CLE_ELM", stock) ~ "S_Puget_Sound",
        grepl("QUINA", stock) ~ "Washington_Coast",
        grepl("ABERN", stock) ~ "L_Columbia_R_fa",
        grepl("WOSS", stock) ~ "ECVI",
        grepl("QUATSE", stock) ~ "ECVI",
        grepl("HORSEY", stock) ~ "UPFR",
        grepl("HORSEF", stock) ~ "MUFR",
        grepl("OKANAG", stock) ~ "U_Columbia_R_su/fa",
        grepl("WENATCH", stock) ~ "U_Columbia_R_su/fa",
        grepl("COLE", stock) ~ "N_California/S_Oregon_Coast",
        grepl("KWIN", stock) ~ "Nass",
        grepl("WANN", stock) ~ "NOMN",
        grepl("HIRS", stock) ~ "NOMN",
        grepl("DEAN", stock) ~ "NOMN",
        grepl("CRAIG", stock) ~ "Stikine",
        grepl("DUDI", stock) ~ "Taku",
        grepl("ANDR", stock) ~ "SSE_Alaska_Stikine_R",
        grepl("SWEET", stock) ~ "Skeena Mid",
        grepl("MORI", stock) ~ "Skeena Bulkley",
        grepl("L_KAL", stock) ~ "Skeena Lower",
        grepl("DESC", stock) ~ "U_Columbia_R_su/fa",
        grepl("LITTLEC", stock) ~ "N_Puget_Sound",
        grepl("STILLAG", stock) ~ "C_Puget_Sound",
        grepl("TYNE", stock) ~ "N_Puget_Sound",
        grepl("SERP", stock) ~ "LWFR-F",
        grepl("ELK", stock) ~ "MUFR",
        grepl("CHILL", stock) ~ "LWFR-F",
        grepl("BIG_Q", stock) ~ "ECVI",
        grepl("EUCH", stock) ~ "N_Oregon_Coast",
        grepl("QUEE", stock) ~ "Washington_Coast",
        grepl("WHITE", stock) ~ "S_Puget_Sound",
        grepl("NOOK", stock) ~ "N_Puget_Sound",
        grepl("SLOQUET", stock) ~ "LWFR-Sp",
        grepl("RAFT", stock) ~ "NOTH",
        grepl("PUNT", stock) ~ "ECVI",
        grepl("SHAK", stock) ~ "Stikine",
        grepl("LEMI", stock) ~ "NOTH",
        grepl("SIUS", stock) ~ "Mid_Oregon_Coast",
        grepl("TOUL", stock) ~ "Central_Valley_fa",
        grepl("ENTI", stock) ~ "Mid_and_Upper_Columbia_R_sp",
        grepl("STAN", stock) ~ "Central_Valley_fa",
        grepl("TAKIA", stock) ~ "NOMN",
        grepl("PORTE", stock) ~ "SOMN",
        grepl("BULK", stock) ~ "Skeena Bulkley",
        grepl("ELWHA", stock) ~ "Washington_Coast",
        grepl("NANA", stock) ~ "ECVI",
        stock == "LITTLE" ~ "SOTH",
        grepl("TRASK", stock) ~ "N_Oregon_Coast",
        grepl("CHEHA", stock) ~ "UPFR",
        grepl("EAGL", stock) ~ "SOTH",
        grepl("NEST", stock) ~ "N_Oregon_Coast",
        grepl("COWEE", stock) ~ "L_Columbia_R_fa",
        grepl("LOBST", stock) ~ "N_California/S_Oregon_Coast",
        grepl("PIST", stock) ~ "N_California/S_Oregon_Coast",
        grepl("HUNTER", stock) ~ "N_California/S_Oregon_Coast",
        grepl("THOMP", stock) ~ "SOTH",
        grepl("SANDY", stock) ~ "L_Columbia_R_fa",
        grepl("KENNET", stock) ~ "UPFR",
        grepl("WILLOW", stock) ~ "UPFR",
        grepl("SLIM", stock) ~ "UPFR",
        grepl("YUBA", stock) ~ "Central_Valley_fa",
        grepl("NEHA", stock) ~ "N_Oregon_Coast",
        grepl("WINCHUK", stock) ~ "N_California/S_Oregon_Coast",
        grepl("CLEEL", stock) ~ "Mid_Oregon_Coast",
        grepl("TRAPP", stock) ~ "Taku_R",
        stock == "EEL_F" ~ "California_Coast",
        grepl("SQUIN", stock) ~ "Skeena Upper",
        grepl("ATN", stock) ~ "NOMN",
        grepl("CHEWU", stock) ~ "Mid_and_Upper_Columbia_R_sp",
        grepl("CHIW", stock) ~ "Mid_and_Upper_Columbia_R_sp",
        grepl("NEVI", stock) ~ "UPFR",
        grepl("PTARM", stock) ~ "UPFR",
        grepl("NAHAT", stock) ~ "LWFR-F",
        grepl("TUY", stock) ~ "Stikine",
        region == "WCVI" ~ "WCVI",
        region == "SKEENA MID" ~ "Skeena Mid",
        region == "SKEENA LOWER" ~ "Skeena Lower",
        region == "MUFR" ~ "MUFR",
        region == "Alaska" ~ "SSE_Alaska",
        grepl("COLUMBIA-SP", region) ~ "Mid_and_Upper_Columbia_R_sp",
        grepl("UPPER COLUMBIA-SU", region) ~ "U_Columbia_R_su/fa",
        grepl("TRINITY", region) ~ "Klamath_R",
        region == "SNAKE-SP/SU" ~ "Snake_R_sp/su",
        region == "MID COL-SP" ~ "Mid_and_Upper_Columbia_R_sp",
        region == "UPPER WILLAMETTE" ~ "L_Columbia_R_fa",
        region == "UP WILLAMETTE" ~ "L_Columbia_R_fa",
        region == "CENTRAL VALLEY-F" ~ "Central_Valley_fa",
        region == "CENT VAL-F" ~ "Central_Valley_fa",
        region == "CENTRAL VALLEY-SP" ~ "Central_Valley_sp",
        TRUE ~ as.character(Region1Name)
      )
    ) %>%
  select(-region) %>% 
  distinct() %>% 
  # glimpse()
  left_join(., 
            stockKey1 %>% 
              select(Region1Name:Region3Name) %>% 
              distinct(), 
            by = "Region1Name") %>% 
  distinct() %>% 
  # glimpse()
  # add higher level regional aggregates
  mutate(Region2Name =
           case_when(
             stock %in% c("SKAGIT_SU", "SKYKOMISH_SU") ~ "Puget Sound Summer",
             stock == "NOOKSACK_SP@KE" ~ "Puget Sound Spring",
             Region1Name == "N_Puget_Sound" ~ "North Puget Sound Fall",
             Region1Name == "C_Puget_Sound" ~ "Central Puget Sound Fall",
             Region1Name == "S_Puget_Sound" ~ "South Puget Sound Fall",
             grepl("GOLD", stock) ~ "West Coast Hatchery",
             grepl("TOQUA", stock) ~ "West Coast Hatchery",
             Region1Name %in% c("CONUMA", "ROBERTSON", "NITINAT", "THORNTON") ~ 
               "West Coast Hatchery",
             Region1Name == "WCVI" ~ "West Coast Wild",
             grepl("QUINS", stock) ~ "Upper Strait of Georgia",
             grepl("NIMP", stock) ~ "North/Central BC",
             Region1Name == "ECVI" ~ 
               "Lower Strait of Georgia/Lower GS Hatchery",
             stock %in% c("CAPILANO") ~ "Fraser Late",
             stock %in% c("BUTE", "DEVEREUX", "HOMATHKO", "KLINAKLINI") ~ 
               "Upper Strait of Georgia",
             Region1Name == "SOMN" ~ 
               "Lower Strait of Georgia/Lower GS Hatchery",
             Region1Name %in% c("Mid_Oregon_Coast", 
                                "N_California/S_Oregon_Coast",
                                "N_Oregon_Coast") ~ 
               "Oregon Coastal North Migrating",
             Region1Name == "Washington_Coast" ~ 
               "Washington Coast/Juan de Fuca",
             Region1Name == "Snake_R_fa" ~ "Snake Fall",
             Region1Name == "L_Columbia_R_fa" ~ "Fall Cowlitz/Willamette",
             Region1Name %in% c("Central_Valley_sp",
                                "Klamath_R",
                                "California_Coast") ~ "California",
             Region1Name == "Taku_R" ~ "Alaska South SE",
             Region2Name == "Califormia" ~ "California",
             TRUE ~ as.character(Region2Name)
           ),
         Region3Name = 
           case_when(
             grepl("Puget", Region1Name) ~ "Puget Sound",
             grepl("Oregon", Region1Name) ~ "Oregon/California",
             Region1Name == "Washington_Coast" ~ "Washington Coast",
             grepl("NIMP", stock) ~ "North/Central BC",
             Region1Name == "ECVI" ~ "SOG",
             stock %in% c("CAPILANO") ~ "Fraser River",
             grepl("CHEAK", stock) ~ "SOG",
             Region1Name == "SOMN" ~ "SOG",
             grepl("Snake", Region1Name) ~ "Snake",
             grepl("Columbia", Region1Name) ~ "Columbia",
             grepl("California", Region2Name) ~ "Oregon/California",
             Region2Name == "Alaska South SE" ~ "Alaska South SE",
             TRUE ~ as.character(Region3Name)
           ),
         #even higher level aggreagtes for preliminary modeling
         Region4Name =
           case_when(
             Region3Name %in% c("Washington Coast", "Oregon/California") ~ 
               "CoastUS",
             Region3Name %in% c("SOG", "Puget Sound", "Fraser River") ~ 
               "SalSea",
             Region3Name %in% c("Columbia", "Snake") ~ "ColR",
             Region3Name == "North/Central BC" ~ "NBC",
             Region3Name == "Alaska South SE" ~ "SEAK",
             TRUE ~ Region3Name
           )) %>% 
  distinct() %>% 
  arrange(Region4Name, Region3Name, Region2Name, Region1Name, stock)

# check for gaps
stockKeyOut %>%
  select(stock, Region1Name, Region3Name) %>%
  filter(is.na(Region1Name)) %>%
  distinct()
stockKeyOut %>%
  select(stock, Region1Name, Region3Name) %>%
  filter(is.na(Region3Name)) %>%
  distinct()
# stockKeyOut %>%
#   select(stock, Region1Name:Region3Name) %>%
#   filter(is.na(stock) | is.na(Region1Name) | is.na(Region2Name) |
#            is.na(Region3Name))

saveRDS(stockKeyOut, here::here("data", "finalStockList.rds"))
write.csv(stockKeyOut, here::here("data", "finalStockList.csv"),
          row.names = FALSE)

