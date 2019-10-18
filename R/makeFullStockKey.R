## Short script to build a more complete stock key
# October 18, 2019

scStocks <- read.csv(here::here("data", "southCoastStockKey.csv"))
stockKey <- readRDS(here::here("data", "tempStockList.rds"))

# Associate misspelled and unknown stocks with higher level regions
stockKey <- stockKey %>% 
  mutate(
    stock = case_when(
      stock == "BIGQUL@LANG" ~ "BIG_QUALICUM",
      TRUE ~ as.character(stock)
    ),
    Region1Name = case_when(
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
      grepl("SKYKOMISH", stock) ~ "N_Puget_Sound",
      stock == "BUTE" ~ "SOMN",
      grepl("SILMIL", stock) ~ "U_Columbia_R_su/fa",
      grepl("HOH_RI", stock) ~ "Washington_Coast",
      grepl("CLE_ELM", stock) ~ "S_Puget_Sound",
      grepl("QUINA", stock) ~ "Washington_Coast",
      grepl("ABERN", stock) ~ "L_Columbia_R_fa",
      grepl("WOSS", stock) ~ "ECVI",
      grepl("QUATSE", stock) ~ "ECVI",
      grepl("HORSE", stock) ~ "MUFR",
      grepl("OKANAG", stock) ~ "U_Columbia_R_su/fa",
      grepl("WENATCH", stock) ~ "U_Columbia_R_su/fa",
      grepl("COLE", stock) ~ "N_California/S_Oregon_Coast",
      grepl("KWIN", stock) ~ "NASS",
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
      grepl("STILLAG", stock) ~ "N_Puget_Sound",
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
      
      
      region == "WCVI" ~ "WCVI",
      region == "SKEENA MID" ~ "Skeena Mid",
      grepl("COLUMBIA-SP", region) ~ "Mid_and_Upper_Columbia_R_sp",
      region == "SNAKE-SP/SU" ~ "Snake_R_sp/su",
      region == "UPPER WILLAMETTE" ~ "L_Columbia_R_fa",
      region == "CENTRAL VALLEY-F" ~ "Central_Valley_fa",
      region == "CENTRAL VALLEY-SP" ~ "Central_Valley_sp",
      TRUE ~ as.character(Region1Name)
    )
  )

stockKey %>% 
  # filter(is.na(region)) %>% 
  select(stock, region, Region1Name) %>% 
  filter(is.na(Region1Name)) %>% 
  distinct() %>% 
  head()

stockKey %>% 
  select(stock, region, Region1Name) %>% 
  filter(grepl("SIUS", stock)) %>% 
  distinct()
