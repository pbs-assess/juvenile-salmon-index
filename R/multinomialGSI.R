## Fit GSI data to multinomial model
# Nov. 28, 2019
# Note: Currently assumes certain assignment

library(tidyverse)

# CLEAN DATA -------------------------------------------------------------------

# identify focal stations from different dataset (FIX EVENTUALLY)
juv <- readRDS(here::here("data", "juvCatchGSI_reg4.rds")) %>% 
  filter(stableStation == "Y") %>% 
  mutate(week = lubridate::week(date))
gsi_long_agg <- readRDS(here::here("data", "longGSI_reg4.rds")) %>% 
  filter(age == "J", 
         station_id %in% juv$station_id) %>% 
  mutate(jdayZ = as.vector(scale(jday)[,1]),
         present = 1,
         #consolidate northern aggregates because rel. rare
         agg = case_when(
           Region4Name %in% c("NBC", "SEAK", "CoastUS") ~ "Other",
           TRUE ~ Region4Name),
         season = as.factor(case_when(
             month %in% c("2", "3") ~ "winter",
             month %in% c("5", "6", "7", "8") ~ "summer",
             month %in% c("9", "10", "11" , "12") ~ "fall")),
         year = as.factor(year)
  ) %>% 
  select(-Region4Name, -ship_fl, -c(xUTM_start:age), -c(aggProb:maxProb)) %>%
  mutate(season = fct_relevel(season, "fall", after = 1)) %>% 
  left_join(., juv %>% select(station_id, week), by = "station_id") %>% 
  distinct()

year_aggs <- expand.grid(year = unique(gsi_long_agg$year), 
                         agg = unique(gsi_long_agg$agg), 
                         present = 1)

## Replace year/aggregate combinations with 0 catches with 1s and spread to wide
# format
gsi_wide <- full_join(gsi_long_agg, year_aggs, 
                      by = c("year", "agg", "present")) %>% 
  pivot_wider(., names_from = agg, values_from = present) %>%
  mutate_if(is.numeric, ~replace_na(., 0)) 

glimpse(gsi_wide)

# PRELIMINARY VIS --------------------------------------------------------------

comp <- gsi_long_agg %>% 
  group_by(agg, season, year) %>% 
  tally() %>% 
  group_by(season, year) %>% 
  mutate(total = sum(n), 
         prop = n / total)

ggplot(comp) +
  geom_bar(aes(x = as.factor(year), y = prop, fill = agg), 
           stat = "identity") +
  facet_wrap(~season) +
  ggsidekick::theme_sleek()

comp_week <- gsi_long_agg %>% 
  group_by(agg, week) %>% 
  tally() %>% 
  group_by(week) %>% 
  mutate(total = sum(n), 
         prop = n / total)

ggplot(comp_week) +
  geom_bar(aes(x = as.factor(week), y = prop, fill = agg), 
           stat = "identity") 

# FIT MODEL --------------------------------------------------------------------

library(TMB)
compile("R/multinomialPractice/multinomial_generic.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic"))

## Data and parameters
year_vec <- seq(from = "2000", to = "2009", by = 1)
dum <- gsi_wide %>% 
  filter(year %in% year_vec) %>%
  mutate(year = factor(year))

# gsi_long_agg %>% 
#   filter(year %in% dum$year) %>% 
#   group_by(year, agg) %>% 
#   tally()

y_obs <- dum  %>% 
  select(Other:WCVI) %>% 
  as.matrix()

# X <- dum$jdayZ

.X <- model.matrix(~ year, dum)
# .X <- cbind(1, X) #predictor with intercept
# .X <- as.matrix(rep(1, length = nrow(y_obs)), ncol = 1) #int. only predictor

data <- list(cov = .X, y_obs = y_obs)
parameters <- list(betas = matrix(data = 0, nrow = ncol(.X), 
                                  ncol = ncol(y_obs) - 1))

## Make a function object
obj <- MakeADFun(data, parameters, DLL="multinomial_generic")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

ssdr <- summary(sdr)
ssdr

r <- obj$report()
r$probs
r$log_odds
r$logit_probs

# PLOT PREDICTIONS -------------------------------------------------------------
k <- ncol(y_obs) # number of stocks
stk_names <- colnames(y_obs)
N <- nrow(y_obs)

logit_probs_mat <- ssdr[rownames(ssdr) %in% "logit_probs", ] 
pred_ci <- data.frame(stock = rep(stk_names, each = N),
                      season = rep(dum$season, times = k),
                      year = as.factor(rep(dum$year, times = k)),
                      logit_prob_est = logit_probs_mat[ , "Estimate"],
                      logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est + (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est + (qnorm(0.975) * logit_prob_se))
         ) %>%
  distinct()

ggplot(pred_ci) +
  # geom_point(aes(x = year, y = pred_prob)) +
  geom_pointrange(aes(x = year, y = pred_prob, ymin = pred_prob_low, 
                      ymax = pred_prob_up)) +
  # geom_ribbon(aes(x = jday, ymin = pred_prob_low, ymax = pred_prob_up), 
  #             fill = "#bfd3e6") +
  # geom_line(aes(x = jday, y = pred_prob), col = "#810f7c", size = 1) +
  facet_grid( ~ stock) +
  labs(y = "Probability", x = "Year")

