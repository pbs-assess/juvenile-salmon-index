# Simulate and fit unordered multinomial response models with addition of 
# random intercepts; extension of multinomial-regression_randomInt that also
# includes two additional fixed effects

library(tidyverse)

set.seed(42)
# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 250
sd_site <- 0.5
N <- n_sites * n_obs_per_site

# fixed intercepts
group_ints <- c(0.3, -1.4)
k <- length(group_ints) #number of groups - 1
reg2_ints <- c(1.5, -0.25) #i.e. how region 2 differs from reference for each group
reg3_ints <- c(0.27, 0.49) #i.e. how region 3 differs from reference for each group
fac2_ints <- c(0, -2) #i.e. how factor 2 differs from reference for each group
fix_ints <- matrix(c(group_ints, reg2_ints, reg3_ints, fac2_ints
                 ), ncol = k, byrow = T)

reg_dat <- data.frame(reg = as.factor(sample(c(1, 2, 3), size = N, replace = T)),
                      fac = as.factor(sample(c(1, 2), size = N, replace = T)))

# model matrix for fixed effects
fix_mm <- model.matrix(~ (reg + fac), reg_dat)

# function to simulate multinomial data and fit hierarchical models
f_sim <- function(trial = 1) {
  
  #calculate fixed effects using parameter matrix and model matrix
  n_fix_cov <- nrow(fix_ints)
  sum_fix_eff <- matrix(NA, nrow = N, ncol = k)
  for (kk in 1:k) {
    fix_eff <- matrix(nrow = N, ncol = n_fix_cov)
    for (i in 1:N) {
      for (j in 1:n_fix_cov) {
        # fix_eff[i, j] <- group_ints[kk] + (fix_ints[j, kk] * fix_mm[i, j])  
        fix_eff[i, j] <- fix_ints[j, kk] * fix_mm[i, j]  
      }
    }
    sum_fix_eff[ , kk] <- apply(fix_eff, 1, sum)
  }

  # generate random effects, then combine with fixed
  site_mean_a <- rnorm(mean = 0,
                       sd = sd_site,
                       n = n_sites)
  site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                       sd = 0, 
                       n = N)
  # site_mean_a <- runif(n = n_sites, min = -1, max = 1)
  # site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
  #                      sd = sd_site, 
  #                      n = N)
  site_dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1), 
                                         each = n_obs_per_site)),
                         site_mean = rep(site_mean_a, each = n_obs_per_site),
                         site_obs = site_obs_ia,
                         sd_site = sd_site
                    ) %>% 
    cbind(reg_dat) %>%  #bind regional fixed effects data 
    # arrange(site, reg) %>% 
    mutate(g1_fe = sum_fix_eff[, 1],
           g2_fe = sum_fix_eff[, 2])
  
  # generate log_odds based on fixed and random intercepts
  datf <- site_dat %>% 
    mutate(site = as.factor(site),
           log_odds_1_3 = g1_fe + site_obs,
           log_odds_2_3 = g2_fe + site_obs,
           # probability of category 1, 2, and 3:
           p1 = exp(log_odds_1_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
           p2 = exp(log_odds_2_3) / (1 + exp(log_odds_1_3) + exp(log_odds_2_3)),
           p3 = 1 - (p1 + p2)
    )
  
  # create the observations:
  y <- numeric(length = nrow(datf))
  p <- datf %>% 
    select(p1:p3) %>%
    as.matrix()
  for (i in seq_along(y)) {
    temp <- rmultinom(1, 1, p[i, ])
    y[i] <- which(temp == 1)
  }
  datf <- cbind(datf, y)
  # ggplot(datf, aes(x = y)) +
  #   geom_histogram() + #negative values in site 2 inc. probability of third cat
  #   ggsidekick::theme_sleek() +
    # facet_wrap(~reg, nrow = 2)
  # facet_wrap(reg~site, nrow = 3)
  
  # matrix version for dmultinom():
  Y <- matrix(ncol = 3, nrow = N, data = 0)
  for (i in seq_along(y)) Y[i, y[i]] <- 1
  
  # true values vector
  tv <- c(as.vector(fix_ints), log(sd_site), site_mean_a, sd_site)
  
  return(list("trial" = trial, "obs" = Y, "fixed_fac" = site_dat$reg, 
              "rand_fac" = datf$site, "trim_data" = tv, 
              "full_data" = datf))
}

dat_list <- vector(mode = "list", length = 100)
for (i in seq_along(dat_list)) {
  dat_list[[i]] <- f_sim(trial = i)
}


library(TMB)
# compile("R/sim_practice/multinomialPractice/multinomial_generic_randInt_fixInt.cpp")
# dyn.load(dynlib("R/sim_practice/multinomialPractice/multinomial_generic_randInt_fixInt"))
# version that integrates out random effects for predictions
# (otherwise same as above)
compile("src/multinomial_generic_randInt2.cpp")
dyn.load(dynlib("src/multinomial_generic_randInt2"))
# fixed effects version for comparison
# compile("R/sim_practice/multinomialPractice/multinomial_generic_fixInt.cpp")
# dyn.load(dynlib("R/sim_practice/multinomialPractice/multinomial_generic_fixInt"))

## Data and parameters
y_obs <- dat_list[[1]]$obs
#vector of random intercept ids
rfac <- as.numeric(dat_list[[1]]$rand_fac) - 1 #subtract for indexing by 0
fac_dat <- dat_list[[1]]$full_data %>% 
  mutate(facs = as.factor(paste(reg, fac, site, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(reg, fac, site, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)

data <- list(y_obs = y_obs,
             rfac = rfac,
             fx_cov = fix_mm, #model matrix from above
             n_rfac = length(unique(rfac)),
             all_fac = fac_dat$facs_n, # vector of factor combinations
             fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
             )
parameters <- list(z_rfac = rep(0, times = length(unique(rfac))),
                   z_ints = matrix(0, nrow = ncol(fix_mm),
                                   ncol = ncol(y_obs) - 1),
                   log_sigma_rfac = 0)

## Make a function object
obj <- MakeADFun(data, parameters, random = c("z_rfac"),
                 DLL = "multinomial_generic_randInt2")
# obj <- MakeADFun(data, parameters, 
#                  DLL = "multinomial_generic_fixInt")

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
sdr_small <- sdr

ssdr <- summary(sdr)
ssdr
ssdr_small <- ssdr

# make dataframe of true predictions and covariates
dumm <- dat_list[[1]]$full_data %>% 
  mutate(facs = as.factor(paste(reg, fac, site, sep = "_")))
obs_dat <- dumm %>% 
  pivot_longer(p1:p3, names_to = "cat", names_prefix = "p",
               values_to = "true_prob") %>% 
  select(cat, site, reg, fac, facs, true_prob) %>%
  distinct() %>% 
  left_join(., fac_key %>% select(facs, facs_n), by = "facs") %>% 
  arrange(facs_n) %>% 
  rbind(., .) #stack to match below

# generate predictions with and without fixed effects, then combine dataframes
logit_probs_mat <- ssdr[rownames(ssdr) %in% "logit_probs_out", ]
pred_ci_pool <- data.frame(cat = as.character(rep(1:(k + 1), 
                                                each = length(unique(fac_key$facs_n)))), 
                         logit_prob_est = logit_probs_mat[ , "Estimate"],
                         logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>% 
  mutate(ests = "pool",
         facs_n = rep(fac_key$facs_n, times = k + 1))
logit_probs_mat_fe <- ssdr[rownames(ssdr) %in% "logit_probs_out_fe", ]
pred_ci <- data.frame(cat = as.character(rep(1:(k + 1), 
                                             each = length(unique(fac_key$facs_n)))), 
                      logit_prob_est = logit_probs_mat_fe[ , "Estimate"],
                      logit_prob_se =  logit_probs_mat_fe[ , "Std. Error"]) %>%
  mutate(ests = "fix",
         facs_n = rep(fac_key$facs_n, times = k + 1)) %>% 
  rbind(pred_ci_pool, .) %>% 
  mutate(pred_prob = plogis(logit_prob_est),
         pred_prob_low = plogis(logit_prob_est +
                                  (qnorm(0.025) * logit_prob_se)),
         pred_prob_up = plogis(logit_prob_est +
                                 (qnorm(0.975) * logit_prob_se))) %>%
  left_join(., obs_dat, by = c("cat", "facs_n")) %>% 
  select(-logit_prob_est, -logit_prob_se) 
  
ggplot(pred_ci %>% filter(fac == "1", reg == "1")) +
  geom_boxplot(aes(x = as.factor(cat), y = true_prob)) +
  geom_pointrange(aes(x = as.factor(cat), y = pred_prob, ymin = pred_prob_low,
                      ymax = pred_prob_up), col = "red") +
  facet_wrap(ests ~ site, nrow = 2) +
  ggsidekick::theme_sleek()


ests <- map(dat_list, function(x) {
  y_obs <- x$obs
  rfac <- as.numeric(x$rand_fac) - 1 #subtract for indexing by 0
  n_rfac <- length(unique(rfac))
  fac_dat <- x$full_data %>% 
    mutate(facs = as.factor(paste(reg, fac, site, sep = "_")),
           facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
    select(reg, fac, site, facs, facs_n)
  fac_key <- fac_dat %>% 
    distinct() %>% 
    arrange(facs_n)
  
  data <- list(y_obs = y_obs,
               rfac = rfac,
               fx_cov = fix_mm, #model matrix from above
               n_rfac = n_rfac,
               all_fac = fac_dat$facs_n, # vector of factor combinations
               fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
  )
  parameters <- list(z_rfac = rep(0, times = length(unique(rfac))),
                     z_ints = matrix(0, nrow = ncol(fix_mm),
                                     ncol = k),
                     log_sigma_rfac = 0)
  ## Make a function object
  obj <- MakeADFun(data, parameters, random = c("z_rfac"),
                   DLL = "multinomial_generic_randInt2")
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr) 
  
  fix_eff <- data.frame(
    par = "beta",
    par_n = seq(1, length(fix_ints), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(fix_ints)
  )
  rand_eff <- data.frame(
    par = c(rep("alpha", n_rfac), "sigma_a"),
    par_n = c(seq(1, n_rfac, by = 1), 1),
    est = c(ssdr[rownames(ssdr) %in% "z_rfac" , 1], 
            ssdr[rownames(ssdr) %in% "sigma_rfac" , 1]),
    true = c(unique(x$full_data$site_mean), sd_site)
  )
  
  rbind(fix_eff, rand_eff) %>% 
    mutate(trial = x$trial) 
}) %>% 
  bind_rows()

coef_dat_hier_m <- ests %>%
  bind_rows() %>% 
  mutate(
    sq_err = (true - est)^2
  ) %>% 
  group_by(par, par_n) %>% 
  mutate(rmse = sqrt(mean(sq_err))) %>% 
  ungroup()

par_est_box <- ggplot(coef_dat_hier_m) +
  geom_boxplot(aes(x = as.factor(par_n), y = est)) +
  geom_point(aes(x = as.factor(par_n), y = true), color = "red", shape = 21) +
  ggsidekick::theme_sleek() +
  facet_wrap(~par, scales = "free")

par_rmse_box <- ggplot(coef_dat_hier_m) +
  geom_boxplot(aes(x = as.factor(par_n), y = rmse)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~par, scales = "free")

pdf(here::here("figs", "hier_multi_performance.pdf")) 
par_est_box
par_rmse_box
dev.off()


# ests %>% 
#   filter(par %in% c("g1_int", "g2_int", "reg2_g1", "reg2_g2", "reg3_g1", 
#                     "reg3_g2", "fac2_g1", "fac2_g2", "sigma_fac1")) %>%
#   ggplot(.) +
#   geom_boxplot(aes(x = par, y = Estimate)) +
#   geom_point(aes(x = par, y = vals), colour = "red")
# 
# r_effs <- ests %>% 
#   filter(grepl("z_rfac", par)) %>%
#   mutate(p_err = (Estimate - vals) / vals) 
# 
# r_effs %>% 
#   select(par, est = Estimate, true = vals) %>% 
#   pivot_longer(-par, names_to = "type", values_to = "value") %>% 
#   ggplot(.) +
#   geom_boxplot(aes(x = type, y = value)) +
#   facet_wrap(~par)
# 
# ggplot(r_effs) +
#   geom_boxplot(aes(x = par, y = p_err))


# MISC -------------------------------------------------------------------------

# # experiment w/ increasing the size of dataset while keeping proportions the same
# dd <- vector(mode = "list", length = 20)
# y_obs <- NULL
# reg_dat_long <- NULL
# rfac_long <- NULL
# fac_dat_long <- NULL
# fac_key_long <- NULL
# for(i in seq_along(dd)) {
#   y_obs <- rbind(y_obs, dat_list[[1]]$obs)
#   reg_dat_long <- rbind(reg_dat_long, reg_dat)
#   
#   rfac_long <- c(rfac_long, as.numeric(dat_list[[1]]$rand_fac) - 1) 
#   fac_dat_long <- rbind(fac_dat_long, fac_dat)
# }
# 
# fix_mm <- model.matrix(~ (reg + fac), reg_dat_long)
# data <- list(y_obs = y_obs,
#              rfac = rfac_long,
#              fx_cov = fix_mm, #model matrix from above
#              n_rfac = length(unique(rfac_long)),
#              all_fac = fac_dat_long$facs_n, # vector of factor combinations
#              fac_key = fac_key$facs_n #ordered unique factor combos in fac_vec
# )
# unsurprisingly reduces SE substantially