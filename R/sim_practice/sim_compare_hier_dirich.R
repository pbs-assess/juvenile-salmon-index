# Simulate unordered dirichlet response models then fit TMB dirichlet model,
# DirichletReg model, and TMB multinomial and generate predictions to evaluate
# performance
# May 14, 2020

library(dirmult)
library(TMB)
library(tidyverse)
library(DirichletReg)


## Simulate data ---------------------------------------------------------------

set.seed(42)

# random intercepts for groups (e.g. sites)
n_sites <- 5
n_obs_per_site <- 100
sd_site <- 0.5
n <- n_sites * n_obs_per_site

#covdat
reg_dat <- data.frame(reg = as.factor(sample(c(1, 2, 3), size = n, 
                                             replace = T)),
                      fac = as.factor(sample(c(1, 2), size = n, replace = T))
)
# model matrix for fixed effects
fix_mm <- model.matrix(~ (reg + fac), reg_dat)

# fixed intercepts (P x K)
P <- ncol(fix_mm) - 1
K <- 4
betas <- matrix(rnorm((P+1)*K), (P+1), K)

# function to simulate dirichlet data based on parameters and matrix
f_sim <- function(trial = 1, scalar = 100) {
  
  #calculate fixed effects using parameter matrix and model matrix
  n_fix_cov <- P
  fix_eff <- fix_mm %*% betas
  sum_fix_eff <- apply(fix_eff, 1, sum)
  
  # generate random effects, then combine with fixed (not incorporated yet)
  site_mean_a <- rnorm(mean = 0,
                       sd = sd_site,
                       n = n_sites)
  site_obs_ia <- rnorm(mean = rep(site_mean_a, each = n_obs_per_site),
                       sd = 0,
                       n = n)
  site_dat <- data.frame(site = as.factor(rep(seq(1, 5, by = 1),
                                              each = n_obs_per_site)),
                         site_mean = rep(site_mean_a, each = n_obs_per_site),
                         site_obs = site_obs_ia,
                         sd_site = sd_site
  ) %>%
    cbind(., reg_dat, sum_fix_eff)  #bind regional fixed effects data
  
  #calc gamma
  Gamma = exp(fix_eff) #fixed effects
  Gamma_plus = apply(Gamma, 1, sum) #sum of fixed_effects
  theta = 1 / (Gamma_plus + 1)
  pi = apply(Gamma, 2, function(x) {x / theta})
  pi_prob = pi / rowSums(pi)
  
  #number of values to simulate per row of input
  # N <- sample(c(1:1000), n, replace = T)
  N <- rep(100, n)
  
  # Simulate the responses Y from package 'dirmult' and make 4 versions for 
  # different model formats
  # set.seed(123 + trial)
  Y = simPop(J = 1, n = N[1], pi = pi_prob[1,], theta = theta[1])$data
  for(jj in c(2:n)){
    Y = rbind(Y, simPop(J = 1, n = N[jj], pi = pi_prob[jj,], theta = theta[jj])$data)
  }
  
  #switch to proportions and add small value
  Y_ppn <- Y / rowSums(Y)  
  Y_ppn <- Compositional::zeroreplace(Y_ppn, a = 2/3)  
  #scalar is necessary to produce realistic SEs
  Y_adj <- Y_ppn * scalar
  
  # switch to binary data based on max value
  maxs <- apply(Y_adj, 1, max)
  Y_int <- ifelse(Y_adj == maxs, 1, 0)
    
  return(list("trial" = trial, "obs_raw" = Y, "obs_adj" = Y_adj, 
              "obs_int" = Y_int
  ))
}

# generate data
n_trials <- 5
sim_list <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  sims_dat <- f_sim(trial = i, scalar = 100)
  sims_raw <- list(trial = i, obs = sims_dat$obs_raw, trans = "raw")
  sims_adj <- list(trial = i, obs = sims_dat$obs_adj, trans = "adj")
  sims_int <- list(trial = i, obs = sims_dat$obs_int, trans = "int")
  sim_list[[i]] <- list(raw = sims_raw, adj = sims_adj, int = sims_int)
}

## Set up model inputs ---------------------------------------------------------

#Dirichlet TMB inputs
#initial parameter values
beta_in_d <- matrix(rnorm(n = (P + 1)*K, mean = 0), (P+1), K)

#prediction model matrix
pred_dat_d <- reg_dat %>%
  select(reg, fac) %>%
  distinct() %>%
  arrange(reg, fac)
pred_mm_d <- model.matrix(~ reg + fac, pred_dat_d)

compile(here::here("src", "dirichlet_fixInt.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_fixInt")))


# Multinomial TMB inputs
#initial parameter values
beta_in_m = matrix(rnorm(n = (P + 1) * (K - 1), mean = 0), nrow = P + 1,
                   ncol = K - 1)

fac_dat <- reg_dat %>% 
  mutate(facs = as.factor(paste(reg, fac, sep = "_")),
         facs_n = (as.numeric(facs) - 1)) %>% #subtract for indexing by 0 
  select(reg, fac, facs, facs_n)
fac_key <- fac_dat %>% 
  distinct() %>% 
  arrange(facs_n)

compile("src/multinomial_fixInt.cpp")
dyn.load(dynlib("src/multinomial_fixInt"))


# Fit models -------------------------------------------------------------------
make_coef_df <- function(betas, ssdr, tran, trial) {
  data.frame(
    beta = seq(1, length(betas), by = 1),
    est = ssdr[rownames(ssdr) %in% "z_ints" , 1],
    true = as.vector(betas),
    tran = tran,
    trial = trial)
} 


fit_list <- map(sim_list, function(sims_in) {
  
  sims_in <- sim_list[[1]]
  
  ## Fit dirichlet TMB wtih raw data
  Y_d1 <- sims_in$raw$obs
  obj <- MakeADFun(data=list(fx_cov = fix_mm,
                             y_obs = Y_d1,
                             pred_cov = pred_mm_d),
                   parameters=list(z_ints = beta_in_d),
                   DLL="dirichlet_fixInt")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr_dir <- sdreport(obj)
  ssdr_dir <- summary(sdr_dir)
  coef_dat_dir <- make_coef_df(betas, ssdr = ssdr_dir, tran = sims_in$raw$trans,
                               trial = sims_in$raw$trial)
  
  ## Fit dirichlet TMB with adjusted data (proportions multiplied by scalar)
  Y_d2 <- sims_in$adj$obs
  obj <- MakeADFun(data=list(fx_cov = fix_mm,
                             y_obs = Y_d2,
                             pred_cov = pred_mm_d),
                   parameters=list(z_ints = beta_in_d),
                   DLL="dirichlet_fixInt")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr_dir2 <- sdreport(obj)
  ssdr_dir2 <- summary(sdr_dir2)
  coef_dat_dir2 <- make_coef_df(betas, ssdr = ssdr_dir2, 
                                tran = sims_in$adj$trans,
                                trial = sims_in$adj$trial)
  
  ## Fit multinomial TMB with integer data 
  Y_m <- sims_in$int$obs
  obj <- MakeADFun(data = list(fx_cov = fix_mm,
                               y_obs = Y_m,
                               all_fac = fac_dat$facs_n,
                               fac_key = fac_key$facs_n),
                   parameters = list(z_ints = beta_in_m),
                   DLL="multinomial_fixInt")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr_m <- sdreport(obj)
  ssdr_m <- summary(sdr_m)
  
  ## Fit dirichlet regression from DirichletReg
  temp <- reg_dat
  temp$comp <- DR_data(sims_in$adj$obs / 100)
  
  fit_dr <- DirichletReg::DirichReg(comp ~ reg, data = temp)
  # betas_out <- coef(m1) %>%
  #   do.call(c, .) %>%
  #   as.vector()
  
  temp <- list(ssdr_dir, ssdr_dir2, ssdr_m, fit_dr, coef_dat_dir, coef_dat_dir2)
})

# function to generate dataframe of predictions for each model type
make_pred_df <- function(model, y_obs, ppn_mat) {
  # create true data values
  colnames(y_obs) <- seq(1, K, by = 1)
  row_sum <- apply(y_obs, 1, sum)
  true_df <- reg_dat %>% 
    cbind(., y_obs, row_sum) %>% 
    pivot_longer(., cols = 3:(3 - 1 + K), names_to = "cat", 
                 values_to = "prob") %>%  
    mutate(model = model) %>% 
    ungroup() %>% 
    group_by(reg, fac, cat) %>%
    mutate(
      true_prob = case_when(
        model %in% c("dir_adj", "dir_raw") ~ mean(prob) / row_sum
        )
      ) %>%
    select(reg, fac, cat, true_val) %>% 
    distinct()
  pred_dat <- expand.grid()
}



pred_prob = plogis(logit_prob_est),
pred_prob_low = plogis(logit_prob_est +
                         (qnorm(0.025) * logit_prob_se)),
pred_prob_up = plogis(logit_prob_est +
                        (qnorm(0.975) * logit_prob_se))

# create data frame of predictions
ppn <- ssdr[rownames(ssdr) %in% "pred_pi_prop" , ]
ppn <- ssdr_m[rownames(ssdr_m) %in% "logit_probs_out" , ]
pred_dat <- expand.grid(reg = pred_dat_d$reg,
                        fac = pred_dat_d$fac,
                        cat = as.character(seq(1, K, by = 1))) %>% 
  distinct() %>%
  arrange(cat, reg) %>%
  left_join(., true_values, by = c("reg", "fac", "cat")) %>% 
  mutate(mu = ppn[ , "Estimate"],
         se = ppn[ , "Std. Error"],
         lo = mu + (qnorm(0.025) * se),
         up = mu + (qnorm(0.975) * se),
         covered = lo < true & up > true) 

true_df %>% 
  group_by(reg, fac) %>% 
  summarize(sum(true_val))



true_values_dir <- reg_dat %>% 
  cbind(., Y_d1) %>%
  # glimpse()
  pivot_longer(., cols = `1`:`4`, names_to = "cat", values_to = "prob") %>% 
  group_by(reg, fac, cat) %>%
  summarize(true = mean(prob) / 100)


if (model == "dir_raw") {
  y_obs <- data_list$raw$obs
} else if (model == "dir_adj") {
  y_obs <- data_list$adj$obs
} else if (model == "multi") {
  y_obs <- data_list$int$obs
} else if (model == "DR") {
  y_obs <- data_list$adj$obs / 100
}
