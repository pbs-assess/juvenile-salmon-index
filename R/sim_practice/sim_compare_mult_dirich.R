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
  
  ### TRY NEW DISTRIBUTION
  # set.seed(123)
  # dirmult::rdirichlet(1, pi[1, ] * (1 - theta[1]) / theta[1])
  
  #number of values to simulate per row of input
  # N <- sample(c(1:1000), n, replace = T)
  N <- rep(100, n)
  
  # Simulate the responses Y from package 'dirmult' and make 4 versions for 
  # different model formats
  # set.seed(123 + trial)
  Y = simPop(J = 1, n = N[1], pi = pi_prob[1,], theta = theta[1])$data
  for(jj in c(2:n)){
    Y = rbind(Y, simPop(J = 1, n = N[jj], pi = pi_prob[jj,], 
                        theta = theta[jj])$data)
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
n_trials <- 100
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

# helper function to generate comparison of coefficients to true values 
make_coef_df <- function(betas, est_betas, tran, trial) {
  data.frame(
    beta = seq(1, length(betas), by = 1),
    est = est_betas,
    true = as.vector(betas),
    tran = tran,
    trial = trial)
} 


fit_list <- map(sim_list, function(sims_in) {
  ## Fit dirichlet TMB wtih raw data (i.e. integers for each category)
  Y_d1 <- sims_in$raw$obs
  obj <- MakeADFun(data=list(fx_cov = fix_mm,
                             y_obs = Y_d1,
                             pred_cov = pred_mm_d),
                   parameters=list(z_ints = beta_in_d),
                   DLL="dirichlet_fixInt")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr_dir <- sdreport(obj)
  ssdr_dir <- summary(sdr_dir)
  coef_dat_dir <- make_coef_df(betas, 
                               est_betas = ssdr_dir[rownames(ssdr_dir) %in% "z_ints", 1], 
                               tran = sims_in$raw$trans,
                               trial = sims_in$raw$trial)
  
  ## Fit dirichlet TMB with adjusted data (i.e. proportions multiplied by scalar)
  Y_d2 <- sims_in$adj$obs
  obj <- MakeADFun(data=list(fx_cov = fix_mm,
                             y_obs = Y_d2,
                             pred_cov = pred_mm_d),
                   parameters=list(z_ints = beta_in_d),
                   DLL="dirichlet_fixInt")
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  sdr_dir2 <- sdreport(obj)
  ssdr_dir2 <- summary(sdr_dir2)
  coef_dat_dir2 <- make_coef_df(betas, 
                                est_betas = ssdr_dir2[rownames(ssdr_dir2) %in% "z_ints", 1], 
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
  # comp_dat <- sims_in$raw$obs / rowSums(sims_in$raw$obs) 
  temp$comp <- DR_data(sims_in$adj$obs / 100)

  fit_dr <- DirichletReg::DirichReg(comp ~ reg + fac, data = temp)
  betas_dr <- coef(fit_dr) %>%
    do.call(c, .) %>%
    as.vector()
  coef_dat_DR <- make_coef_df(betas,
                              est_betas = betas_dr,
                              tran = "DR",
                              trial = sims_in$adj$trial)
  
  # output lists
  obs_list <- list(dir1 = Y_d1, 
                   dir2 = Y_d2, 
                   multi = Y_m, 
                   DR = (Y_d2 / 100)
                   )
  fits_list <- list(dir1 = ssdr_dir, 
                    dir2 = ssdr_dir2, 
                    multi = ssdr_m, 
                    DR = fit_dr
                    )
  coef_list <- list(dir1 = coef_dat_dir, 
                    dir2 = coef_dat_dir2, 
                    DR = coef_dat_DR
                    )
  return(list(obs = obs_list, fits = fits_list, coefs = coef_list))
})
saveRDS(fit_list, here::here("data", "modelFits", "multi_dir_comp.RDS"))


# Compare coefficient estimates ------------------------------------------------
## excludes multinomial because estimated parameters are different
fit_list <- readRDS(here::here("data", "modelFits", "multi_dir_comp.RDS"))

coef_dat <- map(fit_list, function (x) {
  x$coefs %>% 
    bind_rows()
}) %>% 
  bind_rows()

ggplot(coef_dat) + 
  geom_boxplot(aes(x = as.factor(tran), y = est)) +
  geom_hline(aes(yintercept = true)) +
  facet_wrap(~beta, scales = "free_y")


# Compare predictions ----------------------------------------------------------

# function to generate dataframe of true vales for each model type based on
# different input structures
make_true_df <- function(model, trial, y_obs) {
  colnames(y_obs) <- seq(1, K, by = 1)
  #calculate observed row sum to convert mean observations to proportions
  row_sum <- apply(y_obs, 1, sum)
  reg_dat %>% 
    #merge common covariate data to unique observation data
    cbind(., y_obs, row_sum) %>% 
    pivot_longer(., cols = 3:(3 - 1 + K), names_to = "cat", 
                 values_to = "prob") %>%  
    mutate(#model = model,
           trial = trial) %>% 
    group_by(reg, fac, cat) %>%
    mutate(
      true_prob = mean(prob / row_sum)
      # true_prob = case_when(
      #   model %in% c("dir_adj", "dir_raw") ~ mean(prob / row_sum),
      #   model == "multi" ~ sum(prob) / sum(row_sum)
      #   )
      ) %>%
    select(trial, reg, fac, cat, true_prob) %>% 
    distinct()
}


# function to generate dataframe of predictions relative to true values
make_pred_df <- function(model, true_dat, ssdr_in) {
  dum <- expand.grid(reg = pred_dat_d$reg,
              fac = pred_dat_d$fac,
              cat = as.character(seq(1, K, by = 1))) %>% 
    distinct() %>%
    arrange(cat, reg) %>%
    left_join(., true_dat, by = c("reg", "fac", "cat")) %>% 
    mutate(model = model)
  
  if (model == "multi") {
    out <- dum %>% 
      mutate(
        logit_prob_est = ssdr_in[ , "Estimate"],
        logit_prob_se = ssdr_in[ , "Std. Error"],
        mu = plogis(logit_prob_est),
        lo = plogis(logit_prob_est + (qnorm(0.025) * logit_prob_se)),
        up = plogis(logit_prob_est + (qnorm(0.975) * logit_prob_se)),
        covered = lo < true_prob & up > true_prob
      ) %>% 
      select(-logit_prob_est, -logit_prob_se)
  } 
  
  if (model %in% c("dir_adj", "dir_raw")) {
    out <- dum %>% 
      mutate(mu = ssdr_in[ , "Estimate"],
             se = ssdr_in[ , "Std. Error"],
             lo = mu + (qnorm(0.025) * se),
             up = mu + (qnorm(0.975) * se),
             covered = lo < true_prob & up > true_prob) %>% 
      select(-se)
  }
  return(out)
}


# generate predictions for different elements of fit_list (excludes DR for now
# due to different structure)
pred_list <- map(fit_list, function (x) {
  # pull trial from coef dataframe 
  trial_n <- unique(x$coefs$dir1$trial)
  true_df_in <- make_true_df(model = "dir_raw", trial = trial_n,
                             y_obs = x$obs$dir1)
  
  preds_dir1 <- x$fits$dir1[rownames(x$fits$dir1) %in% "pred_pi_prop" , ]
  pred_df_dir1 <- make_pred_df(model = "dir_raw", true_dat = true_df_in,
                               ssdr_in = preds_dir1)
  
  preds_dir2 <- x$fits$dir2[rownames(x$fits$dir2) %in% "pred_pi_prop" , ]
  pred_df_dir2 <- make_pred_df(model = "dir_adj", true_dat = true_df_in,
                               ssdr_in = preds_dir2)
  
  preds_m <- x$fits$multi[rownames(x$fits$multi) %in% "logit_probs_out" , ]
  pred_df_m <- make_pred_df(model = "multi", true_dat = true_df_in,
                            ssdr_in = preds_m)
  
  rbind(pred_df_dir1, pred_df_dir2, pred_df_m)
})

pred_out <- pred_list %>% 
  bind_rows() 
pred_out$sq_err <- (pred_out$mu - pred_out$true_prob)^2

# calculate rmse and average coverage of predicted proportions per trial for 
# each model
fig_summ <- pred_out %>% 
  group_by(model, trial) %>% 
  summarize(coverage = mean(covered),
            rmse = mean(sq_err))  %>% 
  pivot_longer(., cols = c("coverage", "rmse"), names_to = "perf_ind", 
               values_to = "est")

ggplot(fig_summ) +
  geom_boxplot(aes(x = model, y = est)) +
  facet_wrap(~perf_ind, scales = "free_y")
