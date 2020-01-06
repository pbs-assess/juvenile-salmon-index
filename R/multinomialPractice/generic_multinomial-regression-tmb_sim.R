## Simulation test of model code in generic_multinomial-regression-tmb.R

# Model function ---------------------------------------------------------------

set.seed(1)

# Compile model
library(TMB)
compile("R/multinomialPractice/multinomial_generic.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic"))

#random fields simulation function
.rf_sim <- function(model, x, y) {
  out <- sdmTMB:::rf_sim(model, x, y)
  out - mean(out)
}

rf_pars <- list(sig = 0.3, kappa = 1.7)

fitModel <- function(X, N, k, ints, betas, rf_pars = NULL) {
  log_odds <- matrix(NA, nrow = N, ncol = (k - 1))
  sp_list <- vector(mode = "list", length = k)
  for (h in 1:(k - 1)) {
    # normal calculation of log_odds w/out random fields
    if (is.null(rf_pars)) {
      log_odds[ , h] <- ints[h] + betas[h] * X
    }
    
    # with simulated draws from random field 
    if (!is.null(rf_pars)) {
      rf_model <- RandomFields::RMmatern(nu = 1, var = rf_pars$sig^2, 
                                         scale = 1/rf_pars$kappa)
      #assumes relatively coarse grid size currently
      rf_dat <- expand.grid(x = seq(0, 1, length.out = 20), 
                            y = seq(0, 1, length.out = 20)) %>% 
        mutate(rf_effect = .rf_sim(model = rf_model, .$x, .$y)) %>% 
        sample_n(., size = N)

      #### TEMP LIST TO LOOK AT SPATIAL OBSERVATIONS ####
      sp_list[[h]] <- rf_dat      
      
      log_odds[ , h] <- ints[h] + betas[h] * X + rf_dat$rf_effect
    }
  }
  exp_log_odds <- exp(log_odds)

  denominator <- rep(NA, length.out = N)
  for (i in 1:N) {
    denominator[i] <- 1 + sum(exp_log_odds[i, ])
  }
  
  probs <- matrix(NA, nrow = N, ncol = k)
  for (h in 1:k) {
    if (h < k) {
      probs[ , h] <- exp_log_odds[ , h] / denominator
    } 
    else if (h == k) {
      for (i in 1:N) {
        probs[i, h] <- 1 - sum(probs[i, 1:(k-1)]) 
      }
    }
  }
  
  #generate observations
  y <- numeric(length = N)
  for (i in seq_along(y)) {
    temp <- rmultinom(1, 1, probs[i, ])
    y[i] <- which(temp == 1)
  }

  #matrix of observations
  y_obs <- matrix(ncol = k, nrow = N, data = 0)
  for (i in seq_along(y)) {
    y_obs[i, y[i]] <- 1
  }
  
  ## temporary vis
  tt <- cbind(y_obs[ , 1], sp_list[[1]][, 1:2])
  
  ggplot(sp_list[[1]], aes(x, y, fill = rf_effect)) +
    geom_raster() +
    scale_fill_gradient2() 
  ggplot(tt, aes(x, y, fill = y_obs[, 1])) +
    geom_raster() +
    scale_fill_gradient2() 
  
    
  ## Data and parameters
  .X <- cbind(1, X) #predictor with intercept
  data <- list(cov=.X, y_obs = y_obs)
  parameters <- list(betas = matrix(data = 0, nrow = ncol(.X), 
                                    ncol = ncol(y_obs) - 1))
  
  ## Make a function object
  obj <- MakeADFun(data, parameters, DLL="multinomial_generic", silent = TRUE)
  
  ## Call function minimizer
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  ## Get parameter uncertainties and convergence diagnostics
  sdr <- sdreport(obj)
  ssdr <- summary(sdr)
  # r <- obj$report()
  outlist <- list(probs,
                  sdreport(obj), 
                  summary(sdr), 
                  obj$report())
  names(outlist) <- c("true_probs", "sdr", "ssdr", "r")
  return(outlist)
}

# Run Model --------------------------------------------------------------------

library(tidyverse)

N <- 300 # number of observations
k <- 4 #number of groups
ints <- c(0.3, -1.4, 0.3) #intercepts
betas <- c(-3, 4, -1) #slopes
trials <- 2 #number of times to run function

# keep predictor as separate list so it can be used when making DF
pred_list <- lapply(seq_len(trials), function(k) {
  runif(N)
})
names(pred_list) <- seq_len(trials)

sim_list <- map(pred_list, fitModel, N = N, k = 4, ints = ints, 
                                       betas = betas)

# Collapse model runs into dataframe 
coef_list <- lapply(seq_along(sim_list), function(x) {
  dum <- sim_list[[x]]
  
  # estimated vs. true coefficients
  coef_mat <- dum$ssdr[rownames(dum$ssdr) %in% "betas", ]
  coefs <- data.frame(cat = as.factor(rep(1:(k - 1), each = 2)),
                      trial = seq(1, trials, by = 1)[x],
                      pars = rep(c("int", "slope"), times = (k - 1)),
                      est = coef_mat[ , "Estimate"],
                      se = coef_mat[ , "Std. Error"]) %>% 
    arrange(pars) %>% 
    mutate(low = est + (qnorm(0.025) * se),
           up = est + (qnorm(0.975) * se),
           true_val = c(ints, betas))
  
  # predicted vs. true probabilities
  pred <- pred_list[[x]]
  logit_probs_mat <- dum$ssdr[rownames(dum$ssdr) %in% "logit_probs", ] 
  pred_ci <- data.frame(cat = rep(1:k, each = N),
                        trial = seq(1, trials, by = 1)[x],
                        X = rep(pred, times = 4),
                        prob = as.vector(dum$true_probs),
                        logit_prob_est = logit_probs_mat[ , "Estimate"],
                        logit_prob_se =  logit_probs_mat[ , "Std. Error"]) %>% 
    mutate(pred_prob = plogis(logit_prob_est),
           pred_prob_low = plogis(logit_prob_est + (qnorm(0.025) * logit_prob_se)),
           pred_prob_up = plogis(logit_prob_est + (qnorm(0.975) * logit_prob_se)))
  
  #pred list
  fit_list <- list(coefs, pred_ci)
  names(fit_list) <- c("coefs", "pred_ci")
  return(fit_list)
})

coef_dat <- map(coef_list, "coefs") %>% bind_rows()
pred_ci_dat <- map(coef_list, "pred_ci") %>% bind_rows()

# Plot
ggplot(coef_dat %>%
         filter(pars == "int"), 
       aes(x = cat, y = est)) +
  geom_pointrange(aes(ymin = low, ymax = up)) +
  geom_point(aes(x = cat, y = true_val), shape = 3, color = "red") +
  labs(y = "Intercept", x = "Category") +
  facet_wrap(~trial)
ggplot(coef_dat %>%
         filter(pars == "slope"), 
       aes(x = cat, y = est)) +
  geom_pointrange(aes(ymin = low, ymax = up)) +
  geom_point(aes(x = cat, y = true_val), shape = 3, color = "red") +
  labs(y = "Slope", x = "Category") +
  facet_wrap(~trial)

ggplot(pred_ci_dat) +
  geom_ribbon(aes(x = X, ymin = pred_prob_low, ymax = pred_prob_up), 
              fill = "#bfd3e6") +
  geom_line(aes(x = X, y = pred_prob), col = "#810f7c", size = 1) +
  geom_line(aes(x = X, y = prob), size = 1) +
  facet_grid(cat~trial) +
  labs(y = "Probability", x = "Random Predictor")
  
# Coverage calculations --------------------------------------------------------

coef_dat %>%
  group_by(pars) %>%
  mutate(covered = low < true_val & up > true_val) %>%
  summarize(coverage = mean(covered))

coef_dat %>%
  group_by(pars) %>%
  mutate(relative_error = (est - true_val) / true_val) %>%
  summarize(mean_relative_error = mean(relative_error))
