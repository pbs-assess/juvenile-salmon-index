# TMB ------------------------------------------------------

set.seed(1)
N <- 2000 # number of observations
k <- 4 #number of groups
X <- runif(N) #predictor
ints <- c(0.3, -1.4, 0.3)
betas <- c(-3, 4, -1)

# append with reference category
int <- c(ints, NaN) 
beta <- c(betas, NaN)

log_odds <- matrix(NA, nrow = N, ncol = k)
for (i in 1:k) {
  if (i < k) {
    log_odds[ , i] <- int[i] + beta[i] * X
  } else if (i == k) {
    log_odds[ , i] <- NA
  }
}
exp_log_odds <- exp(log_odds)

denominator <- rep(NA, length.out = N)
for (i in 1:N) {
  denominator[i] <- 1 + sum(exp_log_odds[i, 1:(k-1)])
}  

probs <- matrix(NA, nrow = N, ncol = k)
for (i in 1:k) {
  if (i < k) {
    probs[ , i] <- exp_log_odds[ , i] / denominator
  } 
  else if (i == k) {
    for (h in 1:N) {
      probs[h, i] <- 1 - sum(probs[h, 1:(k-1)]) 
    }
  }
}

#generate observations
y <- numeric(length = N)
for (i in seq_along(y)) {
  temp <- rmultinom(1, 1, probs[i, ])
  y[i] <- which(temp == 1)
}
plot(X, jitter(y, 0.3))
head(cbind(x = X, y))

#matrix of observations
y_obs <- matrix(ncol = k, nrow = N, data = 0)
for (i in seq_along(y)) {
  y_obs[i, y[i]] <- 1
}

library(TMB)
compile("R/multinomialPractice/multinomial_generic.cpp")
dyn.load(dynlib("R/multinomialPractice/multinomial_generic"))

## Data and parameters
.X <- cbind(1, X) #predictor with intercept
data <- list(cov=.X, y_obs = y_obs)
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


library(tidyverse)

# Visualize
trim_ssdr <- ssdr[rownames(ssdr) %in% "logit_probs", ] 
pred_ci <- data.frame(cat = rep(1:k, each = N),
                      est = trim_ssdr$Estimate,
                      se =  trim_ssdr$`Std. Error`) %>% 
  mutate(pred_est = plogis(est),
         pred_low = plogis(est - (2 * se)),
         pred_up = plogis(est + (2 * se)))

pred_ci %>% 
  group_by(cat) %>% 
  summarize(meanProb = mean(pred_est))
