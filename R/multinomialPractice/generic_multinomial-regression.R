## Simulate and fit unordered multinomial response models with k categories
# Based on multinomial-regression.R

library(tidyverse)

# structure inputs
set.seed(42)
N <- 1000 # number of observations
k <- 3 #number of groups
# X <- runif(N) #cont. predictor
yrCov <- data.frame(yr = sample(c("yr1", "yr2"), size = N, replace = T)) #fixed predictor representing two cats
.X <- model.matrix(~yr, yrCov)

int <- c(0.3, -1.4)
beta <- c(0.5, -0.5)

#combine intercepts and betas into one matrix that can be multiplied by model
#matrix; columns are parameters for each group (excluding ref. cat.)
betas <- matrix(c(int, beta), nrow = ncol(.X), ncol = k - 1, byrow = T)
m <- nrow(betas) #number of covariates (i.e. 1 int. + slopes)

log_odds <- matrix(NA, nrow = N, ncol = k)
for (h in 1:k) {
  if (h < k) {
    cov_effects <- matrix(NA, nrow = N, ncol = m)
    for (i in seq_len(N)) {
      for (j in seq_len(m)) {
        cov_effects[i, j] <- betas[j , h] * .X[i, j]  
      }
      log_odds[i, h] <- sum(cov_effects[i, ])
    }
  } else if (h == k) {
    log_odds[ , h] <- NA
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

plot(.X[ , 2], jitter(y, 0.3))

#matrix of observations
Yobs <- matrix(ncol = k, nrow = N, data = 0)
for (i in seq_along(y)) {
  Yobs[i, y[i]] <- 1
}

# saveRDS(Yobs, here::here("R", "multinomialPractice", "exDat.RDS"))
Yobs <- readRDS(here::here("R", "multinomialPractice", "exDat.RDS"))


set.seed(42)
yrCov <- data.frame(yr = sample(c("yr1", "yr2"), size = nrow(Yobs), 
                                replace = T)) #fixed predictor representing two cats
.X <- model.matrix(~yr, yrCov)

nll <- function(pars, Y, covMatrix) {
  k <- ncol(Y)
  N <- nrow(Y) # number of rows in observation matrix
  
  # define intercepts and betas based on pars vector and add reference 
  # category to vector of parameters
  betas <- matrix(pars, nrow = ncol(covMatrix), ncol = k - 1, byrow = T)
  m <- nrow(betas) #number of covariates (i.e. 1 int. + slopes)
  
  log_odds <- matrix(NA, nrow = N, ncol = k)
  for (h in 1:k) {
    if (h < k) {
      cov_effects <- matrix(NA, nrow = N, ncol = m)
      for (i in seq_len(N)) {
        for (j in seq_len(m)) {
          cov_effects[i, j] <- betas[j , h] * covMatrix[i, j]  
        }
        log_odds[i, h] <- sum(cov_effects[i, ])
      }
    } else if (h == k) {
      log_odds[ , h] <- NA
    }
  }
  exp_log_odds <- exp(log_odds)
  
  denominator <- rep(NA, length.out = N)
  for (i in 1:N) {
    denominator[i] <- 1 + sum(exp_log_odds[i, 1:(k-1)])
  }  
  
  probs <- matrix(NA, nrow = N, ncol = k)
  for (h in 1:k) {
    if (h < k) {
      probs[ , h] <- exp_log_odds[ , h] / denominator
    } 
    else if (h == k) {
      for (i in 1:N) {
        #not sure if we can vectorize this without apply...
        probs[i, h] <- 1 - sum(probs[i, 1:(k-1)]) 
      }
    }
  }
  
  #log-likelihood
  nll <- vector(length = N)
  for (i in seq_len(N)) { # not vectorized
    nll[i] <- -dmultinom(Y[i, ], size = 1, prob = probs[i, ], log = TRUE)
  }
  sum(nll)
}

set.seed(42)
# .X <- matrix(data = runif(N), nrow = N, ncol = 1) #predictor
yrCov <- data.frame(yr = sample(c("yr1", "yr2"), size = nrow(Yobs), 
                                     replace = T)) #fixed predictor representing two cats
.X <- model.matrix(~yr, yrCov)

# Fit with initial values
m2 <- nlminb(rep(0, 4), nll, Y = Yobs, covMatrix = .X)
m2

