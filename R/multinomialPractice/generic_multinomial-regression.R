## Simulate and fit unordered multinomial response models with k categories
# Based on multinomial-regression.R

library(tidyverse)

# structure inputs
set.seed(42)
N <- 100 # number of observations
k <- 3 #number of groups
X <- runif(N) #predictor
ints <- c(0.3, -1.4)
betas <- c(-3, 4)

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
      #not sure if we can vectorize this without apply...
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
Yobs <- matrix(ncol = k, nrow = N, data = 0)
for (i in seq_along(y)) {
  Yobs[i, y[i]] <- 1
}

saveRDS(Yobs, here::here("R", "multinomialPractice", "exDat.RDS"))
Yobs <- readRDS(here::here("R", "multinomialPractice", "exDat.RDS"))

### Generic function 2 - estimate probabilities from input matrix (i.e. main
# function)
#Gen function 2 inputs 
# swap out ints and betas for pure vector if necessary
pars <- c(ints, betas)  
pars <- c(0.3, -1.4, -3, 4, 1, 1, 2, 2)

set.seed(42)
N <- nrow(Yobs) # number of observations
k <- ncol(Yobs) # number of groups
covMatrix <- matrix(data = runif(N), nrow = N, ncol = 1) #predictor


nll <- function(pars, Y, k, covMatrix) {
  N <- nrow(Y) # number of rows in observation matrix
  
  # define intercepts and betas based on pars vector and add  reference 
  # category to vector of parameters
  numBetas <- ncol(covMatrix)
  int <- c(pars[1:(k-1)], NA)
  betasVec <- pars[k:length(pars)]
  beta <- matrix(NA, nrow = numBetas, ncol = k) 
  for (j in 1:numBetas) {
    beta[j , ] <- c(betasVec[(j * 2 - 1):(j * 2)], NA)
  }
  ### adjusted to allow for multiple covariate effects, but not very elegant...
  # calculate log-odds for each category 
  log_odds <- matrix(NA, nrow = N, ncol = k)
  for (h in 1:k) {
    if (h < k) {
      #calculate cumulative covariate effects
      covEffects <- matrix(NA, nrow = N, ncol = numBetas)
      for (j in 1:numBetas) {
        covEffects[ , j] <- beta[j, h] * covMatrix[ , j]
      }
      for (i in 1:N) {
        log_odds[i , h] <- int[h] + sum(covEffects[i, ])
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
X <- matrix(data = runif(N), nrow = N, ncol = 1) #predictor

# Fit with initial values
m2 <- nlminb(rep(0, 4), nll, Y = Yobs, k = 3, covMatrix = X)
m2

