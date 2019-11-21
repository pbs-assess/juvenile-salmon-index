## Practice fitting TMB version of 3-category, then generic multinomial models


Yobs <- readRDS(here::here("R", "multinomialPractice", "exDat.RDS"))
N <- nrow(Yobs) # number of observations
k <- ncol(Yobs) # number of groups
covMatrix <- matrix(data = runif(N), nrow = N, ncol = 1) #predictor

# Data 
data <- list()
data$Yobs <- Yobs
data$cov <- covMatrix
data$k <- k

# Define parameters
nBetas <- ncol(covMatrix)
int <- rep(0, times = k - 1)
betas <- matrix(0, nrow = nBetas, ncol = k)

parameters = list(
  ints,
  betas
)




library(Rcpp)
pars <- c(0.3, -1.4, -3, 4, 1, 1, 2, 2)
ints <- c(0.3, -1.4)
betas <- matrix(c(-3, 4), nrow = 1, ncol = 2)

logOdds <- function(ints, betas, covMatrix, N, k) {
  log_odds <- matrix(NA, nrow = N, ncol = (k-1))
  numBetas <- ncol(covMatrix)
  for (h in 1:(k-1)) {
    #calculate cumulative covariate effects
    covEffects <- matrix(NA, nrow = N, ncol = numBetas)
    for (j in 1:numBetas) {
      covEffects[ , j] <- betas[j, h] * covMatrix[ , j]
    }
    for (i in 1:N) {
      log_odds[i , h] <- ints[h] + sum(covEffects[i, ])
    }
  }
  return(log_odds)
}

logOdds(ints, betas, covMatrix, N, k)

cppFunction('NumericMatrix logOddsC(NumericVector ints, NumericMatrix betas,
int k) {
  NumericMatrix mat(k, k);
  
  for(int i = 0, )
  
  return(mat);
}')

cppFunction('NumericMatrix matrixC(int k) {
            NumericMatrix mat(k, k);
            return(mat);
}')









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
nll <- vector(length = length(y))
for (i in seq_along(y)) { # not vectorized
  nll[i] <- -dmultinom(Y[i, ], size = 1, prob = probs[i, ], log = TRUE)
}
sum(nll)