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
