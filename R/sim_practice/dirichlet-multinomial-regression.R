# Simualte dirichlet-multinomial data and recover pars using log-likelihood 
# adapted from:
#https://github.com/carriegu/STAT-840/blob/master/DMRegression/R/DM_likelihood.R


library(tidyverse)
library(dirmult) #useful for simulating data

J <- 4 #n categories (e.g. stocks)
P <- 3 #n strata (or covariates - 1 if some continuous)
n <- 50 #n observations (e.g. total sampling events)

# input data frame and design matrix
dat <- data.frame(strata = sample(1:P, n, replace = T)) %>% 
  mutate(strata_f = as.factor(strata))
X <- model.matrix(~strata_f, dat)
N <- sample(c(5:100), n) #sample size per observation

# simulate initial beta matrix to simulate Y
beta0 <- matrix(rnorm((ncol(X)) * J), ncol(X), J)
# reparametrization for data generation
Gamma = exp(X %*% beta0) #fixed effects
Gamma_plus = apply(Gamma, 1, sum) #sum of fixed_effects
theta = 1/(Gamma_plus + 1)
pi = apply(Gamma, 2, function(x) {x / theta})

# Simulate the responses Y from package 'dirmult'
set.seed(123)
Y = simPop(J = 1, n = N[1], pi = pi[1,], theta = theta[1])$data
for(jj in c(2:n)){
  Y = rbind(Y, simPop(J = 1, n = N[jj], pi = pi[jj,], theta = theta[jj])$data)
}

# dir-mult neg loglikelihood
DM_negloglikelihood = function(beta, X, Y){
  
  ll = 0
  
  # Extract some parameters
  n = nrow(Y)
  J = ncol(Y)
  P = ncol(X) - 1
  
  # Reshape beta to a matrix
  beta = matrix(beta, nrow = P+1, ncol = J)
  
  Gamma = exp(X %*% beta)
  n = apply(Y, 1, sum)
  Gamma_plus = apply(Gamma, 1, sum)
  for(i in c(1:nrow(X))){
    ll = ll + lgamma(n[i]+1) + lgamma(Gamma_plus[i]) - lgamma(n[i]+Gamma_plus[i])
    for(j in c(1:ncol(Y))){
      ll = ll + lgamma(Y[i,j] + Gamma[i,j]) - lgamma(Gamma[i,j]) - lgamma(Y[i,j]+1)
    }
  }
  return(-ll)
}

beta_in <- rep(0, ncol(X) * J)
m2 <- nlminb(beta_in, DM_negloglikelihood, X = X, Y = Y)
m2
