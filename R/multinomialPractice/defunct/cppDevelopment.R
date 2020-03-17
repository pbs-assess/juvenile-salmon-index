## Practice coding skeleton of multinomial TMB model in CPP

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
N <- nrow(Yobs) # number of observations
k <- ncol(Yobs) # number of groups

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
  
  return(probs)
}

## Example off cpp function that calculates probabilities but not NLL
cppFunction('NumericMatrix logOddsC(NumericVector ints, NumericMatrix betas,
NumericMatrix cov, int N, int k, int m) {
  NumericMatrix exp_log_odds(N, (k - 1));
  NumericVector denom(N);
  NumericMatrix probs(N, k);
  
  for(int h = 0; h < (k - 1); ++h) {
    NumericMatrix covEffects(N, m);

    for(int j = 0; j < m; ++j) {
      for(int i = 0; i < N; ++i) {
        covEffects(i, j) = betas(j, h) * cov(i, j);
      }
    }
    
    for(int i = 0; i < N; ++i) {
      double sumCovEff = 0;
      for(int j = 0; j < m; ++j) {
        sumCovEff += covEffects(i, j); 
      }
      exp_log_odds(i, h) = exp(ints[h] + sumCovEff);
    }
  }
  
  for(int i = 0; i < N; ++i) {
    double sumExpLogOdds = 0;
    for(int h = 0; h < (k - 1); ++h) {
      sumExpLogOdds += exp_log_odds(i, h);
    }
    denom[i] = 1 + sumExpLogOdds;
  }
  
  for(int g = 0; g < k; ++g) {
    if (g < (k - 1)) {
      for(int i = 0; i < N; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom[i];
      }
    } else if (g == (k - 1)) {
      for(int i = 0; i < N; ++i) {
        double summedProbs = 0;
        for (int h = 0; h < (k - 1); ++h) {
          summedProbs += probs(i, h);
        }
        probs(i, g) = 1 - summedProbs;
      } 
    }
  }
  
  return probs;
}')

dum2 <- logOddsC(ints, betas, cov = covMatrix, N = N,  k = k, m = 1)
dum <- logOdds(ints, betas, covMatrix, N, k)


## -----------------------------------------------------------------------------

library(Rcpp)

logOddsC2(y_obs, fac1, beta1 = c(b0, b2), z_fac1 = site_mean_a,
          log_sigma = sd_global, n_obs = nrow(y_obs), n_cat = ncol(y_obs),
          n_fac = length(site_mean_a))

log_odds <- matrix(nrow = nrow(y_obs), ncol = ncol(y_obs) - 1)
for (k in seq_len(ncol(y_obs) - 1)) {  
  for (i in seq_len(nrow(y_obs))) {
    log_odds[i, k] = 1 + site_mean_a[fac1[i]]
  }
}

cppFunction('NumericVector logOddsC2(NumericMatrix y_obs, NumericVector fac1,
NumericVector beta1, NumericVector z_fac1, double log_sigma, 
int n_obs, int n_cat, int n_fac) {
  
  // Matrices for storing intermediate objects
  NumericMatrix log_odds(n_obs, (n_cat - 1));
  NumericMatrix exp_log_odds(n_obs, (n_cat - 1));
  NumericMatrix probs(n_obs, n_cat);
  NumericMatrix logit_probs(n_obs, n_cat);
  NumericVector denom(n_obs);
  NumericVector ests(n_fac);
  
  for (int h = 0; h < n_fac; h++) {
    ests(h) = z_fac1(h);
  }
  
  return(ests);
  
  // Calculate log-odds, then probabilities
  //for (int k = 0; k < (n_cat - 1); ++k) {
  // for (int i = 0; i < n_obs; ++i) {
  //   log_odds(i, k) = beta1(k) + (z_fac1(fac1(i)) * exp(log_sigma));
  //  exp_log_odds(i, k) = exp(log_odds(i, k));
  // }
  //}
//
  //return exp_log_odds;
}')


# ------------------------------------------------------------------------------
# Match and fill vectors

# dummy simulated data from hierarchical multinomial model
dumm <- readRDS(here::here("R", "multinomialPractice", "temp_data.RDS"))
probs_mat <- dumm %>% 
  select(p1:p3) %>% 
  as.matrix()
fac_vec <- dumm %>% 
  pull(facs) %>% 
  as.numeric() 
fac_vec <- fac_vec - 1
fac_key <- sort(unique(fac_vec))

# Fill using R
pred_out <- matrix(NA,
                   nrow = length(fac_key),
                   ncol = (k + 1))
for (h in 1:nrow(pred_out)) {
  indexx <- min(which(fac_key[h] == fac_vec)) 
  for (j in 1:ncol(probs_mat)) {
    pred_out[h, j] <- probs_mat[indexx, j]
  }
}
pred_out

sub_pred(probs_mat, fac_vec, fac_key)

library(Rcpp)
cppFunction('NumericMatrix sub_pred(NumericMatrix probs_mat, 
  NumericVector fac_vec, NumericVector fac_key) {
    int n_fac = fac_key.size();
    int cats = probs_mat.cols();
    int n_obs = fac_vec.size();
    
    NumericMatrix pred_out(n_fac, cats);
    NumericVector indexx(n_fac);    

    for (int f = 0; f < n_fac; ++f) {
      for (int i = 0; i < n_obs; ++i) {
        if (fac_vec(i) == fac_key(f)) {
          indexx(f) = i;
          break;
        }
      }
      for (int j = 0; j < cats; ++j) {
        pred_out(f, j) = probs_mat(indexx(f), j);
      }
    }

    return(pred_out);
}')
