#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(y_obs);
  DATA_IVECTOR(fac1); // vector of factor levels
  
  int n_obs = y_obs.rows(); // number of observations
  int n_cat = y_obs.cols(); // number of categories
  int n_fac = fac1.size(); // number of factor levels

  // Parameters
  PARAMETER_VECTOR(beta1); // intercepts for k-1 categories
  PARAMETER_VECTOR(z_fac1); // vector of random intercepts
  PARAMETER(log_sigma); // global SD
  PARAMETER(log_sigma_fac1); // among random intercept SD

  // Matrices for storing intermediate objects
  matrix<Type> log_odds(n_obs, (n_cat - 1));
  matrix<Type> exp_log_odds(n_obs, (n_cat - 1));
  matrix<Type> probs(n_obs, n_cat);
  matrix<Type> logit_probs(n_obs, n_cat);
  vector<Type> denom(n_obs);

  // Calculate log-odds, then probabilities
  for (int k = 0; k < (n_cat - 1); ++k) {
    for (int i = 0; i < n_obs; ++i) {
      log_odds(i, k) = beta1(k) + (z_fac1(fac1(i)) * exp(log_sigma));
      exp_log_odds(i, k) = exp(log_odds(i, k));
    }
  }

  for (int i = 0; i < n_obs; ++i) {
    Type sum_exp_log_odds = 0.;
    for (int k = 0; k < (n_cat - 1); ++k) {
      sum_exp_log_odds += exp_log_odds(i, k);
    }
    denom(i) = 1. + sum_exp_log_odds;
  }

  for (int g = 0; g < n_cat; ++g) {
    if (g < (n_cat - 1)) {
      for (int i = 0; i < n_obs; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom(i);
      }
    } else if (g == (n_cat - 1)) {
      for (int i = 0; i < n_obs; ++i) {
        Type summed_probs = 0;
        for (int k = 0; k < (n_cat - 1); ++k) {
          summed_probs += probs(i, k);
        }
        probs(i, g) = 1. - summed_probs;
      }
    }
    for (int i = 0; i < n_obs; ++i) {
      logit_probs(i, g) = logit(probs(i, g));
    }
  }

  Type jnll = 0.; //initialize joint negative log likelihood

  for (int i = 0; i < n_obs; i++) {
    jnll -=
        dmultinom(vector<Type>(y_obs.row(i)), vector<Type>(probs.row(i)), true);
  }

  // Probability of random coefficients
  for (int h = 0; h < n_fac; h++) {
    jnll -= dnorm(z_fac1(h), Type(0.0), exp(log_sigma_fac1), true);
  }

  Type sigma = exp(log_sigma);
  Type sigma_fac1 = exp(log_sigma_fac1);
  REPORT(sigma);
  REPORT(sigma_fac1);

  REPORT(log_odds);
  REPORT(probs);
  REPORT(logit_probs);
  ADREPORT(logit_probs);
  
  return jnll;
}
