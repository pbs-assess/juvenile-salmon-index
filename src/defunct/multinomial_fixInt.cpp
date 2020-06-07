#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(y_obs);
  DATA_MATRIX(fx_cov);   // model matrix for fixed effects
  DATA_IVECTOR(all_fac); // vector of combined factor levels (fix and rand.)
  DATA_IVECTOR(fac_key); // vector of unique combined factor levels (fix and rand.)

  int n_obs = y_obs.rows();           // number of observations
  int n_cat = y_obs.cols();           // number of categories
  int n_fix_cov = fx_cov.cols();      // number of types of fixed covariates in mm
  int n_fac_comb = fac_key.size();    // number of unique factor combinations

  // Parameters
  PARAMETER_MATRIX(z_ints);  // matrix of fixed intercepts (rows = fixed cov, cols = k-1)
  
  // Matrices for storing intermediate objects
  matrix<Type> fx_eff(n_obs, n_fix_cov);
  matrix<Type> log_odds(n_obs, (n_cat - 1));
  matrix<Type> exp_log_odds(n_obs, (n_cat - 1));
  vector<Type> denom(n_obs);
  matrix<Type> probs(n_obs, n_cat);
  matrix<Type> logit_probs(n_obs, n_cat);
  
  // Calculate log-odds, then probabilities for each group
  for (int k = 0; k < (n_cat - 1); ++k) {
    for (int j = 0; j < n_fix_cov; ++j) {
      for (int i = 0; i < n_obs; ++i) {
        fx_eff(i, j) = z_ints(j, k) * fx_cov(i, j); // calculate fixed effects
      }
    }
    for (int i = 0; i < n_obs; ++i) {
      Type sum_fix_eff = 0;
      for (int j = 0; j < n_fix_cov; ++j) {
        sum_fix_eff += fx_eff(i, j);
      }
      log_odds(i, k) = sum_fix_eff;
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


  // Populate output matrix with unique combinations of factor levels for 
  // predictions
  matrix<Type> logit_probs_out(n_fac_comb, n_cat);
  vector<int> temp_index(n_fac_comb);    

  for (int j = 0; j < n_fac_comb; ++j) {
    for (int i = 0; i < n_obs; ++i) {
      if (all_fac(i) == fac_key(j)) {
        temp_index(j) = i;
        break;
      }
    }
    for (int g = 0; g < n_cat; ++g) {
      logit_probs_out(j, g) = logit_probs(temp_index(j), g);
    }
  }
  ADREPORT(logit_probs_out);
 
  return jnll;
}
