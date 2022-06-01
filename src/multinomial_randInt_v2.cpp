#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(y_obs);
  DATA_MATRIX(fx_cov);   // model matrix for fixed effects
  DATA_IVECTOR(rfac);    // vector of random factor levels
  DATA_INTEGER(n_rfac);  // number of random factor levels
  DATA_MATRIX(pred_cov);    // model matrix of predictions
  
  int n_obs = y_obs.rows();           // number of observations
  int n_cat = y_obs.cols();           // number of categories
  int n_preds = pred_cov.rows();         // number of predictions
  
  // Parameters
  PARAMETER_VECTOR(z_rfac);  // vector of random intercepts
  PARAMETER_MATRIX(z_ints);  // matrix of fixed intercepts (rows = fixed cov, cols = k-1)
  PARAMETER(log_sigma_rfac); // among random intercept SD

  // Matrices for storing intermediate objects
  matrix<Type> log_odds(n_obs, (n_cat - 1));
  matrix<Type> exp_log_odds(n_obs, (n_cat - 1));
  vector<Type> denom(n_obs);
  matrix<Type> probs(n_obs, n_cat);
  matrix<Type> logit_probs(n_obs, n_cat);
  
  // Calculate log-odds, then probabilities for each group
  matrix<Type> fx_eff = fx_cov * z_ints;
  for (int k = 0; k < (n_cat - 1); ++k) {
    for (int i = 0; i < n_obs; ++i) {
      log_odds(i, k) = fx_eff(i, k) + z_rfac(rfac(i)); // add random intercept here
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

  // Probability of random intercepts
  for (int h = 0; h < n_rfac; h++) {
    jnll -= dnorm(z_rfac(h), Type(0.0), exp(log_sigma_rfac), true);
  }

  Type sigma_rfac = exp(log_sigma_rfac);
  ADREPORT(sigma_rfac);

   // Calculate log-odds, then probabilities for each fixed effects group
  matrix<Type> pred_log_odds = pred_cov * z_ints;
  matrix<Type> pred_exp_log_odds = exp(pred_log_odds.array());
  
  vector<Type> pred_denom(n_preds);
  for (int ii = 0; ii < n_preds; ++ii) {
    Type sum_exp_log_odds = 0.;
    for (int k = 0; k < (n_cat - 1); ++k) {
      sum_exp_log_odds += pred_exp_log_odds(ii, k);
    }
    pred_denom(ii) = 1. + sum_exp_log_odds;
  }

  matrix<Type> pred_probs(n_preds, n_cat);
  for (int g = 0; g < n_cat; ++g) {
    if (g < (n_cat - 1)) {
      for (int ii = 0; ii < n_preds; ++ii) {
        pred_probs(ii, g) = pred_exp_log_odds(ii, g) / pred_denom(ii);
      }
    } else if (g == (n_cat - 1)) {
      for (int ii = 0; ii < n_preds; ++ii) {
        Type summed_probs = 0;
        for (int k = 0; k < (n_cat - 1); ++k) {
          summed_probs += pred_probs(ii, k);
        }
        pred_probs(ii, g) = 1. - summed_probs;
      }
    }
  }
  REPORT(pred_probs);
  ADREPORT(pred_probs);

  return jnll;
}
