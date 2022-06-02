#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // data inputs
  DATA_MATRIX(y_obs);    // response matrix (n-by-k matrix)
  DATA_MATRIX(fx_cov);   // model matrix for covariates
  DATA_MATRIX(pred_cov);    // model matrix of predictions

  // parameter inputs
  PARAMETER_MATRIX(z_ints); // parameter matrix
  PARAMETER(log_phi); // precision parameter


  // The dimensions of data
  int n_obs = y_obs.rows();         // number of observations
  int n_groups = y_obs.cols();      // number of categories
  int n_preds = pred_cov.rows();         // number of predictions
  
  // calculate fixed effects
  matrix<Type> log_y = log(y_obs.array());
  matrix<Type> fx_eff = fx_cov * z_ints;
  matrix<Type> exp_fx_eff = exp(fx_eff.array());

  // add column of 1s (intercepts) to exp_fix_eff matrix
  matrix<Type> mu1(n_obs, n_groups);
  for(int i = 0; i < n_obs; i++) {
    mu1(i, 0) = 1;
    for(int k = 0; k < (n_groups - 1); k++) {
      mu1(i, (k + 1)) = exp_fx_eff(i, k);
    }
  }
  vector<Type> mu1_sum = mu1.rowwise().sum();
  
  // calc rowwise sums
  matrix<Type> ba(n_obs, n_groups);
  matrix<Type> ma(n_obs, n_groups);
  Type phi = exp(log_phi);

  for (int i = 0; i < n_obs; i++) {
    for (int k = 0; k < n_groups; k++) {
      ma(i, k) = mu1(i, k) / mu1_sum(i);
      ba(i, k) = phi * ma(i, k);
    }
  }  

  // log-likelihood
  Type jnll = 0;
  jnll = -1 * n_obs * lgamma(phi);
  for (int i = 0; i < n_obs; i++) {
    for (int k = 0; k < n_groups; k++) {
      jnll += lgamma(ba(i, k));
      jnll -= log_y(i, k) * (ba(i, k) - 1);
    }
  }

  // calculate predictions
  matrix<Type> pred_mu1(n_preds, n_groups);
  matrix<Type> pred_est(n_preds, n_groups);

  matrix<Type> pred_eff = pred_cov * z_ints;
  matrix<Type> exp_pred_eff = exp(pred_eff.array());
  
  // add column of 1s (intercepts) to prediction model matrix
  for(int j = 0; j < n_preds; j++) {
    pred_mu1(j, 0) = 1;
    for(int k = 0; k < (n_groups - 1); k++) {
      pred_mu1(j, (k + 1)) = exp_pred_eff(j, k);
    }
  }
  vector<Type> pred_mu1_sum = pred_mu1.rowwise().sum();
  
  for (int j = 0; j < n_preds; j++) {
    for (int k = 0; k < (n_groups); k++) {
      pred_est(j, k) = pred_mu1(j, k) / pred_mu1_sum(j);
    }
  } 

  REPORT(pred_est);
  ADREPORT(pred_est);
  
  // Return negative loglikelihood
  return jnll;
  
}

