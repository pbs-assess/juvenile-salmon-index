#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // data inputs
  DATA_MATRIX(y_obs);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(fx_cov);   // covariate matrix (a n-by-j matrix, first column is 1 to account for intercept)
  DATA_MATRIX(pred_cov);    // model matrix of predictions

  // parameter inputs
  PARAMETER_MATRIX(z_ints); // parameter matrix
  PARAMETER(log_phi); // precision parameter


  // The dimensions of data
  int n_obs = y_obs.rows();         // number of observations
  // int n_cat = y_obs.cols();         // number of categories
  int n_fix_cov = fx_cov.cols();    // number of types of fixed covariates in m
  // int n_levels = pred_cov.rows();   // number of covariates to make predictions on

  
  // calculate fixed effects
  matrix<Type> log_y = log(y_obs.array());
  matrix<Type> fx_eff = fx_cov * z_ints;
  matrix<Type> exp_fx_eff = exp(fx_eff.array());

  // add column of 1s (intercepts) to exp_fix_eff matrix
  matrix<Type> mu1(n_obs, (n_fix_cov + 1));
  for(int i = 0; i < n_obs; i++) {
    mu1(i, 0) = 1;
    for(int m = 0; m < n_fix_cov; m++){
      mu1(i, (m + 1)) = exp_fx_eff(i, m);
    }
  }
  vector<Type> mu1_sum = mu1.rowwise().sum();
  
  // calc rowwise sums
  matrix<Type> ma(n_obs, (n_fix_cov + 1));
  matrix<Type> ba(n_obs, (n_fix_cov + 1));
  Type phi = exp(log_phi);

  for (int i = 0; i < n_obs; i++) {
    for (int m = 0; m < (n_fix_cov + 1); m++) {
      ma(i, m) = mu1(i, m) / mu1_sum(i);
      ba(i, m) = phi * ma(i, m);
    }
  }  

  // log-likelihood
  Type jnll = 0;
  jnll = -1 * n_obs * lgamma(phi);
  for (int i = 0; i < n_obs; i++) {
    for (int m = 0; m < (n_fix_cov + 1); m++) {
      jnll += lgamma(ba(i, m));
      jnll -= log_y(i, m) * (ba(i, m) - 1);
    }
  }

  // calculate predictions
  matrix<Type> pred_eff = pred_cov * z_ints;
  matrix<Type> exp_pred_eff = exp(pred_eff.array());

  // add column of 1s (intercepts) to prediction model matrix
  matrix<Type> pred_mu1(n_obs, (n_fix_cov + 1));
  for(int i = 0; i < n_obs; i++) {
    pred_mu1(i, 0) = 1;
    for(int m = 0; m < n_fix_cov; m++){
      pred_mu1(i, (m + 1)) = exp_pred_eff(i, m);
    }
  }
  vector<Type> pred_mu1_sum = pred_mu1.rowwise().sum();
  
  for (int i = 0; i < n_obs; i++) {
    for (int m = 0; m < (n_fix_cov + 1); m++) {
      pred_est(i, m) = pred_mu1(i, m) / pred_mu1_sum(i);
    }
  } 

  REPORT(pred_gamma);
  ADREPORT(pred_gamma);
  
  Return negative loglikelihood
  return jnll;
  
}

