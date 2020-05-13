#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () 
{
  // data inputs
  DATA_MATRIX(y_obs);    // response matrix (a n-by-k matrix)
  DATA_MATRIX(fx_cov);   // covariate matrix (a n-by-j matrix, first column is 1 to account for intercept)
  // DATA_IVECTOR(rfac);    // vector of random factor levels
  // DATA_INTEGER(n_rfac);  // number of random factor levels
  // DATA_IVECTOR(all_fac); // vector of combined factor levels (fix and rand.)
  // DATA_IVECTOR(fac_key); // vector of unique combined factor levels (fix and rand.)
  DATA_MATRIX(pred_cov);    // model matrix fo predictions

  // parameter inputs
  PARAMETER_MATRIX(z_ints); // parameter matrix
  // PARAMETER_VECTOR(z_rfac);  // vector of random intercepts
  // PARAMETER(log_sigma_rfac); // among random intercept SD

  // The dimensions of data
  int n_obs = y_obs.rows();         // number of observations
  int n_cat = y_obs.cols();         // number of categories
  // int n_fix_cov = fx_cov.cols();    // number of types of fixed covariates in mm
  int n_levels = pred_cov.rows();   // number of covariates to make predictions on

  // calculate fixed effects
  matrix<Type> fx_eff = fx_cov * z_ints;
  matrix<Type> gamma = exp(fx_eff.array());
  vector<Type> n_plus = y_obs.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  
  Type jll = 0; // initialize joint log-likelihood
  for(int i = 0; i <= (n_obs - 1); i++){
    jll = jll + lgamma((n_plus(i) + 1));
    jll = jll + lgamma(gamma_plus(i));
    jll = jll - lgamma((n_plus(i) + gamma_plus(i)));
    for(int k = 0; k <= (n_cat - 1); k++){
      jll += lgamma((y_obs(i, k) + gamma(i, k)));
      jll -= lgamma(gamma(i, k));
      jll -= lgamma((y_obs(i, k) + 1));
    }
  }
  
  // calculate predictions
  matrix<Type> pred_eff(n_levels, n_cat);    //pred effects on log scale
  matrix<Type> pred_gamma(n_levels, n_cat);  //transformed pred effects 
  vector<Type> pred_gamma_plus(n_levels);        
  vector<Type> pred_theta(n_levels); 
  matrix<Type> pred_pi(n_levels, n_cat);      // predicted counts in real 
  vector<Type> pred_n_plus(n_levels); 
  matrix<Type> pred_pi_prop(n_levels, n_cat); // predicted counts as ppn.

  pred_eff = pred_cov * z_ints; 
  pred_gamma = exp(pred_eff.array());
  pred_gamma_plus = pred_gamma.rowwise().sum();
  pred_theta = 1 / (pred_gamma_plus + 1);
  for(int m = 0; m < n_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
    }
  }
  pred_n_plus = pred_pi.rowwise().sum();
  for(int m = 0; m < n_levels; m++) {
    for(int k = 0; k < n_cat; k++) {
      pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
    }
  }

  REPORT(pred_gamma);
  ADREPORT(pred_gamma);
  REPORT(pred_pi);
  ADREPORT(pred_pi);
  REPORT(pred_pi_prop);
  ADREPORT(pred_pi_prop);
  
  Type jnll;
  jnll = -jll;
  
  // Return negative loglikelihood
  return jnll;
  
}

