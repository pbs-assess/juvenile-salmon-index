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
  // parameter inputs
  PARAMETER_MATRIX(z_ints); // parameter matrix
  // PARAMETER_VECTOR(z_rfac);  // vector of random intercepts
  // PARAMETER(log_sigma_rfac); // among random intercept SD

  // The dimensions of data
  // Type k;
  // k = y_obs.array().cols();
  // REPORT(k);
  // Type j;
  // j = fx_cov.array().cols();
  // j = j - 1;
  // Type n;
  // n = X.array().rows();
  
  int n_obs = y_obs.rows();           // number of observations
  int n_cat = y_obs.cols();           // number of categories
  int n_fix_cov = fx_cov.cols();      // number of types of fixed covariates in mm
  int n_fix_cov_m1 = n_fix_cov - 1;   // number of fixed covs minus intercept

  // calculate fixed effects, first row is intercept parameter)
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
  
  
  Type negloglikelihood;
  negloglikelihood = -jll;
  
  // Return negative loglikelihood
  return negloglikelihood;
  
}

