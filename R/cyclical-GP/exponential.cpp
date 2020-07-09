#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  PARAMETER(b0);
  PARAMETER(log_gp_sigma);
  PARAMETER(log_gp_theta);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(eta)
  DATA_INTEGER(flag);
  
  int n = y.size(); 
  
  Type gp_sigma = exp(log_gp_sigma);
  Type gp_theta = exp(log_gp_theta);
  Type sigma = exp(log_sigma);
  REPORT(sigma);
  REPORT(gp_sigma);
  REPORT(gp_theta);
  
  matrix<Type> C(n, n);
  C.setZero();  
  for(int i=0; i<C.rows(); i++)
    for(int j=0; j<C.cols(); j++) {
      Type dist_sq = (x[j] - x[i]) * (x[j] - x[i]);
      C(i, j) = gp_sigma*gp_sigma * exp(- 1 / (2 * gp_theta * gp_theta) * dist_sq);
    }
    
  Type jnll = density::MVNORM(C)(eta);
  if (flag == 0) return jnll;
  jnll-= dnorm(y, vector<Type>(b0 + eta), sigma, true).sum();
  
  jnll-= dnorm(gp_sigma, Type(0), Type(2));
  jnll-= dnorm(gp_theta, Type(0), Type(3));
  jnll-= dnorm(b0, Type(0), Type(10));
  jnll-= dnorm(sigma, Type(0), Type(1));
  jnll-= log_gp_sigma;
  jnll-= log_gp_theta;
  jnll-= log_sigma;
  
  return jnll;
}
