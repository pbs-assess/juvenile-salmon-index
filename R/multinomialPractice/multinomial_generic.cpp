#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Yobs);
  DATA_MATRIX(cov);
  DATA_VECTOR(k);

  PARAMETER_MATRIX(betas);
  PARAMETER_VECTOR(ints);
  
  Type nll = 0;
  
  matrix<Type> log_odds = beta0 + x*beta1;
  nll -= sum(dnorm(NL,mu,sigma,true));

  return(nll);
}
