#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);

  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_sigma);

  int n1 = y1_i.size();
  
  Type jnll = 0.0; // initialize joint negative log likelihood

  // Linear predictor
  vector<Type> linear_predictor1_i(n1);
  linear_predictor1_i = X1_ij * b1_j;

  for(int i = 0; i < n1; i++){
    jnll -= dnorm(y1_i(i),
        linear_predictor1_i(i),
        exp(log_sigma), true);
  }
 
  return jnll;
}
