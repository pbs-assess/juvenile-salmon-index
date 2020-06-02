#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data inputs
  DATA_MATRIX(Y); // response matrix (a n-by-J matrix)
  DATA_MATRIX(X); // covariate matrix (a n-by-(P+1) matrix, first column is 1 to account for intercept)
  // parameter inputs
  PARAMETER_MATRIX(beta); // parameter matrix
  
  // The dimensions of data
  Type J;
  J = Y.array().cols();
  REPORT(J);
  Type P;
  P = X.array().cols();
  P = P - 1;
  Type n;
  n = X.array().rows();
  
  // calculate negative loglikelihood (a (P+1)-by-J matrix, first row is intercept parameter)
  matrix<Type> XB = X * beta;
  matrix<Type> gamma = exp(XB.array());
  vector<Type> n_plus = Y.rowwise().sum(); // row sum of response
  vector<Type> gamma_plus = gamma.rowwise().sum(); // row sum of gamma
  
  Type ll = 0;
  int i;
  int j;
  for(i=0; i<=(n-1); i++){
    ll = ll + lgamma((n_plus(i)+1));
    ll = ll + lgamma(gamma_plus(i));
    ll = ll - lgamma((n_plus(i)+gamma_plus(i)));
    for(j=0; j<=(J-1); j++){
      ll += lgamma((Y(i,j)+gamma(i,j)));
      ll -= lgamma(gamma(i,j));
      ll -= lgamma((Y(i,j)+1));
    }
  }
  
  
  Type negloglikelihood;
  negloglikelihood = -ll;
  
  // Return negative loglikelihood
  return negloglikelihood;
  
}

