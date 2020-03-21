#include <TMB.hpp>

// template<class Type>
// Type invlogit_p1(Type x){
//   return 1.0 / (1.0 + exp(-x)) + 1.0;
// }

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);

  // DATA_IVECTOR(factor1k_i);
  // DATA_INTEGER(nk1);

  PARAMETER_VECTOR(b1_j);
  // PARAMETER(log_phi);
  // PARAMETER(logit_p);
  PARAMETER(log_sigma);

  // DATA_MATRIX(X1_pred_ij);

  int n1 = y1_i.size();
  // int n_cov = b1_j.size();

  Type jnll = 0.0; // initialize joint negative log likelihood

  // Linear predictor
  vector<Type> linear_predictor1_i(n1);
  linear_predictor1_i = X1_ij * b1_j;

  // for (int i = 0; i < n1; ++i) {
  //   Type dum_fix = 0;
  //   for (int j = 0; j < n_cov; ++j) {
  //     dum_fix += b1_j(j) * X1_ij(i, j);
  //   }
  //   linear_predictor1_i(i) = dum_fix;
  // }

  for(int i = 0; i < n1; i++){
    jnll -= dnorm(y1_i(i),
        linear_predictor1_i(i),
        exp(log_sigma), true);
  }
  // vector<Type> log_prediction(X1_pred_ij.rows());
  // log_prediction = X1_pred_ij * b1_j;

  // REPORT(log_prediction);
  // ADREPORT(log_prediction);

  return jnll;
}
