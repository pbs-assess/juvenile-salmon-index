#include <TMB.hpp>

template<class Type>
Type invlogit_p1(Type x){
  return 1.0 / (1.0 + exp(-x)) + 1.0;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Abundance data
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(factor1k_i);
  DATA_INTEGER(nk1);
  DATA_MATRIX(X1_pred_ij);
  // Composition data
  DATA_MATRIX(y2_ig);		// matrix of observed distribuons for g (groups - 1)
  DATA_MATRIX(X2_ij);      	// model matrix for fixed effects
  DATA_IVECTOR(factor2k_i); // vector of random factor levels
  DATA_INTEGER(nk2); 		// number of factor levels
  DATA_IVECTOR(m2_all_fac); // combined factor levels (fix only)
  DATA_IVECTOR(m2_fac_key); // unique combined factor levels (fix only)

  // Abundance parameters
  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_phi);
  PARAMETER(logit_p);
  PARAMETER_VECTOR(z1_k);
  PARAMETER(log_sigma_zk1);
  // Composition parameters
  PARAMETER_MATRIX(b2_jg); 	 // matrix of fixed int. (rows = fixed cov, cols = g)
  PARAMETER_VECTOR(z2_k); 	 // vector of random int.
  PARAMETER(log_sigma_zk2);  // among random int SD

  // Abundance intermediate variables
  int n1 = y1_i.size();
  vector<Type> linear_predictor1_i(n1);
  // Composition intermediate variables
  int n2 = y2_ig.rows(); 		// number of observations
  int n_cat = y2_ig.cols(); 		// number of categories
  int n_fix_cov = X2_ij.cols(); 	// number of types of fixed covariates 
  int n_fac_comb = m2_fac_key.size();  // number of unique factor combinations
  matrix<Type> fx_eff(n2, n_fix_cov);
  matrix<Type> log_odds(n2, (n_cat - 1));
  matrix<Type> exp_log_odds(n2, (n_cat - 1));
  vector<Type> denom(n2);
  matrix<Type> probs(n2, n_cat);
  matrix<Type> exp_log_odds_fe(n2, (n_cat - 1));   // identical to above for FE pred
  vector<Type> denom_fe(n2);
  matrix<Type> probs_fe(n2, n_cat);
  matrix<Type> logit_probs_fe(n2, n_cat);


  // Abundance linear predictor
  linear_predictor1_i = X1_ij * b1_j;

  // Composition linear predictor
  // Calculate log-odds, then probabilities for each group
  for (int g = 0; g < (n_cat - 1); ++g) {
    for (int j = 0; j < n_fix_cov; ++j) {
      for (int i = 0; i < n2; ++i) {
        fx_eff(i, j) = b2_jg(j, g) * X2_ij(i, j); // calculate fixed effects
      }
    }
    for (int i = 0; i < n2; ++i) {
      Type linear_predictor2_ig = 0;
      for (int j = 0; j < n_fix_cov; ++j) {
        linear_predictor2_ig += fx_eff(i, j);
      }
      log_odds(i, g) = linear_predictor2_ig + z2_k(factor2k_i(i)); // add random intercept here
      exp_log_odds(i, g) = exp(log_odds(i, g));
      exp_log_odds_fe(i, g) = exp(linear_predictor2_ig); // fixed effects only
    }
  }

  for (int i = 0; i < n2; ++i) {
    Type sum_exp_log_odds = 0.;
    Type sum_exp_log_odds_fe = 0.;
    for (int g = 0; g < (n_cat - 1); ++g) {
      sum_exp_log_odds += exp_log_odds(i, g);
      sum_exp_log_odds_fe += exp_log_odds_fe(i, g);
    }
    denom(i) = 1. + sum_exp_log_odds;
    denom_fe(i) = 1. + sum_exp_log_odds_fe;
  }

  for (int gg = 0; gg < n_cat; ++gg) {
    if (gg < (n_cat - 1)) {
      for (int i = 0; i < n2; ++i) {
        probs(i, gg) = exp_log_odds(i, gg) / denom(i);
        probs_fe(i, gg) = exp_log_odds_fe(i, gg) / denom_fe(i);
      }
    } else if (gg == (n_cat - 1)) {
      for (int i = 0; i < n2; ++i) {
        Type summed_probs = 0;
        Type summed_probs_fe = 0;
        for (int g = 0; g < (n_cat - 1); ++g) {
          summed_probs += probs(i, g);
          summed_probs_fe += probs_fe(i, g);
        }
        probs(i, gg) = 1. - summed_probs;
        probs_fe(i, gg) = 1. - summed_probs_fe;
      }
    }
    for (int i = 0; i < n2; ++i) {
      logit_probs_fe(i, gg) = logit(probs_fe(i, gg)); 
    }
  }


  Type jnll = 0.0; // initialize joint negative log likelihood

  // Abundance nll
  for(int i = 0; i < n1; i++){
    jnll -= dtweedie(y1_i(i),
        exp(linear_predictor1_i(i) + z1_k(factor1k_i(i))),
        exp(log_phi), invlogit_p1(logit_p), true);
  }
  // Composition nll
  for (int i = 0; i < n2; i++) {
    jnll -=
        dmultinom(vector<Type>(y2_ig.row(i)), vector<Type>(probs.row(i)), true);
  }
  // Probability of random composition coefficients
  for (int k = 0; k < nk1; k++){
    jnll -= dnorm(z1_k(k), Type(0.0), exp(log_sigma_zk1), true);
  }
  for (int k = 0; k < nk2; k++) {
    jnll -= dnorm(z2_k(k), Type(0.0), exp(log_sigma_zk2), true);
  }

  
  // Abundance fixed effect predictions
  vector<Type> log_pred_abund(X1_pred_ij.rows());
  log_pred_abund = X1_pred_ij * b1_j;

  REPORT(log_pred_abund);
  ADREPORT(log_pred_abund);

  // Composition fixed effect predictions
  matrix<Type> logit_pred_prob(n_fac_comb, n_cat);
  vector<int> temp_index(n_fac_comb);    

  for (int m = 0; m < n_fac_comb; ++m) {
    for (int i = 0; i < n2; ++i) {
      if (m2_all_fac(i) == m2_fac_key(m)) {
        temp_index(m) = i;
        break;
      }
    }
    for (int g = 0; g < n_cat; ++g) {
      logit_pred_prob(m, g) = logit_probs_fe(temp_index(m), g);
    }
  }
  
  REPORT(logit_pred_prob);
  ADREPORT(logit_pred_prob);

  // Combined predictions
  matrix<Type> pred_abund_mg(n_fac_comb, n_cat);

  for (int m = 0; m < n_fac_comb; ++m) {
  	for (int g = 0; g < n_cat; ++g) {
  		pred_abund_mg(m, g) = exp(log_pred_abund(m)) * invlogit(logit_pred_prob(m, g));
  	}
  }

  REPORT(pred_abund_mg);
  ADREPORT(pred_abund_mg);

  return jnll;
}
