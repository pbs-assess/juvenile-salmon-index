#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  // Abundance 
  DATA_VECTOR(y1_i);
  DATA_MATRIX(X1_ij);
  DATA_IVECTOR(factor1k_i);
  DATA_INTEGER(nk1);
  DATA_MATRIX(X1_pred_ij);
  // Composition 
  DATA_MATRIX(y2_ig);		// matrix of observed distribuons for g (groups - 1)
  DATA_MATRIX(X2_ij);      	// model matrix for fixed effects
  DATA_IVECTOR(factor2k_i); // vector of random factor levels
  DATA_INTEGER(nk2); 		// number of factor levels
  DATA_MATRIX(X2_pred_ij);    // prediction matrix for compositon
  // DATA_IVECTOR(m2_all_fac); // combined factor levels (fix only)
  // DATA_IVECTOR(m2_fac_key); // unique combined factor levels (fix only)


  // PARAMETERS 
  // Abundance 
  PARAMETER_VECTOR(b1_j);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(z1_k);
  PARAMETER(log_sigma_zk1);
  // Composition
  PARAMETER_MATRIX(b2_jg); 	 // matrix of fixed int. (rows = fixed cov, cols = g)
  PARAMETER_VECTOR(z2_k); 	 // vector of random int.
  PARAMETER(log_sigma_zk2);  // among random int SD


  // DEFINE INTERMEDIATES
  // Abundance 
  int n1 = y1_i.size();
  vector<Type> linear_predictor1_i(n1);
  // Composition 
  int n2 = y2_ig.rows(); 		// number of observations
  int n_cat = y2_ig.cols(); 		// number of categories
  int n_preds = X2_pred_ij.rows();         // number of predictions
  // int n_fix_cov = X2_ij.cols(); 	// number of types of fixed covariates 
  // int n_fac_comb = m2_fac_key.size();  // number of unique factor combinations
  // matrix<Type> fx_eff(n2, n_fix_cov);
  matrix<Type> log_odds(n2, (n_cat - 1));
  matrix<Type> exp_log_odds(n2, (n_cat - 1));
  vector<Type> denom(n2);
  matrix<Type> probs(n2, n_cat);
  matrix<Type> logit_probs(n_obs, n_cat);
  // matrix<Type> exp_log_odds_fe(n2, (n_cat - 1));   // identical to above for FE pred
  // vector<Type> denom_fe(n2);
  // matrix<Type> probs_fe(n2, n_cat);
  // matrix<Type> logit_probs_fe(n2, n_cat);


  // LINEAR PREDICTORS
  // Abundance
  linear_predictor1_i = X1_ij * b1_j;
  
  // Composition linear
  // Calculate log-odds, then probabilities for each group
  matrix<Type> fx_eff = fx_cov * z_ints;
  for (int k = 0; k < (n_cat - 1); ++k) {
    for (int i = 0; i < n2; ++i) {
      log_odds(i, k) = fx_eff(i, k) + z_rfac(rfac(i)); // add random intercept here
      exp_log_odds(i, k) = exp(log_odds(i, k)); 
    }
  }

  for (int i = 0; i < n2; ++i) {
    Type sum_exp_log_odds = 0.;
    for (int k = 0; k < (n_cat - 1); ++k) {
      sum_exp_log_odds += exp_log_odds(i, k);
    }
    denom(i) = 1. + sum_exp_log_odds;
  }

  for (int g = 0; g < n_cat; ++g) {
    if (g < (n_cat - 1)) {
      for (int i = 0; i < n2; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom(i);
      }
    } else if (g == (n_cat - 1)) {
      for (int i = 0; i < n2; ++i) {
        Type summed_probs = 0;
        for (int k = 0; k < (n_cat - 1); ++k) {
          summed_probs += probs(i, k);
        }
        probs(i, g) = 1. - summed_probs;
      }
    }
    for (int i = 0; i < n2; ++i) {
      logit_probs(i, g) = logit(probs(i, g)); 
    }
  }


  // LOG-LIKELIHOOD
  Type jnll = 0.0; // initialize joint negative log likelihood

  // Abundance nll
  Type s1, s2;
  for(int i = 0; i < n1; i++){
    s1 = linear_predictor1_i(i) + z1_k(factor1k_i(i)); //mu
    s2 = 2.0 * (s1 - log_phi); //scale
    jnll -= dnbinom_robust(y1_i(i), s1, s2, true);
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

  
  // FIXED EFFECTS PREDICTIONS
  // Abundance 
  vector<Type> log_pred_abund(X1_pred_ij.rows());
  log_pred_abund = X1_pred_ij * b1_j;

  REPORT(log_pred_abund);
  ADREPORT(log_pred_abund);

  // Composition 
  matrix<Type> pred_log_odds = X2_pred_ij * b2_jg;
  matrix<Type> pred_exp_log_odds = exp(pred_log_odds.array());
  
  vector<Type> pred_denom(n_preds);
  for (int ii = 0; ii < n_preds; ++ii) {
    Type sum_exp_log_odds = 0.;
    for (int k = 0; k < (n_cat - 1); ++k) {
      sum_exp_log_odds += pred_exp_log_odds(ii, k);
    }
    pred_denom(ii) = 1. + sum_exp_log_odds;
  }

  matrix<Type> pred_probs(n_preds, n_cat);
  for (int g = 0; g < n_cat; ++g) {
    if (g < (n_cat - 1)) {
      for (int ii = 0; ii < n_preds; ++ii) {
        pred_probs(ii, g) = pred_exp_log_odds(ii, g) / pred_denom(ii);
      }
    } else if (g == (n_cat - 1)) {
      for (int ii = 0; ii < n_preds; ++ii) {
        Type summed_probs = 0;
        for (int k = 0; k < (n_cat - 1); ++k) {
          summed_probs += pred_probs(ii, k);
        }
        pred_probs(ii, g) = 1. - summed_probs;
      }
    }
  }
  REPORT(pred_probs);
  ADREPORT(pred_probs);

  // Combined predictions
  matrix<Type> pred_abund_mg(n_fac_comb, n_cat);

  for (int m = 0; m < n_fac_comb; ++m) {
  	for (int g = 0; g < n_cat; ++g) {
  		pred_abund_mg(m, g) = exp(log_pred_abund(m)) * pred_prob(m, g);
  	}
  }

  REPORT(pred_abund_mg);
  ADREPORT(pred_abund_mg);

  return jnll;
}
