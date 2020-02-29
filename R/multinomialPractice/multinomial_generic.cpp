#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(y_obs);
  DATA_MATRIX(cov);

  PARAMETER_MATRIX(betas);

  int N = y_obs.rows(); // number of observations
  int k = y_obs.cols(); // number of categories
  int m = cov.cols(); // number of covariates

  Type jnll = 0.;

  matrix<Type> exp_log_odds(N, (k - 1));
  matrix<Type> log_odds(N, (k - 1));
  matrix<Type> cov_effects(N, m);
  matrix<Type> probs(N, k);
  matrix<Type> logit_probs(N, k);
  vector<Type> denom(N);

  for (int h = 0; h < (k - 1); ++h) {
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < N; ++i) {
        cov_effects(i, j) = betas(j, h) * cov(i, j);
      }
    }
    for (int i = 0; i < N; ++i) {
      Type sum_cov_effects = 0;
      for (int j = 0; j < m; ++j) {
        sum_cov_effects += cov_effects(i, j);
      }
      log_odds(i, h) = sum_cov_effects;
      exp_log_odds(i, h) = exp(log_odds(i, h));
    }
  }

  for (int i = 0; i < N; ++i) {
    Type sum_exp_log_odds = 0.;
    for (int h = 0; h < (k - 1); ++h) {
      sum_exp_log_odds += exp_log_odds(i, h);
    }
    denom(i) = 1. + sum_exp_log_odds;
  }

  for (int g = 0; g < k; ++g) {
    if (g < (k - 1)) {
      for (int i = 0; i < N; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom(i);
      }
    } else if (g == (k - 1)) {
      for (int i = 0; i < N; ++i) {
        Type summed_probs = 0;
        for (int h = 0; h < (k - 1); ++h) {
          summed_probs += probs(i, h);
        }
        probs(i, g) = 1. - summed_probs;
      }
    }
    for (int i = 0; i < N; ++i) {
      logit_probs(i, g) = logit(probs(i, g));
    }
  }

  for (int i = 0; i < N; i++) {
    jnll -=
        dmultinom(vector<Type>(y_obs.row(i)), vector<Type>(probs.row(i)), true);
  }

  REPORT(log_odds);
  REPORT(probs);
  REPORT(logit_probs);
  ADREPORT(logit_probs);
  
  return jnll;
}
