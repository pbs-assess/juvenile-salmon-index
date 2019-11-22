#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_MATRIX(Yobs);
  DATA_MATRIX(cov);

  PARAMETER_MATRIX(betas);

  int N = Yobs.rows();
  int k = Yobs.cols();
  int m = cov.cols();

  Type jnll = 0.;

  matrix<Type> exp_log_odds(N, (k - 1));
  matrix<Type> log_odds(N, (k - 1));
  matrix<Type> covEffects(N, m);
  matrix<Type> probs(N, k);
  vector<Type> denom(N);

  for (int h = 0; h < (k - 1); ++h) {
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < N; ++i) {
        log_odds(i, j) = betas(j, h) * cov(i, j);
        exp_log_odds(i, j) = exp(log_odds(i, j));
      }
    }
  }

  for (int i = 0; i < N; ++i) {
    Type sumExpLogOdds = 0.;
    for (int h = 0; h < (k - 1); ++h) {
      sumExpLogOdds += exp_log_odds(i, h);
    }
    denom(i) = 1. + sumExpLogOdds;
  }

  for (int g = 0; g < k; ++g) {
    if (g < (k - 1)) {
      for (int i = 0; i < N; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom(i);
      }
    } else if (g == (k - 1)) {
      for (int i = 0; i < N; ++i) {
        Type summedProbs = 0;
        for (int h = 0; h < (k - 1); ++h) {
          summedProbs += probs(i, h);
        }
        probs(i, g) = 1. - summedProbs;
      }
    }
  }

  for (int i = 0; i < N; i++) {
    jnll -=
        dmultinom(vector<Type>(Yobs.row(i)), vector<Type>(probs.row(i)), true);
  }

  REPORT(probs);
  REPORT(log_odds);
  ADREPORT(log_odds);
  return jnll;
}
