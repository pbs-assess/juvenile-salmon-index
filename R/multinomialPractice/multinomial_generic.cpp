#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Yobs);
  DATA_MATRIX(cov);
  DATA_NUMERIC(k);

  PARAMETER_MATRIX(betas);
  PARAMETER_VECTOR(ints);
  
  Type nll = 0;
  
  NumericMatrix log_odds(N, (k-1));
  for(int i = 0; i < k; ++i) {
    log_odds( , i)
  }

  matrix<Type> log_odds = beta0 + x*beta1;
  nll -= sum(dnorm(NL,mu,sigma,true));

  return(nll);
}


NumericMatrix exp_log_odds(N, (k - 1));
  NumericVector denom(N);
  NumericMatrix probs(N, k);
  
  for(int h = 0; h < (k - 1); ++h) {
    NumericMatrix covEffects(N, m);

    for(int j = 0; j < m; ++j) {
      for(int i = 0; i < N; ++i) {
        covEffects(i, j) = betas(j, h) * cov(i, j);
      }
    }
    
    for(int i = 0; i < N; ++i) {
      double sumCovEff = 0;
      for(int j = 0; j < m; ++j) {
        sumCovEff += covEffects(i, j); 
      }
      exp_log_odds(i, h) = exp(ints[h] + sumCovEff);
    }
  }
  
  for(int i = 0; i < N; ++i) {
    double sumExpLogOdds = 0;
    for(int h = 0; h < (k - 1); ++h) {
      sumExpLogOdds += exp_log_odds(i, h);
    }
    denom[i] = 1 + sumExpLogOdds;
  }
  
  for(int g = 0; g < k; ++g) {
    if (g < (k - 1)) {
      for(int i = 0; i < N; ++i) {
        probs(i, g) = exp_log_odds(i, g) / denom[i];
      }
    } else if (g == (k - 1)) {
      for(int i = 0; i < N; ++i) {
        double summedProbs = 0;
        for (int h = 0; h < (k - 1); ++h) {
          summedProbs += probs(i, h);
        }
        probs(i, g) = 1 - summedProbs;
      } 
    }
  }
  
  return probs;