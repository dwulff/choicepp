#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

//' Probability of choosing EV
//'
//' Function calculates the probability of choosing EV for a given a
//' binary problem (between options with maximally 2 outcomes) and sample size based
//' on a strategy that chooses the options with the higher experienced mean.
//'
//' @param p numeric vector containing the outcomes of A, the respective probabilities
//'   of the outcomes of A, the outcomes of B, and the respective probabilities of the
//'   outcomes of B.
//' @param n integer specifying the sample size in terms of the number of samples taken
//'   from \emph{each} of the options.
//'
//' @return a numeric value specifying the probability of choosing the higher EV option.
//'
//' @export
// [[Rcpp::export]]
double pevbyn_nm(NumericVector p, int n) {
  double oA1 = p[0], oA2 = p[2], pA1 = p[1], pA2 = p[3], oB1 = p[4], oB2 = p[6], pB1 = p[5], pB2 = p[7];
  double lik, m0, m1, ev0, ev1, pA, pEV, lik0 = 0, lik1 = 0;
  for(int i = 0; i <= n; ++i){
    for(int j = 0; j <= n; ++j){
      lik = dbinomC200(i,n,pA1) * dbinomC200(j,n,pB1);
      m0 = i * oA1 + (n - i) * oA2;
      m1 = j * oB1 + (n - j) * oB2;
      if(m0 != m1){
        if(m0 > m1){
          lik0 += lik;
          } else{
          lik1 += lik;
          };
        }
      }
    }
  pA = lik0 / (lik0 + lik1);
  ev0 = oA1 * pA1 + oA2 * pA2;
  ev1 = oB1 * pB1 + oB2 * pB2;
  if(ev0 != ev1){
    if(ev0 > ev1){
      pEV = pA;
      } else {
      pEV = 1-pA;
      }
    } else {
    pEV = NA_REAL;
    }
  return pEV;
  }

