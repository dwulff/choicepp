#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double submean(std::vector<double> x, int a, int b) {
  double m = 0;
  int n = x.size();
  for(int i = a; i < b && i < n; ++i){
    m += x[i];
    }
  return m / double(b - a);
  }

////////////////////////////////////////////////////////////////////////
//
//      RECENCY: WITHIN OPTION
//
////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericVector recency_wos(GenericVector oo, int choice){
  NumericVector res(2);
  std::vector<double> s0 = oo[0];
  std::vector<double> s1 = oo[1];
  int n0  = s0.size();
  int n1  = s1.size();
  int nh0 = std::ceil(n0 / 2.);
  int nh1 = std::ceil(n1 / 2.);
  double m00,m01,m10,m11;
  if(n0 % 2 == 0){
    m00 = submean(s0,0,nh0);
    m01 = submean(s0,nh0,n0);
    } else {
    m00 = submean(s0,0,nh0);
    m01 = submean(s0,nh0-1,n0);
    }
  if(n1 % 2 == 0){
    m10 = submean(s1,0,nh1);
    m11 = submean(s1,nh1,n1);
    } else {
    m10 = submean(s1,0,nh1);
    m11 = submean(s1,nh1-1,n1);
    }
  int pri, rec;
  if(m00 != m10 && choice != -1){
    if((m00 > m10 && choice == 0) || (m00 < m10 && choice == 1)){
      pri = 1;
    } else {
      pri = 0;
    }
  } else {
    pri = -1;
  }
  if(m01 != m11 && choice != -1){
    if((m01 > m11 && choice == 0) || (m01 < m11 && choice == 1)){
      rec = 1;
    } else {
      rec = 0;
    }
  } else {
    rec = -1;
  }
  res[0] = pri;
  res[1] = rec;
  return res;
}


//' Compute recency effect
//'
//' @export
// [[Rcpp::export]]
NumericMatrix recency(GenericVector ss, NumericVector choices) {
  int np = ss.size();
  NumericMatrix res(np,2);
  for(int p = 0; p < np; p++){
    res(p,_) = recency_wos(ss[p],choices[p]);
    }
  return res;
  }
