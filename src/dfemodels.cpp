#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


//' Compute natural mean
//'
//' @export
// [[Rcpp::export]]
int nm(GenericVector oo) {
  std::vector<double> s0 = oo[0];
  std::vector<double> s1 = oo[1];
  double m0 = 0, m1 = 0;
  std::vector<double>::const_iterator it;
  for(it = s0.begin(); it != s0.end(); ++it) m0 += *it;
  for(it = s1.begin(); it != s1.end(); ++it) m1 += *it;
  int choice;
  if(m0 != m1){
    if(m0 > m1){
      choice = 0;
      } else {
      choice = 1;
      }
    } else {
    choice = -1;
    }
  return choice;
  }

//' Compute natural mean
//'
//' @export
// [[Rcpp::export]]
int nm_rand(GenericVector oo, double phi) {
  std::vector<double> s0 = oo[0];
  std::vector<double> s1 = oo[1];
  double m0 = 0, m1 = 0, p, r;
  std::vector<double>::const_iterator it;
  for(it = s0.begin(); it != s0.end(); ++it) m0 += *it;
  for(it = s1.begin(); it != s1.end(); ++it) m1 += *it;
  p = choice_rule(m0,m1,phi);
  r = double(std::rand()) / RAND_MAX;
  int choice;
  if(r < p){
    choice = 0;
    } else {
    choice = 1;
    }
  return choice;
  }


//' Compute natural mean ignoring first experiences
//'
//' @export
// [[Rcpp::export]]
int nm_rec_rand(GenericVector oo, double phi, int ignore) {
  std::vector<double> s0 = oo[0];
  std::vector<double> s1 = oo[1];
  double m0 = 0, m1 = 0, p, r;
  int ignore0 = ignore, ignore1 = ignore;
  if(ignore > s0.size()) ignore0 = s0.size();
  if(ignore > s1.size()) ignore1 = s1.size();
  std::vector<double>::const_iterator it;
  for(it = s0.begin() + ignore0; it != s0.end(); ++it) m0 += *it;
  for(it = s1.begin() + ignore1; it != s1.end(); ++it) m1 += *it;
  p = choice_rule(m0,m1,phi);
  r = double(std::rand()) / RAND_MAX;
  int choice;
  if(r < p){
    choice = 0;
  } else {
    choice = 1;
  }
  return choice;
}


//' Compute natural mean
//'
//' @export
// [[Rcpp::export]]
NumericVector nms(GenericVector ss, double phi) {
  int np = ss.size();
  NumericVector res(np);
  for(int p = 0; p < np; p++){
    res[p] = nm_rand(ss[p], phi);
    }
  return res;
  }

//' Compute natural mean
//'
//' @export
// [[Rcpp::export]]
NumericVector nms_rec(GenericVector ss, double phi, std::vector<int> ignore) {
  int np = ss.size();
  NumericVector res(np);
  if(ignore.size() != np){
    std::vector<int> ignore_n(np);
    for(int i = 0; i < np; ++i) ignore_n[i] = ignore[1];
    ignore  = ignore_n;
    }
  for(int p = 0; p < np; p++){
    res[p] = nm_rec_rand(ss[p], phi, ignore[p]);
  }
  return res;
}
