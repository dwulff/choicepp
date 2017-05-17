#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


// [[Rcpp::export]]
double dbinomC200(int k, int n, double p){
  return rootChooseLookup200(n, k) * std::pow(p,double(k)) * std::pow(1-p,double(n-k));
  }


// [[Rcpp::export]]
double dbinomC(int k, int n, double p, int nlookup, double root){
  std::vector<double> tab = lookup(nlookup,root);
  return rootChooseLookup(n, k, tab) * std::pow(p,double(k)) * std::pow(1-p,double(n-k));
  }
