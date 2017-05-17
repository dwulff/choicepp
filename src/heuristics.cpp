#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector reduce_option(NumericVector opt, bool has_n = false){
  int i, no ;
  if(has_n == true){
    no = (opt.size() - 2)/2;
    } else {
    no = (opt.size() - 1)/2;
    }
  std::vector<int> indices;
  for(i = 0; i < no; ++i){
    if(opt[no + i] > 0) indices.push_back(i);
    }
  int new_no = indices.size();
  NumericVector new_opt(new_no * 2);
  for(i = 0; i < new_no; ++i){
    new_opt[i]          = opt[indices[i]];
    new_opt[i + new_no] = opt[indices[i] + no];
    }
  return new_opt;
  }


// [[Rcpp::export]]
std::vector<double> sort_by_p(NumericVector opt, bool has_n = true){
  opt = reduce_option(opt, has_n);
  int no = opt.size() / 2;
  std::vector< std::pair<double, double> > option;
  std::pair<double, double> event;
  for(int i = 0; i < no; i++){
    event.first  = opt[i];
    event.second = opt[no + i];
    option.push_back(event);
    }
  option = mysort(option, true);
  std::vector< std::pair<double, double> >::const_iterator it;
  std::vector<double> sorted_out;
  for (it = option.begin(); it != option.end(); ++it) sorted_out.push_back(it->first);
  return sorted_out;
  }



// [[Rcpp::export]]
int bayes_mean(NumericVector opt_A, NumericVector opt_B){
  int n_A = opt_A[opt_A.size()-1], n_B = opt_B[opt_B.size()-1];
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double A = 0, B = 0;
  double prior_A = 1.0 / n_A, prior_B = 1.0 / n_B;
  for(int i = 0; i < no_A; ++i){
    A += opt_A[i] * ((opt_A[no_A + i] + prior_A) / (1 + no_A * prior_A));
    }
  for(int i = 0; i < no_B; ++i){
    B += opt_B[i] * ((opt_B[no_B + i] + prior_B) / (1 + no_B * prior_B));
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int natural_mean(NumericVector opt_A, NumericVector opt_B){
  int no_A = (opt_A.size() - 2)/2, no_B = (opt_B.size() - 2)/2;
  double A = 0, B = 0;
  for(int i = 0; i < no_A; ++i){
    A += opt_A[i] * opt_A[no_A + i];
    }
  for(int i = 0; i < no_B; ++i){
    B += opt_B[i] * opt_B[no_B + i];
    }
  int choice;
  if(A == B){
    return -1;
  } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }


// [[Rcpp::export]]
int equiprobable(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double A = 0, B = 0;
  for(int i = 0; i < no_A; ++i) A += opt_A[i];
  for(int i = 0; i < no_B; ++i) B += opt_B[i];
  A = A / no_A;
  B = B / no_B;
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int minimax(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double A = opt_A[0], B = opt_B[0];
  if(no_A > 1){
    for(int i = 1; i < no_A; ++i) if(opt_A[i] < A) A = opt_A[i];
    }
  if(no_B > 1){
    for(int i = 1; i < no_B; ++i) if(opt_B[i] < B) B = opt_B[i];
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int maximax(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double A = opt_A[0], B = opt_B[0];
  if(no_A > 1){
    for(int i = 1; i < no_A; ++i) if(opt_A[i] > A) A = opt_A[i];
    }
  if(no_B > 1){
    for(int i = 1; i < no_B; ++i) if(opt_B[i] > B) B = opt_B[i];
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int minimax_regret(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double max_A = opt_A[0], max_B = opt_B[0], min_A = opt_A[0], min_B = opt_A[0];
  if(no_A > 1){
    for(int i = 1; i < no_A; ++i){
      if(opt_A[i] > max_A) max_A = opt_A[i];
      if(opt_A[i] < min_A) min_A = opt_A[i];
      }
    }
  if(no_B > 1){
    for(int i = 1; i < no_B; ++i){
      if(opt_B[i] > max_B) max_B = opt_B[i];
      if(opt_B[i] < min_B) min_B = opt_B[i];
      }
    }
  double A = min_A - max_A;
  double B = min_B - max_B;
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A < B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int payoffelimination(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double A = 0, B = 0;
  int min_no = no_A;
  if(no_A > no_B) min_no = no_B;
  std::vector<double> outs_A, outs_B;
  for(int i = 0; i < min_no; ++i){
    if(opt_A[i] != opt_B[i]){
      if(opt_A[i] > opt_B[i]){
        A = 1;
        } else {
        B = 1;
        }
      break;
      }
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }


// [[Rcpp::export]]
int betterthanaverage(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  int A = 0, B = 0;
  double grand_avg = 0;
  for(int i = 0; i < no_A; ++i) grand_avg += opt_A[i];
  for(int i = 0; i < no_B; ++i) grand_avg += opt_B[i];
  grand_avg = grand_avg / (no_A + no_B);
  for(int i = 0; i < no_A; ++i) if(opt_A[i] > grand_avg) A++;
  for(int i = 0; i < no_B; ++i) if(opt_B[i] > grand_avg) B++;
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }


// [[Rcpp::export]]
int mostlikely(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double A = opt_A[0], B = opt_B[0], pA = opt_A[no_A], pB = opt_B[no_B];
  if(no_A > 1){
    for(int i = 1; i < no_A; ++i) if(opt_A[no_A + i] > pA) A = opt_A[i];
    }
  if(no_B > 1){
    for(int i = 1; i < no_B; ++i) if(opt_B[no_B + i] > pB) B = opt_B[i];
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int leastlikely(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double out_A = opt_A[0], out_B = opt_B[0], A = opt_A[no_A], B = opt_B[no_B];
  if(no_A > 1){
    for(int i = 1; i < no_A; ++i) if(opt_A[i] < out_A) A = opt_A[no_A + i];
    }
  if(no_B > 1){
    for(int i = 1; i < no_B; ++i) if(opt_B[i] < out_B) B = opt_B[no_B + i];
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A < B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int probable(NumericVector opt_A, NumericVector opt_B, double p_limit = .5){
  int no_A = (opt_A.size() - 2)/2, no_B = (opt_B.size() - 2)/2;
  double A = 0, B = 0;
  for(int i = 0; i < no_A; ++i){
    if(opt_A[no_A + i] >= p_limit) A += opt_A[i] * opt_A[no_A + i];
    }
  for(int i = 0; i < no_B; ++i){
    if(opt_B[no_B + i] >= p_limit) B += opt_B[i] * opt_B[no_B + i];
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }

// [[Rcpp::export]]
int exppayoffelimination(NumericVector opt_A, NumericVector opt_B){
  opt_A = reduce_option(opt_A, true), opt_B = reduce_option(opt_B, true);
  int no_A = opt_A.size()/2, no_B = opt_B.size()/2;
  double exp_A, exp_B, A = 0, B = 0;
  int min_no = no_A;
  if(no_A > no_B) min_no = no_B;
  std::vector<double> outs_A, outs_B;
  for(int i = 0; i < min_no; ++i){
    exp_A = opt_A[i] * opt_A[no_A + i];
    exp_B = opt_B[i] * opt_B[no_B + i];
    if(exp_A != exp_B){
      if(exp_A > exp_B){
        A = 1;
        } else {
        B = 1;
        }
      break;
      }
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }


// [[Rcpp::export]]
int lexicographic(NumericVector opt_A, NumericVector opt_B){
  std::vector<double> out_A = sort_by_p(opt_A), out_B = sort_by_p(opt_B);
  int no_A = out_A.size(), no_B = out_B.size();
  int min_no = out_A.size();
  if(no_A > no_B) min_no = no_B;
  int A = 0, B = 0;
  for(int i = 0; i < min_no; ++i){
    if(out_A[i] != out_B[i]){
      if(out_A[i] > out_B[i]){
        A = 1;
        } else {
        B = 1;
        }
      break;
      }
    }
  int choice;
  if(A == B){
    return -1;
    } else{
    if(A > B){
      choice = 0;
      } else {
      choice = 1;
      }
    }
  return choice;
  }


// [[Rcpp::export]]
NumericMatrix toolbox(GenericVector prob){
  NumericMatrix As = prob[0];
  NumericMatrix Bs = prob[1];
  int np = As.nrow();
  NumericMatrix preds(np,13);
  for(int p = 0; p < np; p++){
    NumericVector pred(13);
    pred[0] = bayes_mean(As(p,_),Bs(p,_));
    pred[1] = natural_mean(As(p,_),Bs(p,_));
    pred[2] = equiprobable(As(p,_),Bs(p,_));
    pred[3] = maximax(As(p,_),Bs(p,_));
    pred[4] = minimax(As(p,_),Bs(p,_));
    pred[5] = minimax_regret(As(p,_),Bs(p,_));
    pred[6] = payoffelimination(As(p,_),Bs(p,_));
    pred[7] = betterthanaverage(As(p,_),Bs(p,_));
    pred[8] = mostlikely(As(p,_),Bs(p,_));
    pred[9] = leastlikely(As(p,_),Bs(p,_));
    pred[10] = probable(As(p,_),Bs(p,_));
    pred[11] = exppayoffelimination(As(p,_),Bs(p,_));
    pred[12] = lexicographic(As(p,_),Bs(p,_));
    preds(p,_) = pred;
    }
  return preds;
  }




