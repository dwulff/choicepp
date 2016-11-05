#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////
//
//    TAX HELPERS
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//' Value function of TAX
//'
//' \code{u}, the TAX value function, transforms the magnitute of
//'   an outcome into utility space.
//'
//' @param o numeric specifying the magnitude of the to be transformed
//'   outcome
//' @param alpha numeric specifying the exponent of the utility function for
//'   the gain domain.
//'
//' @return a utility
//'
//' @export
// [[Rcpp::export]]// [[Rcpp::export]]
double u(double o,  double alpha){
  if(o < 0) return -1 * pow(std::abs(o), alpha);
  return pow(std::abs(o), alpha);
  }


//////////////////////////////////////////////////////////////////////////////
//' Weighting function of TAX
//'
//' \code{t} is the weighting function proposed by Tverksy & Kahneman
//'   (1992). It transforms probabilities into decision weights.
//'
//' @param p numeric specifying the to be transformed probability.
//' @param gamma numeric specifying the exponential factor
//'   (sensitivity) of the probability weighting function.
//'
//' @return a decision weight
//'
//' @references Birnbaum, M. H. (1999). Testing critical properties of
//'   decision making on the Internet. Psychological Science, 10(5), 399-407.
//'
//' @export
// [[Rcpp::export]]
double t(double p, double gamma){
  return pow(std::abs(p), gamma);
  }


//////////////////////////////////////////////////////////////////////////////
//' Transfer function of TAX
//'
//' \code{t} is the weighting function proposed by Tverksy & Kahneman
//'   (1992). It transforms probabilities into decision weights.
//'
//' @param pi numeric specifying the to be transformed probability at i.
//' @param pk numeric specifying the to be transformed probability at k.
//' @param delta numeric specifying the transfer factor.
//' @param gamma numeric specifying the exponential factor (sensitivity) of
//'   the probability weighting function. See \link{t}.
//'
//' @return a transfer weight
//'
//' @references Birnbaum, M. H. (1999). Testing critical properties of
//'   decision making on the Internet. Psychological Science, 10(5), 399-407.
//'
//' @export
// [[Rcpp::export]]
double w(double pi, double pk, double delta, double gamma, double n){
  if(delta < 0) return (delta*t(pk,gamma))/(n+1);
  else return (delta*t(pk,gamma))/(n+1);
  return 0;
  }

//////////////////////////////////////////////////////////////////////////////
//' Utility function of TAX
//'
//' \code{t} is the weighting function proposed by Tverksy & Kahneman
//'   (1992). It transforms probabilities into decision weights.
//'
//' @param opt numeric vector specifying the outcomes and probabilities of an
//'   option. The function expects a length of (number of outcomes * 2 + 1),
//'   that is ordered according to max_loss, min_loss, max_gain, min_gain, and
//'   whose last entry indicates the number of loss outcomes.
//' @param par numeric vector specifying the parameters of the CPT model. See
//'   \link{u}, \link{t}, and \link{w}.
//' @param type integer specifying the parameterization of TAX.
//'   \code{type = 1} is the one parameter TAX with just a \code{gamma}
//'   parameter. \code{type = 2} is the two parameter TAX with a \code{gamma}
//'   and a \code{alpha} parameter. \code{type = 2} is the two parameter TAX
//'   with a \code{gamma} and a \code{delta} parameter. \code{type = 2} is the
//'   three parameter TAX with a \code{gamma}, a \code{alpha}, and a
//'   \code{delta} parameter.
//'
//' @return a utility
//'
//' @references Birnbaum, M. H. (1999). Testing critical properties of
//'   decision making on the Internet. Psychological Science, 10(5), 399-407.
//'
//' @export
// [[Rcpp::export]]
double utility_tax(NumericVector opt, std::vector<double> par, int type){
  int no = (opt.size() - 1)/2;
  double uxi, nom = 0, den = 0;
  if(type == 0){
    for(int i = 0; i < no; i++){
      uxi = opt[i];
      nom += uxi * t(opt[i + no],par[0]);
      den += t(opt[i + no],par[0]);
      for(int k = 0; k < i; k++){
        nom += (uxi - opt[k]) * w(opt[i + no],opt[k + no],1,par[0],double(no));
        }
      }
    }
  if(type == 1){
    for(int i = 0; i < no; i++){
      uxi = u(opt[i],par[0]);
      nom += uxi * t(opt[i + no],par[1]);
      den += t(opt[i + no],par[1]);
      for(int k = 0; k < i; k++){
        nom += (uxi - u(opt[k],par[0])) * w(opt[i + no],opt[k + no],1,par[1],double(no));
        }
      }
    }
  if(type == 2){
    for(int i = 0; i < no; i++){
      uxi = opt[i];
      nom += uxi * t(opt[i + no],par[0]);
      den += t(opt[i + no],par[0]);
      for(int k = 0; k < i; k++){
        nom += (uxi - opt[k]) * w(opt[i + no],opt[k + no],par[1],par[0],double(no));
        }
      }
    }
  if(type == 3){
    for(int i = 0; i < no; i++){
      uxi = u(opt[i],par[0]);
      nom += uxi * t(opt[i + no],par[1]);
      den += t(opt[i + no],par[1]);
      for(int k = 0; k < i; k++){
        nom += (uxi - u(opt[k],par[0])) * w(opt[i + no],opt[k + no],par[2],par[1],double(no));
        }
      }
    }
  return nom / den;
  }


//////////////////////////////////////////////////////////////////////////////
//
//    TAX FUNS
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//' TAX-based choice probabilities
//'
//' \code{tax_prob} computes for a set of choice problems and a set of
//'   parameters the respective probabilities of choosing option A.
//'
//' @param par numeric vector specifying the parameters of the TAX model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of tax model. See
//'   \link{utility_tax}.
//'
//' @return a vector of choice probabilities
//'
//' @export
// [[Rcpp::export]]
std::vector<double> tax_prob(std::vector<double> par,
                             GenericVector problems,
                             int type = 0){
  NumericMatrix As = problems[0], Bs = problems[1];
  int i, n = As.nrow();
  double utA, utB;
  NumericVector A, B;
  std::vector<double> probs;
  for(i = 0; i < n; i++){
    A = As(i,_);
    B = Bs(i,_);
    utA = utility_tax(A, par, type);
    utB = utility_tax(B, par, type);
    probs.push_back(choice_rule(utA,utB,par.back()));
    }
  return probs;
  }

//////////////////////////////////////////////////////////////////////////////
//' TAX-based likelihood
//'
//' \code{tax_prob} computes for a set of choice problems and a set of
//'   parameters the combined (negative) log-likelihood of a set of choices.
//'
//' @param par numeric vector specifying the parameters of the TAX model.
//' @param problems list as generated by \link{p_arrange}.
//' @param choices numeric vector of 0s and 1s indicating the choices where
//'   0 = A and 0 = B.
//' @param type integer specifying the parameterization of the tax model. See
//'   \link{utility_tax}.
//'
//' @return a vector of choice probabilities
//'
//' @export
// [[Rcpp::export]]
double tax_lik(std::vector<double> par,
               GenericVector problems,
               std::vector<int> choices,
               int type = 000,
               double limit = .0001){
  int choice;
  double pA, llik = 0;
  std::vector<double> probs = tax_prob(par,problems,type);
  for(int i = 0; i < probs.size(); i++){
    //cout << llik << "\n";
    choice = choices[i];
    pA  = probs[i];
    if(pA < limit) pA = limit;
    if((1-pA) < limit) pA = 1 - limit;
    if(choice == 0){
      llik += log(pA);
    } else {
      llik += log(1 - pA);
    }
  }
  return -llik;
}

//////////////////////////////////////////////////////////////////////////////
//' TAX-based deterministic choices
//'
//' \code{tax_choice} produces deterministic choices for a set of choice
//'   problems and parameters.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of the tax model
//'   function. See \link{utility_tax}.
//'
//' @return a vector of choices
//'
//' @export
// [[Rcpp::export]]
std::vector<int> tax_choice(std::vector<double> par,
                            GenericVector problems,
                            int type = 000){
  NumericMatrix As = problems[0], Bs = problems[1];
  int i, n = As.nrow();
  double utA, utB, r;
  NumericVector A, B;
  std::vector<int> choices;
  for(i = 0; i < n; i++){
    //cout << llik << "\n";
    A = As(i,_);
    B = Bs(i,_);
    utA = utility_tax(A, par, type);
    utB = utility_tax(B, par, type);
    if( utA != utB){
      if(utA > utB){
        choices.push_back(0);
      } else {
        choices.push_back(1);
      }
    } else {
      r = double(std::rand()) / RAND_MAX;
      if(r < .5){
        choices.push_back(0);
      } else {
        choices.push_back(1);
      }
    }
  }
  return choices;
}


//////////////////////////////////////////////////////////////////////////////
//' TAX-based probabilistic choices
//'
//' \code{tax_rndchoice} produces probabilistic choices for a set of choice
//'   problems and parameters.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of the tax model
//'   function. See \link{utility_tax}.
//'
//' @return a vector of choices
//'
//' @export
// [[Rcpp::export]]
std::vector<int> tax_rndchoice(std::vector<double> par,
                               GenericVector problems,
                               int type = 000){
  double r, pA;
  std::vector<double> probs = tax_prob(par,problems,type);
  std::vector<int> choices;
  for(int i = 0; i < probs.size(); i++){
    pA  = probs[i];
    r = double(std::rand()) / RAND_MAX;
    //std::cout << pA << "\n";
    if(r < pA){
      choices.push_back(0);
    } else {
      choices.push_back(1);
    }
  }
  return choices;
}
