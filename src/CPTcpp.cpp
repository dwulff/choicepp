#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////
//
//    VALUE FUNS
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//' Standard value function
//'
//' \code{v}, being the standard value function, transforms the magnitute of
//'   an outcome into utility space
//'
//' @param o numeric specifying the magnitude of the to be transformed
//'   outcome
//' @param alpha numeric specifying the exponent of the utility function for
//'   the gain domain.
//' @param beta numeric specifying the exponent of the utility function for
//'   the loss domain.
//' @param lambda numeric specifying the slope of the utility function for
//'   the loss domain relative to the gain domain.
//'
//' @return a utility
//'
//' @export
// [[Rcpp::export]]
double v(double o,  double alpha, double beta, double lambda){
  if(o < 0) return -1 * lambda * pow(std::abs(o), beta);
  return pow(std::abs(o), alpha);
  }


//////////////////////////////////////////////////////////////////////////////
//' Wrapper for value function
//'
//' \code{v_wrapper} is a wrapper for different parameterizations of the value
//'   function.
//'
//' @details
//' The type argument controls the type of value function implemented.
//'
//' xx0 implements only \code{alpha}
//'
//' xx1 implements \code{alpha} and \code{lambda}
//'
//' xx2 implements \code{alpha} and \code{beta}
//'
//' xx3 implements \code{alpha}, \code{beta}, and \code{lambda}
//'
//' @param o numeric specifying the to be transformed outcome.
//' @param type integer specifying the parameterization of the value function.
//'   See Details.
//'
//' @return a utility
//'
//' @export
// [[Rcpp::export]]
double v_wrapper(double o, std::vector<double> par, int type){
  int third = type - int(std::floor(double(type) / 10.) * 10.);
  if(third == 0)      return v(o, par[0], par[0], 1);
  else if(third == 1) return v(o, par[0], par[0], par[1]);
  else if(third == 2) return v(o, par[0], par[1], 1);
  else if(third == 3) return v(o, par[0], par[1], par[2]);
  return 0;
  }


//////////////////////////////////////////////////////////////////////////////
//
//    WEIGHTING FUNS
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//' Weighting function according to Tversky & Kahneman (1992)
//'
//' \code{w_tk} is the weighting function proposed by Tverksy & Kahneman
//'   (1992). It transforms probabilities into decision weights.
//'
//' @param p numeric specifying the to be transformed probability.
//' @param gamma_l,gamma_g numeric specifying the exponential factor
//'   (sensitivity) of the weighting function for the loss (l) and gain
//'   (g) domain.
//'
//' @return a decision weight
//'
//' @references Tversky, A., & Kahneman, D. (1992). Advances in prospect
//'   theory: Cumulative representation of uncertainty. Journal of Risk and
//'   uncertainty, 5(4), 297-323.
//'
//' @export
// [[Rcpp::export]]
double w_tk(double p, double o, double gamma_l, double gamma_g){
  if(o < 0){
    double nom = pow(p, gamma_l);
    double denom = pow(1.0 - p, gamma_l);
    return  nom / pow((nom + denom), 1/gamma_l);
    }
  if(o > 0){
    double nom = pow(p, gamma_g);
    double denom = pow(1.0 - p, gamma_g);
    return  nom / pow((nom + denom), 1/gamma_g);
    }
  return 0;
  }


//////////////////////////////////////////////////////////////////////////////
//' Weighting function according to Goldstein & Einhorn (1987)
//'
//' \code{w_ge} is the weighting function proposed by Goldstein & Einhorn
//'   (1987). It transforms probabilities into decision weights.
//'
//' @param p numeric specifying the to be transformed probability.
//' @param delta_l,delta_g numeric specifying the multiplicative factor
//'   (elevation) of the weighting function for the loss (l) and gain
//'   (g) domain.
//' @param gamma_l,gamma_g numeric specifying the exponential factor
//'   (sensitivity) of the weighting function for the loss (l) and gain
//'   (g) domain.
//'
//' @return a decision weight
//'
//' @references Goldstein, W. M. and Einhorn, H. J. (1987). Expression Theory
//'   and the Preference Reversal Phenomena. Psychological Review 94, 236—254.
//'
//' @export
// [[Rcpp::export]]
double w_ge(double p, double o, double delta_l, double delta_g, double gamma_l, double gamma_g){
  if(o < 0){
    double nom = delta_l * pow(p, gamma_l);
    double denom = pow(1.0 - p, gamma_l);
    return  nom / (nom + denom);
    }
  if(o > 0){
    double nom = delta_g * pow(p, gamma_g);
    double denom = pow(1.0 - p, gamma_g);
    return  nom / (nom + denom);
    }
  return 0.0;
  }


//////////////////////////////////////////////////////////////////////////////
//' Weighting function according to Prelec (1998)
//'
//' \code{w_p} is the weighting function proposed by Tverksy & Kahneman
//'   (1992). It transforms probabilities into decision weights.
//'
//' @param p numeric specifying the to be transformed probability.
//' @param gamma_l,gamma_g numeric specifying the exponential factor
//'   (sensitivity) of the weighting function for the loss (l) and gain
//'   (g) domain.
//'
//' @return a decision weight
//'
//' @references Prelec, D. (1998). The probability weighting function.
//'   Econometrica 66, 497–527.
//'
//' @export
// [[Rcpp::export]]
double w_p(double p, double o, double delta_l, double gamma_l, double delta_g, double gamma_g){
  if(o < 0)      return  exp(-delta_l * pow((-log(p)), gamma_l));
  else if(o > 0) return  exp(-delta_g * pow((-log(p)), gamma_g));
  return 0;
  }

//////////////////////////////////////////////////////////////////////////////
//' Wrapper for weighting function
//'
//' \code{w_wrapper} is a wrapper for different parameterizations of the
//'   weighting function.
//'
//' @details
//' The \code{type} argument controls the type of weighting function implemented.
//'
//' 00x implements the weighting function of Tversky & Kahneman
//' (1992) with one \code{gamma}.
//'
//' 01x implements the weighting function of Tversky & Kahneman
//' (1992) with one \code{gamma} for losses and one \code{gamma} for gains.
//'
//' 10x implements the weighting function of Goldstein & Einhorn
//' (1987) with one \code{delta} and one \code{gamma}.
//'
//' 11x implements the weighting function of Goldstein & Einhorn
//' (1987) with one \code{delta}, one \code{gamma} for losses, and one
//' \code{gamma} for gains.
//'
//' 12x implements the weighting function of Goldstein & Einhorn
//' (1987) with one \code{delta} for losses, one \code{delta} for gains,
//' and one \code{gamma}.
//'
//' 13x implements the weighting function of Goldstein & Einhorn
//' (1987) with one \code{delta} for losses, one \code{delta} for gains,
//' and one \code{gamma} for losses and one \code{gamma} for gains.
//'
//' 20x implements the "one"-parameter weighting function of Prelec
//' (1998) with one \code{gamma}.
//'
//' 21x implements the "one"-parameter weighting function of Prelec
//' (1998) with one \code{gamma} for losses and one \code{gamma} for gains.
//'
//' 22x implements the "two"-parameter weighting function of Prelec
//' (1998) with one \code{delta} and one \code{gamma}.
//'
//' 23x implements the "two"-parameter weighting function of Prelec
//' (1998) with one \code{delta}, one \code{gamma} losses, and one
//' \code{gamma} for gains.
//'
//' 24x implements the "two"-parameter weighting function of Prelec
//' (1998) with one \code{delta} for losses, one \code{delta} gains, and
//' one \code{gamma}.
//'
//' 25x implements the "two"-parameter weighting function of Prelec
//' (1998) with one \code{delta} for losses, one \code{delta} gains, one
//' \code{gamma} for losses and one \code{gamma} for gains.
//'
//' @param p numeric specifying the to be transformed probability.
//' @param o numeric specifying the outcome associated with the probability.
//' @param type integer specifying the parameterization of the weighting
//'   function. See Details.
//'
//' @return a decision weight
//'
//' @export
// [[Rcpp::export]]
double w_wrapper(double p, double o, std::vector<double> par, int type){
  int first = int(std::floor(double(type) / 100.));
  int second = int(std::floor(double(type - first*100) / 10.));
  int third = type - (first*100+second*10);
  if(first == 0){
    if(second == 0){
      if(third == 0)      return w_tk(p, o, par[1], par[1]);
      else if(third == 1) return w_tk(p, o, par[2], par[2]);
      else if(third == 2) return w_tk(p, o, par[2], par[2]);
      else if(third == 3) return w_tk(p, o, par[3], par[3]);
      }
    else if(second == 1){
      if(third == 0)      return w_tk(p, o, par[1], par[2]);
      else if(third == 1) return w_tk(p, o, par[2], par[3]);
      else if(third == 2) return w_tk(p, o, par[2], par[3]);
      else if(third == 3) return w_tk(p, o, par[3], par[4]);
      }
    }
  if(first == 1){
    if(second == 0){
      if(third == 0)      return w_ge(p, o, par[1], par[1], par[2], par[2]);
      else if(third == 1) return w_ge(p, o, par[2], par[2], par[3], par[3]);
      else if(third == 2) return w_ge(p, o, par[2], par[2], par[3], par[3]);
      else if(third == 3) return w_ge(p, o, par[3], par[3], par[4], par[4]);
      }
    else if(second == 1){
      if(third == 0)      return w_ge(p, o, par[1], par[1], par[2], par[3]);
      else if(third == 1) return w_ge(p, o, par[2], par[2], par[3], par[4]);
      else if(third == 2) return w_ge(p, o, par[2], par[2], par[3], par[4]);
      else if(third == 3) return w_ge(p, o, par[3], par[3], par[4], par[5]);
      }
    else if(second == 1){
      if(third == 0)      return w_ge(p, o, par[1], par[2], par[3], par[3]);
      else if(third == 1) return w_ge(p, o, par[2], par[3], par[4], par[4]);
      else if(third == 2) return w_ge(p, o, par[2], par[3], par[4], par[4]);
      else if(third == 3) return w_ge(p, o, par[3], par[4], par[5], par[5]);
      }
    else if(second == 1){
      if(third == 0)      return w_ge(p, o, par[1], par[2], par[3], par[4]);
      else if(third == 1) return w_ge(p, o, par[2], par[3], par[4], par[5]);
      else if(third == 2) return w_ge(p, o, par[2], par[3], par[4], par[5]);
      else if(third == 3) return w_ge(p, o, par[3], par[4], par[5], par[6]);
      }
    }
  else if(first == 2){
    if(second == 0){
      if(third == 0)      return w_p(p, o,      1,      1, par[1], par[1]);
      else if(third == 1) return w_p(p, o,      1,      1, par[2], par[2]);
      else if(third == 2) return w_p(p, o,      1,      1, par[2], par[2]);
      else if(third == 3) return w_p(p, o,      1,      1, par[3], par[3]);
      }
    else if(second == 1){
      if(third == 0)      return w_p(p, o,      1,      1, par[1], par[2]);
      else if(third == 1) return w_p(p, o,      1,      1, par[2], par[3]);
      else if(third == 2) return w_p(p, o,      1,      1, par[2], par[3]);
      else if(third == 3) return w_p(p, o,      1,      1, par[3], par[4]);
      }
    else if(second == 2){
      if(third == 0)      return w_p(p, o, par[1], par[1], par[2], par[2]);
      else if(third == 1) return w_p(p, o, par[2], par[2], par[3], par[3]);
      else if(third == 2) return w_p(p, o, par[2], par[2], par[3], par[3]);
      else if(third == 3) return w_p(p, o, par[3], par[3], par[4], par[4]);
      }
    else if(second == 3){
      if(third == 0)      return w_p(p, o, par[1], par[1], par[2], par[3]);
      else if(third == 1) return w_p(p, o, par[2], par[2], par[3], par[4]);
      else if(third == 2) return w_p(p, o, par[2], par[2], par[3], par[4]);
      else if(third == 3) return w_p(p, o, par[3], par[3], par[4], par[5]);
      }
    else if(second == 4){
      if(third == 0)      return w_p(p, o, par[1], par[2], par[3], par[3]);
      else if(third == 1) return w_p(p, o, par[2], par[3], par[4], par[4]);
      else if(third == 2) return w_p(p, o, par[2], par[3], par[4], par[4]);
      else if(third == 3) return w_p(p, o, par[3], par[4], par[5], par[5]);
      }
    else if(second == 5){
      if(third == 0)      return w_p(p, o, par[1], par[2], par[3], par[4]);
      else if(third == 1) return w_p(p, o, par[2], par[3], par[4], par[5]);
      else if(third == 2) return w_p(p, o, par[2], par[3], par[4], par[5]);
      else if(third == 3) return w_p(p, o, par[3], par[4], par[5], par[6]);
      }
    }
  return 0;
  }


//////////////////////////////////////////////////////////////////////////////
//' Calculate the utility of an option
//'
//' \code{utility} calculates the utility of an option based on a specific
//'   type of prospect theory. For details see \link{v_wrapper} and
//'   \link{w_wrapper}.
//'
//'
//' @param opt numeric vector specifying the outcomes and probabilities of an
//'   option. The function expects a length of (number of outcomes * 2 + 1),
//'   that is ordered according to max_loss, min_loss, max_gain, min_gain, and
//'   whose last entry indicates the number of loss outcomes.
//' @param par numeric vector specifying the parameters of the CPT model. See
//'   \link{v_wrapper} and \link{v_wrapper}.
//' @param type integer specifying the parameterization of the value and
//'   weighting function. See \link{v_wrapper} and \link{v_wrapper}.
//'
//' @return a utility.
//'
//' @export
// [[Rcpp::export]]
double utility(NumericVector opt, std::vector<double> par, int type){
  int i, n = opt.size(), no = (opt.size() - 1)/2;
  int nneg = opt[n - 1];
  int npos = no - nneg;
  double u, w, nw, ut = 0;
  std::vector<double> us, cmp, ps;
  if(nneg > 0){
    for(i = 0; i < nneg; i++){
      u = v_wrapper(opt[i],par,type);
      us.push_back(u);
      ps.push_back(opt[i + no]);
      }
    cmp = cump(ps);
    w   = w_wrapper(cmp[0],opt[0], par, type);
    ut += w * us[0];
    for(i = 1; i < nneg; i++){
      nw = w_wrapper(cmp[i],opt[i], par, type);
      ut += us[i] * (nw - w);
      w = nw;
      }
    us.clear();
    ps.clear();
    }
  if(npos > 0){
    for(i = 0; i < npos; i++){
      u = v_wrapper(opt[i + nneg],par,type);
      us.push_back(u);
      ps.push_back(opt[i + nneg + no]);
      }
    cmp = cump(ps);
    w   = w_wrapper(cmp[0],opt[nneg], par, type);
    ut += w * us[0];
    for(i = 1; i < npos; i++){
      nw = w_wrapper(cmp[i],opt[i + nneg], par, type);
      ut += us[i] * (nw - w);
      w = nw;
      }
    }
  return ut;
  }


//////////////////////////////////////////////////////////////////////////////
//
//    CPT FUNS
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//' CPT-based choice probabilities
//'
//' \code{cpt_prob} computes for a set of choice problems and a set of
//'   parameters the respective probabilities of choosing option A.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of the value and
//'   weighting function. See \link{v_wrapper} and \link{w_wrapper}.
//'
//' @return a vector of choice probabilities
//'
//' @export
// [[Rcpp::export]]
std::vector<double> cpt_prob(std::vector<double> par,
                             GenericVector problems,
                             int type = 000){
  NumericMatrix As = problems[0], Bs = problems[1];
  int i, n = As.nrow();
  double utA, utB;
  NumericVector A, B;
  std::vector<double> probs;
  for(i = 0; i < n; i++){
    A = As(i,_);
    B = Bs(i,_);
    utA = utility(A, par, type);
    utB = utility(B, par, type);
    probs.push_back(choice_rule(utA,utB,par.back()));
    }
  return probs;
  }



//////////////////////////////////////////////////////////////////////////////
//' CPT-based likelihood
//'
//' \code{cpt_prob} computes for a set of choice problems and a set of
//'   parameters the combined (negative) log-likelihood of a set of choices.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems list as generated by \link{p_arrange}.
//' @param choices numeric vector of 0s and 1s indicating the choices where
//'   0 = A and 0 = B.
//' @param type integer specifying the parameterization of the value and
//'   weighting function. See \link{v_wrapper} and \link{w_wrapper}.
//'
//' @return a vector of choice probabilities
//'
//' @export
// [[Rcpp::export]]
double cpt_lik(std::vector<double> par,
               GenericVector problems,
               std::vector<int> choices,
               int type = 000,
               double limit = .0001){
  int choice;
  double pA, llik = 0;
  std::vector<double> probs = cpt_prob(par,problems,type);
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
//' CPT-based deterministic choices
//'
//' \code{cpt_choice} produces deterministic choices for a set of choice
//'   problems and parameters.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of the value and
//'   weighting function. See \link{v_wrapper} and \link{w_wrapper}.
//'
//' @return a vector of choices
//'
//' @export
// [[Rcpp::export]]
std::vector<int> cpt_choice(std::vector<double> par,
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
    utA = utility(A, par, type);
    utB = utility(B, par, type);
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
//' CPT-based probabilistic choices
//'
//' \code{cpt_rndchoice} produces probabilistic choices for a set of choice
//'   problems and parameters.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of the value and
//'   weighting function. See \link{v_wrapper} and \link{w_wrapper}.
//'
//' @return a vector of choices
//'
//' @export
// [[Rcpp::export]]
std::vector<int> cpt_rndchoice(std::vector<double> par,
               GenericVector problems,
               int type = 000){
  double r, pA;
  std::vector<double> probs = cpt_prob(par,problems,type);
  std::vector<int> choices;
  for(int i = 0; i < probs.size(); i++){
    pA  = probs[i];
    r = double(std::rand()) / RAND_MAX;
    if(r < pA){
      choices.push_back(0);
      } else {
      choices.push_back(1);
      }
    }
  return choices;
  }


//////////////////////////////////////////////////////////////////////////////
//' CPT-based probabilistic choices plus certainty effect
//'
//' \code{cpt_rndchoice} produces probabilistic choices for a set of choice
//'   problems and parameters.
//'
//' @param par numeric vector specifying the parameters of the CPT model.
//' @param problems a list as generated by \link{p_arrange}.
//' @param type integer specifying the parameterization of the value and
//'   weighting function. See \link{v_wrapper} and \link{w_wrapper}.
//'
//' @return a vector of choices
//'
//' @export
// [[Rcpp::export]]
std::vector<int> cpt_rndchoice_cert(std::vector<double> par,
                               GenericVector problems,
                               int type = 000){
  double r, pA;
  NumericMatrix As = problems[0], Bs = problems[1];
  int i, n = As.nrow();
  double utA, utB;
  NumericVector A, B;
  std::vector<int> choices;
  for(i = 0; i < n; i++){
    A = As(i,_);
    B = Bs(i,_);
    utA = utility(A, par, type);
    utB = utility(B, par, type);
    int nA = A.size(), nB = B.size();
    bool Acert = false, Bcert = false;
    int npA = (nA/2) - 1;
    int npB = (nB/2) - 1;
    for(int i = 0; i < npA; ++i) if(A[i + npA] == 1.0) Acert = true;
    for(int i = 0; i < npB; ++i) if(B[i + npB] == 1.0) Bcert = true;
    if(Acert) utA += std::abs(((utA + utB) / 2.)) * .5;
    if(Bcert) utB += std::abs(((utA + utB) / 2.)) * .5;
    pA = choice_rule(utA,utB,par.back());
    r = double(std::rand()) / RAND_MAX;
    if(r < pA){
      choices.push_back(0);
      } else {
      choices.push_back(1);
      }
    }
  return choices;
  }

