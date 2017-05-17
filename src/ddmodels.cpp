#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////////////////////
//
//        DD choice rules
//
////////////////////////////////////////////////////////////////////////////////////////////////

//' Exponential choice rule
//'
//' Exponential choice rule is defined as
//' \eqn{p_ss = \frac{1}{1 + e^{\phi * (u_ll - u_ss)}}}
//'
//' @param diff difference in utility between sooner smaller
//'   and larger later options.
//' @param phi sensitivity parameter.
//'
//' @return Probability to choose sooner smaller option \eqn{p_ss}.
//'
//' @export
// [[Rcpp::export]]
double expo_rule(double diff, double phi){
  return 1 / (1 + exp(-phi * diff));
  }

//' Homothetic choice rule
//'
//' Homothetic choice rule is defined as
//' \eqn{p_ss = \frac{u_ss^\phi}{u_ss^\phi + u_ll^\phi}}
//'
//' @param u_ss utility of sooner smaller option.
//' @param u_ll utility of larger later option.
//' @param phi sensitivity parameter.
//'
//' @return Probability to choose sooner smaller option \eqn{p_ss}.
//'
//' @export
// [[Rcpp::export]]
double homo_rule(double u_ss, double u_ll, double phi){
  double a = std::pow(u_ll,phi);
  return a/(a + std::pow(u_ss,phi));
  }


////////////////////////////////////////////////////////////////////////////////////////////////
//
//        EXPONENTIAL MODEL
//
////////////////////////////////////////////////////////////////////////////////////////////////

//' Expo 1
//'
double expo_a(double outcome, double delay, double k){
  return outcome * std::pow(k, delay);
  }

//' Expo 2
//'
double expo_b(double outcome, double delay, double k){
  return outcome * exp(-k * delay);
  }

//' EXPO likelihood
//'
//' Likelihood function for various exponential discounting models.
//'
//' Function takes parameter values, the choice problem, and the actual choices
//' to compute a likelihood of the choices given the parameters.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements \eqn{k^delay} variant with exponential choice rule
//'
//' 01 implements \eqn{k^delay} variant with homothetic choice rule
//'
//' 10 implements \eqn{e^{-k*delay}} variant with exponential choice rule
//'
//' 11 implements \eqn{e^{-k*delay}} variant with homothetic choice rule
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param choices numeric vector containing the choice choice (0 = sooner smaller,
//'   1 = larger later).
//' @param epsilon numeric scalar specifying a boundary on the likelihood of a choice.
//'   Specifically, likelihoods will be bound between epsilon and 1-epsilon.
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return The negative log-liklihood XXX.
//'
//' @export
// [[Rcpp::export]]
double EXPO_lik(std::vector<double> par,
                NumericMatrix problems,
                std::vector<int> choices,
                double epsilon = .0001,
                int type = 00) {
  int n = problems.nrow();
  std::vector<double> ps_ll(n);
  double u_ss, u_ll, p_ss, ll = 0;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else if(type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_b(problems(i,2),problems(i,3),par[0]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else if(type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_b(problems(i,2),problems(i,3),par[0]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, problems[i]);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return -ll;
  }


//' EXPO deterministic choices
//'
//' Generate deterministic choice from various exponential discounting models.
//'
//' Function takes parameter values and the choice problem to generate
//' deterministic choices.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements \eqn{k^delay} variant with exponential choice rule
//'
//' 01 implements \eqn{k^delay} variant with homothetic choice rule
//'
//' 10 implements \eqn{e^{k*delay}} variant with exponential choice rule
//'
//' 11 implements \eqn{e^{k*delay}} variant with homothetic choice rule
//'
//'
//' Note that for deterministic choice 00 = 01 and 10 = 11, i.e., choice rules
//' become irrelevant.
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param type integer specifying the model variant to be used. See details.
//' @param random logical specifiyng whether ties should lead to NA (here -1) or
//'   random choice.
//'
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> EXPO_choice(std::vector<double> par,
                             NumericMatrix problems,
                             int type,
                             bool random = false) {
  int n = problems.nrow();
  std::vector<int> choices(n);
  double u_ss, u_ll;
  if(type == 00 || type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      choices[i] = detchoice(u_ll - u_ss, random);
      }
    } else if(type == 10 || type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),par[0],problems(i,1));
      u_ll = expo_b(problems(i,2),par[0],problems(i,3));
      choices[i] = detchoice(u_ss - u_ll, random);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return choices;
  }

//' EXPO random choices
//'
//' Generate random choice from various exponential discounting models.
//'
//' Function takes parameter values and the choice problem to generate
//' deterministic choices.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements \eqn{k^delay} variant with exponential choice rule
//'
//' 01 implements \eqn{k^delay} variant with homothetic choice rule
//'
//' 10 implements \eqn{e^{k*delay}} variant with exponential choice rule
//'
//' 11 implements \eqn{e^{k*delay}} variant with homothetic choice rule
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param type integer specifying the model variant to be used. See details.
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<int>  EXPO_rndchoice(std::vector<double> par,
                                 NumericMatrix problems,
                                 int type = 00) {
  int n = problems.nrow();
  std::vector<int> choices(n);
  double u_ss, u_ll, p_ss;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_b(problems(i,2),problems(i,3),par[0]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_b(problems(i,2),problems(i,3),par[0]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return choices;
  }

//' EXPO probabilities
//'
//' Generate predicted choice probabilities from exponential discounting model with
//' exponential choice rule.
//'
//' Function takes parameter values and the choice problem to generate
//' probabilistic predictions.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements \eqn{k^delay} variant with exponential choice rule
//'
//' 01 implements \eqn{k^delay} variant with homothetic choice rule
//'
//' 10 implements \eqn{e^{k*delay}} variant with exponential choice rule
//'
//' 11 implements \eqn{e^{k*delay}} variant with homothetic choice rule
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return A vector of choice probabilities.
//'
//' @export
// [[Rcpp::export]]
std::vector<double>  EXPO_prob(std::vector<double> par,
                                 NumericMatrix problems,
                                 int type = 00) {
  int n = problems.nrow();
  std::vector<double> probs(n);
  double u_ss, u_ll;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      probs[i] = expo_rule(u_ll - u_ss, par.back());
      }
    } else if(type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = expo_a(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_a(problems(i,2),problems(i,3),par[0]);
      probs[i] = homo_rule(u_ss, u_ll, par.back());
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_b(problems(i,2),problems(i,3),par[0]);
      probs[i] = expo_rule(u_ll - u_ss, par.back());
      }
    } else if(type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = expo_b(problems(i,0),problems(i,1),par[0]);
      u_ll = expo_b(problems(i,2),problems(i,3),par[0]);
      probs[i] = homo_rule(u_ss, u_ll, par.back());
      }
    } else {
    Rf_error("Non matching type value");
    }
  return probs;
  }



////////////////////////////////////////////////////////////////////////////////////////////////
//
//        HYPER
//
////////////////////////////////////////////////////////////////////////////////////////////////

//' Hyper 1
//'
double hyper_1(double outcome, double delay, double k){
  return outcome * (1.0 / (1.0 + k * delay));
  }

//' Hyper 2
//'
double hyper_2(double outcome, double delay, double k, double s){
  return outcome * std::pow(1.0 / (1.0 + k * delay),s);
  }

//' HYPER likelihood
//'
//' Likelihood function for various hyperbolic discounting models.
//'
//' Function takes parameter values, the choice problem, and the actual choices
//' to compute a likelihood of the choices given the parameters.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements the one parameter variant with exponential choice rule
//'
//' 01 implements the one parameter variant with homothetic choice rule
//'
//' 10 implements the two parameter variant with exponential choice rule
//'
//' 11 implements the two parameter variant with homothetic choice rule
//'
//'
//' The delay of the one parameter variant is modeled as
//' \deqn{1.0 / (1.0 + k * delay)}
//'
//' and that of the two parameter variant is modeled as
//' \eqn{(1.0 / (1.0 + k * delay))^s}
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param choices numeric vector containing the choice choice (0 = sooner smaller,
//'   1 = larger later).
//' @param epsilon numeric scalar specifying a boundary on the likelihood of a choice.
//'   Specifically, likelihoods will be bound between epsilon and 1-epsilon.
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return The negative log-liklihood.
//'
//' @export
// [[Rcpp::export]]
double HYPER_lik(std::vector<double> par,
                 NumericMatrix problems,
                 std::vector<int> choices,
                 double epsilon = .0001,
                 int type = 00) {
  int n = problems.nrow();
  std::vector<double> ps_ll(n);
  double u_ss, u_ll, p_ss, ll = 0;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else if(type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else if(type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll += get_lik(p_ss, choices[i]);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return -ll;
  }


//' HYPER deterministic choices
//'
//' Generate deterministic choice from various hyperbolic discounting models.
//'
//' Function takes parameter values and the choice problem to generate
//' deterministic choices.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements the one parameter variant with exponential choice rule
//'
//' 01 implements the one parameter variant with homothetic choice rule
//'
//' 10 implements the two parameter variant with exponential choice rule
//'
//' 11 implements the two parameter variant with homothetic choice rule
//'
//'
//' The delay of the one parameter variant is modeled as
//' \deqn{1.0 / (1.0 + k * delay)}
//'
//' and that of the two parameter variant is modeled as
//' \eqn{(1.0 / (1.0 + k * delay))^s}
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param type integer specifying the model variant to be used. See details.
//' @param random logical specifiyng whether ties should lead to NA (here -1) or
//'   random choice.
//'
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> HYPER_choice(std::vector<double> par,
                              NumericMatrix problems,
                              int type = 00,
                              bool random = false) {
  int n = problems.nrow();
  std::vector<int> choices(n);
  double u_ss, u_ll;
  if(type == 00 || type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      choices[i] = detchoice(u_ll - u_ss, random);
      }
    } else if(type == 10 || type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      choices[i] = detchoice(u_ll - u_ss, random);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return choices;
  }

//' HYPER random choices
//'
//' Generate random choice from various hyperbolic discounting models.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements the one parameter variant with exponential choice rule
//'
//' 01 implements the one parameter variant with homothetic choice rule
//'
//' 10 implements the two parameter variant with exponential choice rule
//'
//' 11 implements the two parameter variant with homothetic choice rule
//'
//'
//' The delay of the one parameter variant is modeled as
//' \deqn{1.0 / (1.0 + k * delay)}
//'
//' and that of the two parameter variant is modeled as
//' \eqn{(1.0 / (1.0 + k * delay))^s}
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<int>  HYPER_rndchoice(std::vector<double> par,
                                  NumericMatrix problems,
                                  int type = 00) {
  int n = problems.nrow();
  std::vector<int> choices(n);
  double u_ss, u_ll, p_ss;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      p_ss = expo_rule(u_ll - u_ss, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      p_ss = homo_rule(u_ss, u_ll, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return choices;
  }


//' HYPER probabilities
//'
//' Generate predicted choice probabilities from various hyperbolic
//' discounting models.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements the one parameter variant with exponential choice rule
//'
//' 01 implements the one parameter variant with homothetic choice rule
//'
//' 10 implements the two parameter variant with exponential choice rule
//'
//' 11 implements the two parameter variant with homothetic choice rule
//'
//'
//' The delay of the one parameter variant is modeled as
//' \deqn{1.0 / (1.0 + k * delay)}
//'
//' and that of the two parameter variant is modeled as
//' \eqn{(1.0 / (1.0 + k * delay))^s}
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return A vector of choice probabilities.
//'
//' @export
// [[Rcpp::export]]
std::vector<double>  HYPER_prob(std::vector<double> par,
                                  NumericMatrix problems,
                                  int type = 00) {
  int n = problems.nrow();
  std::vector<double> probs(n);
  double u_ss, u_ll;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      probs[i] = expo_rule(u_ll - u_ss, par.back());
      }
    } else if(type == 01){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_1(problems(i,0), problems(i,1), par[0]);
      u_ll = hyper_1(problems(i,2), problems(i,3), par[0]);
      probs[i] = homo_rule(u_ss, u_ll, par.back());
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      probs[i] = expo_rule(u_ll - u_ss, par.back());
      }
    } else if(type == 11){
    for(int i = 0; i < n; ++i){
      u_ss = hyper_2(problems(i,0), problems(i,1), par[0], par[1]);
      u_ll = hyper_2(problems(i,2), problems(i,3), par[0], par[1]);
      probs[i] = homo_rule(u_ss, u_ll, par.back());
      }
    } else {
    Rf_error("Non matching type value");
    }
  return probs;
  }

////////////////////////////////////////////////////////////////////////////////////////////////
//
//        ITCH
//
////////////////////////////////////////////////////////////////////////////////////////////////

//' ITCH 4
//'
//' @export
// [[Rcpp::export]]
double itch_4(NumericVector problem, std::vector<double> par){
  return par[0] * problem[0] + par[1] * problem[1] + par[2] * problem[2] + par[3] * problem[3];
  }

//' ITCH 5
//'
//' @export
// [[Rcpp::export]]
double itch_5(NumericVector problem, std::vector<double> par){
  return par[4] * problem[0] + par[1] * problem[1] + par[2] * problem[2] + par[3] * problem[3];
  }

//' ITCH data
//'
//' Generates matrix containing the four variables needed to compute
//' the ITCH.
//'
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//'
//' @return Four column matrix containing the difference in outcomes, the differences
//' in outcomes normalized by the mean of outcomes, the difference in delays, and the
//' difference in delays normalized by the mean of delays.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix ITCH_pgen(NumericMatrix problems){
  int n = problems.nrow();
  NumericMatrix n_problems(n,problems.ncol());
  for(int i = 0; i < n; ++i){
    NumericVector problem = problems(i,_);
    double o_d = problem[2] - problem[0];
    double d_d = problem[3] - problem[1];
    n_problems(i,0) = o_d;
    n_problems(i,1) = o_d / ((problem[0] + problem[2])/2);
    n_problems(i,2) = d_d;
    n_problems(i,3) = d_d / ((problem[1] + problem[3])/2);
    }
  return n_problems;
  }

//' ITCH deterministic choices
//'
//' Generate deterministic choice from various intertemporal choice heuristics (ITCH).
//'
//' Function takes parameter values and the choice problem to generate
//' deterministic choices.
//'
//' Functions implements different variants as a function of \code{type}:
//'
//' 00 implements the four parameter variant with exponential choice rule
//'
//' 10 implements the five parameter variant with exponential choice rule
//'
//'
//' The four and five parameter variants are defined as
//' \deqn{\beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//' and
//'
//' \deqn{\beta_5 +
//'       \beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//'
//' With o_ss, o_ll, d_ss, d_ll being the outcomes and delays of the
//' sooner smaller and larger later options, respectively, and o* and d*
//' are the arithmetic means of both outcomes and delays, respectively.
//'
//' Note that the ITCH cannot be used with a homothetic choice rule.
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param choices numeric vector containing the choice choice (0 = sooner smaller,
//'   1 = larger later).
//' @param epsilon numeric scalar specifying a boundary on the likelihood of a choice.
//'   Specifically, likelihoods will be bound between epsilon and 1-epsilon.
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return The negative log-liklihood.
//'
//' @export
// [[Rcpp::export]]
double ITCH_lik(std::vector<double> par,
                 NumericMatrix problems,
                 std::vector<int> choices,
                 double epsilon = .0001,
                 int type = 00) {
  int n = problems.nrow();
  std::vector<double> ps_ll(n);
  double u_d, p_ss, ll = 0;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_d  = itch_4(problems(i,_), par);
      p_ss = expo_rule(u_d, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll  += get_lik(p_ss, choices[i]);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_d  = itch_5(problems(i,_), par);
      p_ss = expo_rule(u_d, par.back());
      p_ss = smooth_p(p_ss, epsilon);
      ll  += get_lik(p_ss, choices[i]);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return -ll;
  }

//' ITCH deterministic choices
//'
//' Generate deterministic choice from various intertemporal choice heuristics (ITCH).
//'
//' Function takes parameter values and the choice problem to generate
//' deterministic choices.
//'
//' Functions implements different model variants as a function of \code{type}:
//'
//' 00 implements the four parameter variant with exponential choice rule
//'
//' 10 implements the five parameter variant with exponential choice rule
//'
//'
//' The four and five parameter variants are defined as
//' \deqn{\beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//' and
//'
//' \deqn{\beta_5 +
//'       \beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//'
//' With o_ss, o_ll, d_ss, d_ll being the outcomes and delays of the
//' sooner smaller and larger later options, respectively, and o* and d*
//' are the arithmetic means of both outcomes and delays, respectively.
//'
//' Note that the ITCH cannot be used with a homothetic choice rule.
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param choices numeric vector containing the choice choice (0 = sooner smaller,
//'   1 = larger later).
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> ITCH_choice(std::vector<double> par,
                NumericMatrix problems,
                int type = 00,
                bool random = false) {
  int n = problems.nrow();
  std::vector<int> choices(n);
  double u_d;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_d  = itch_4(problems(i,_), par);
      choices[i] = detchoice(u_d, random);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_d  = itch_5(problems(i,_), par);
      choices[i] = detchoice(u_d, random);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return choices;
}


//' ITCH random choices
//'
//' Generate random choice from various intertemporal choice heuristics (ITCH).
//'
//' Function takes parameter values and the choice problem to generate
//' random choices.
//'
//' Functions implements different model variants as a function of \code{type}:
//'
//' 00 implements the four parameter variant with exponential choice rule
//'
//' 10 implements the five parameter variant with exponential choice rule
//'
//'
//' The four and five parameter variants are defined as
//' \deqn{\beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//' and
//'
//' \deqn{\beta_5 +
//'       \beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//'
//' With o_ss, o_ll, d_ss, d_ll being the outcomes and delays of the
//' sooner smaller and larger later options, respectively, and o* and d*
//' are the arithmetic means of both outcomes and delays, respectively.
//'
//' Note that the ITCH cannot be used with a homothetic choice rule.
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param choices numeric vector containing the choice choice (0 = sooner smaller,
//'   1 = larger later).
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<int> ITCH_rndchoice(std::vector<double> par,
                             NumericMatrix problems,
                             int type = 00) {
  int n = problems.nrow();
  std::vector<int> choices(n);
  double p_ss, u_d;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_d  = itch_4(problems(i,_), par);
      p_ss = expo_rule(u_d, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_d  = itch_5(problems(i,_), par);
      p_ss = expo_rule(u_d, par.back());
      choices[i] = rndchoice(p_ss);
      }
    } else {
    Rf_error("Non matching type value");
    }
  return choices;
  }


//' ITCH probabilities
//'
//' Generate predicted choice probabilities from variants of the intertemporal
//' choice heuristic (ITCH).
//'
//' Function takes parameter values and the choice problem to generate
//' random choices.
//'
//' Functions implements different model variants as a function of \code{type}:
//'
//' 00 implements the four parameter variant with exponential choice rule
//'
//' 10 implements the five parameter variant with exponential choice rule
//'
//'
//' The four and five parameter variants are defined as
//' \deqn{\beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//' and
//'
//' \deqn{\beta_5 +
//'       \beta_0 * (o_ss - o_ll) + \beta_1 * \frac{o_ss - o_ll}{o*} +
//'       \beta_2 * (d_ss - d_ll) + \beta_3 * \frac{d_ss - d_ll}{d*}}
//'
//' With o_ss, o_ll, d_ss, d_ll being the outcomes and delays of the
//' sooner smaller and larger later options, respectively, and o* and d*
//' are the arithmetic means of both outcomes and delays, respectively.
//'
//' Note that the ITCH cannot be used with a homothetic choice rule.
//'
//' @param par numeric vector specifying the discounting factor
//'   and the choice sensitivity (in that order).
//' @param problems a numeric matrix with four columns containing for each decision
//'   the outcome and delay of the sooner smaller option, the outcome and delay
//'   of the larger later option (in that order).
//' @param choices numeric vector containing the choice choice (0 = sooner smaller,
//'   1 = larger later).
//' @param type integer specifying the model variant to be used. See details.
//'
//' @return A vector of choices.
//'
//' @export
// [[Rcpp::export]]
std::vector<double> ITCH_prob(std::vector<double> par,
                                NumericMatrix problems,
                                int type = 00) {
  int n = problems.nrow();
  std::vector<double> probs(n);
  double u_d;
  if(type == 00){
    for(int i = 0; i < n; ++i){
      u_d  = itch_4(problems(i,_), par);
      probs[i] = expo_rule(u_d, par.back());
      }
    } else if(type == 10){
    for(int i = 0; i < n; ++i){
      u_d  = itch_5(problems(i,_), par);
      probs[i] = expo_rule(u_d, par.back());
      }
    } else {
    Rf_error("Non matching type value");
    }
  return probs;
  }

