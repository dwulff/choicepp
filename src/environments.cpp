#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////
//
//    STOCHASTIC DOMINANCE
//
//////////////////////////////////////////////////////////////////////////////

std::vector<double> sort_opt(std::vector<double> o, std::vector<double> p){
  std::vector< std::pair<double, double> > option;
  std::pair<double, double> event;
  int nout = o.size();
  for(int i = 0; i < nout; i++){
    event.first  = o[i];
    event.second = p[i];
    option.push_back(event);
    }
  option = mysort(option, false);
  std::vector< std::pair<double, double> >::const_iterator it;
  std::vector<double> s_opt;
  for (it = option.begin(); it != option.end(); ++it) s_opt.push_back(it->first);
  for (it = option.begin(); it != option.end(); ++it) s_opt.push_back(it->second);
  return s_opt;
  }



//////////////////////////////////////////////////////////////////////////////
//' Test for 1st order stochastic dominance
//'
//' \code{stdom1} tests if one of two problems is stochastically dominant
//'   according to 1st order stochastic dominance.
//'
//' @param A numeric vector containing the outcomes and probabilities
//'   (in that order) of option A
//' @param B numeric vector containing the outcomes and probabilities
//'   (in that order) of option B
//'
//' @return A boolean indicating the status of 1st order stochastic dominance
//'
//' @export
// [[Rcpp::export]]
bool stdom1(std::vector<double> oA,
            std::vector<double> pA,
            std::vector<double> oB,
            std::vector<double> pB){
  int i, nL, nR, nA = oA.size(), nB = oB.size();
  std::vector<double> osA, osB, psA, psB;
  std::vector<double> osL, osR, psL, psR;
  std::vector<double> A = sort_opt(oA,pA);
  std::vector<double> B = sort_opt(oB,pB);
  double mxA = 0, mxB = 0;
  for(i = 0; i < nA; i++){
    if(mxA < A[i]) mxA = A[i];
    osA.push_back(A[i]);
    psA.push_back(A[i+nA]);
    }
  for(i = 0; i < nB; i++){
    if(mxB < B[i]) mxB = B[i];
    osB.push_back(B[i]);
    psB.push_back(B[i+nB]);
    }
  if(mxA < mxB){
    nL = nA;
    nR = nB;
    osL = osA;
    psL = psA;
    osR = osB;
    psR = psB;
    } else {
    nL = nB;
    nR = nA;
    osL = osB;
    psL = psB;
    osR = osA;
    psR = psA;
    }
  std::vector<double> cpsL = cump(psL);
  std::vector<double> cpsR = cump(psR);
  int j = 0;
  std::vector<bool> dom_test;
  if(osL[0] > osR[0]) return false;
  for(i = 1; i < nL; i++){
    std::vector<int> within;
    for(; j < nR; j++){
      if(osR[j] >= osL[i-1] && osR[j] <= osL[i]){
        within.push_back(j);
        } else {
        break;
        }
      }
    for(int k = 0; k < within.size(); k++){
      if(j != 0) if(cpsL[i-1] < cpsR[k]) return false;
      }
    }
  return true;
  }


//////////////////////////////////////////////////////////////////////////////
//' Test for 2nd order stochastic dominance
//'
//' \code{stdom2} tests if one of two problems is stochastically dominant
//'   according to 2nd order stochastic dominance.
//'
//' @param A numeric vector containing the outcomes and probabilities
//'   (in that order) of option A
//' @param B numeric vector containing the outcomes and probabilities
//'   (in that order) of option B
//'
//' @return A boolean indicating the status of 2nd order stochastic dominance
//'
//' @export
// [[Rcpp::export]]
bool stdom2(std::vector<double> oA,
            std::vector<double> pA,
            std::vector<double> oB,
            std::vector<double> pB,
            double d_crit = 0){
  int i, nA = oA.size(), nB = oB.size();
  double evA = 0, evB = 0, varA = 0, varB = 0;
  for(i = 0; i < nA; i++){
    evA += oA[i] * pA[i];
    varA += oA[i] * oA[i] * pA[i];
    }
  varA -= evA * evA;
  for(i = 0; i < nB; i++){
    evB += oB[i] * pB[i];
    varB += oB[i] * oB[i] * pB[i];
    }
  varB -= evB * evB;

  //std::cout << evA << '\t' << evB << '\t' << varA << '\t' << varB << '\n';

  if(evA > evB && varA <= varB) return true;
  if(evA >= evB && varA < varB) return true;
  if(evA < evB && varA >= varB) return true;
  if(evA <= evB && varA > varB) return true;

  // remove equal problems
  if(evA == evB && varA == varB) return true;

  // remove too distinct problems
  if(d_crit != 0) if((std::abs(evA-evB)/sqrt(varA+varB)) > d_crit) return true;

  //std::cout << std::abs(evA-evB)/sqrt(varA+varB) << '\n';

  return false;
  }



//////////////////////////////////////////////////////////////////////////////
//' Generate random decision environments
//'
//' \code{pgen_rnd} generates random environments of twp-option choice
//'   problems without 1st and 2nd order stochastic dominance and an
//'   ecologically
//'
//' @param ns integer vector of length three specifying the number of choice
//'   problems in the domain of Gain, Loss, and Mixed (in that order).
//' @param nA integer specifying the number of outcomes in option A
//' @param nB integer specifying the number of outcomes in option B
//' @param lower integer specifying the lower limit of the outcome range
//' @param upper integer specifying the upper limit of the outcome range
//' @param ecological bool specifies whether the outomes should be re-
//'   arranged such that high outcomes are coupled with low probabilities
//' @param stdom1_test bool specifies wehther the function controls for 1st oder
//'   stochastic dominance. If true, only non-dominant problems will be
//'   returned.
//' @param stdom1_test bool specifies wehther the function controls for 2nd oder
//'   stochastic dominance. If true, only non-dominant problems will be
//'   returned.
//' @param d_crit double specifiying the maximum effect size between the options.
//'     I.e., a d_crit value of .5 implies that options differences in
//'     expected values is at most .5 times the pooled standard deviation.
//'
//' @return A matrix of choice problems
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pgen_rnd(std::vector<int> ns,
                       int nA = 2,
                       int nB = 2,
                       double lower = -100,
                       double upper =  100,
                       bool ecological = true,
                       bool stdom1_test = true,
                       bool stdom2_test = true,
                       double d_crit = .5) {
  std::vector<double> ps;
  int nG = ns[0], nL = ns[1], nM = ns[2];
  NumericMatrix problems(nG + nL + nM, nA * 2 + nB * 2);
  int iG = 0, iL = 0, iM = 0, i = 0;
  while(iG < nG){
    std::vector<double> psA = nrnf_bound(nA);
    std::vector<double> osA = nrnf(nA,0,upper);
    std::vector<double> psB = nrnf_bound(nB);
    std::vector<double> osB = nrnf(nB,0,upper);

    if(ecological){
      std::vector<double> ps, os;
      for(int i = 0; i < nA; i++){
        ps.push_back(psA[i]);
        os.push_back(osA[i]);
        }
      for(int i = 0; i < nB; i++){
        ps.push_back(psB[i]);
        os.push_back(osB[i]);
        }
      std::vector<int> order = sort_index(os);
      std::sort(ps.begin(), ps.end());
      double psumA = 0, psumB = 0;
      for(int i = 0; i < order.size(); i++){
        int ord = order[i];
        if(ord < nA){
          psumA += ps[i];
          psA[ord] = ps[i];
          } else {
          psumB += ps[i];
          psB[ord-nA] = ps[i];
          }
        }
      for(int i = 0; i < nA; i++) psA[i] = psA[i] / psumA;
      for(int i = 0; i < nA; i++) psB[i] = psB[i] / psumB;
      }

    if(stdom1_test == true) if(nA > 1 || nB > 1) if(stdom1(osA,psA,osB,psB)) continue;
    if(stdom2_test == true) if(nA > 1 || nB > 1) if(stdom2(osA,psA,osB,psB,d_crit)) continue;

    for(int j = 0; j < nA; j++){
      problems(i,j) = osA[j];
      problems(i,j + nA) = psA[j];
      }
    ps = nrnf_bound(nB);
    for(int j = 0; j < nB; j++){
      problems(i,j + nA * 2) = osB[j];
      problems(i,j + nA * 2 + nB) = psB[j];
      }
    i++;
    iG++;
    }
  while(iL < nL){
    std::vector<double> psA = nrnf_bound(nA);
    std::vector<double> osA = nrnf(nA,lower,0);
    std::vector<double> psB = nrnf_bound(nB);
    std::vector<double> osB = nrnf(nB,lower,0);

    if(ecological){
      std::vector<double> ps, os;
      for(int i = 0; i < nA; i++){
        ps.push_back(psA[i]);
        os.push_back(osA[i]);
      }
      for(int i = 0; i < nB; i++){
        ps.push_back(psB[i]);
        os.push_back(osB[i]);
      }
      std::vector<int> order = sort_index(os);
      std::sort(ps.begin(), ps.end());
      double psumA = 0, psumB = 0;
      for(int i = 0; i < order.size(); i++){
        int ord = order[i];
        if(ord < nA){
          psumA += ps[i];
          psA[ord] = ps[i];
        } else {
          psumB += ps[i];
          psB[ord-nA] = ps[i];
        }
      }
      for(int i = 0; i < nA; i++) psA[i] = psA[i] / psumA;
      for(int i = 0; i < nA; i++) psB[i] = psB[i] / psumB;
    }


    if(stdom1_test == true) if(nA > 1 || nB > 1)  if(stdom1(osA,psA,osB,psB)) continue;
    if(stdom2_test == true) if(nA > 1 || nB > 1)  if(stdom2(osA,psA,osB,psB,d_crit)) continue;

    for(int j = 0; j < nA; j++){
      problems(i,j) = osA[j];
      problems(i,j + nA) = psA[j];
    }
    ps = nrnf_bound(nB);
    for(int j = 0; j < nB; j++){
      problems(i,j + nA * 2) = osB[j];
      problems(i,j + nA * 2 + nB) = psB[j];
    }
    i++;
    iL++;
  }
  while(iM < nM){
    std::vector<double> psA = nrnf_bound(nA);
    std::vector<double> osA;
    osA.push_back(rnf(lower,upper));
    if(osA[0] < 0) osA.push_back(rnf(0,upper));
    else osA.push_back(rnf(lower,0));
    for(int j = 2; j < nA; j++) osA.push_back(rnf(lower,upper));
    std::vector<double> psB = nrnf_bound(nB);
    std::vector<double> osB;
    osB.push_back(rnf(lower,upper));
    if(osB[0] < 0) osB.push_back(rnf(0,upper));
    else osB.push_back(rnf(lower,0));
    for(int j = 2; j < nA; j++) osB.push_back(rnf(lower,upper));

    if(ecological){
      std::vector<double> ps, os;
      for(int i = 0; i < nA; i++){
        ps.push_back(psA[i]);
        os.push_back(osA[i]);
      }
      for(int i = 0; i < nB; i++){
        ps.push_back(psB[i]);
        os.push_back(osB[i]);
      }
      std::vector<int> order = sort_index(os);
      std::sort(ps.begin(), ps.end());
      double psumA = 0, psumB = 0;
      for(int i = 0; i < order.size(); i++){
        int ord = order[i];
        if(ord < nA){
          psumA += ps[i];
          psA[ord] = ps[i];
        } else {
          psumB += ps[i];
          psB[ord-nA] = ps[i];
        }
      }
      for(int i = 0; i < nA; i++) psA[i] = psA[i] / psumA;
      for(int i = 0; i < nA; i++) psB[i] = psB[i] / psumB;
    }

    if(stdom1_test == true) if(nA > 1 || nB > 1) if(stdom1(osA,psA,osB,psB)) continue;
    if(stdom2_test == true) if(nA > 1 || nB > 1) if(stdom2(osA,psA,osB,psB,d_crit)) continue;

    for(int j = 0; j < nA; j++){
      problems(i,j) = osA[j];
      problems(i,j + nA) = psA[j];
    }
    ps = nrnf_bound(nB);
    for(int j = 0; j < nB; j++){
      problems(i,j + nA * 2) = osB[j];
      problems(i,j + nA * 2 + nB) = psB[j];
    }
    i++;
    iM++;
  }
  return problems;
  }




//////////////////////////////////////////////////////////////////////////////
//
//    ENVIRONMENT ORGANIZERS
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//' Rearrange problems
//'
//' \code{p_arrange} rearranges the matrix of problem to match the cumulative and
//'   decumulative specifications of CPT and TAX
//'
//' @param problems np x (2*nA + 2*nB) matrix containing the np problems each
//'   having nA and nB outcomes.
//' @param nA number of outcomes in option A
//'
//' @return a list of two matrices containing the outcomes and probabilities of
//'   the A and B options rearranged for CPT and TAX.
//'
//' @export
// [[Rcpp::export]]
GenericVector p_arrange(NumericMatrix problems, int nA = 0){
  int i, j, nr = problems.nrow(), nc = problems.ncol();
  double tmp;
  tmp = 0;
  if(nA == 0) {
    nA = nc / 2;
    } else  {
    nA = nA * 2;
    }
  NumericMatrix  As(nr, 1 + nA);
  NumericMatrix  Bs(nr, 1 + nc - nA);
  std::vector<double> A,B,Aar,Bar;
  GenericVector trProblems(2);
  for(i = 0; i < nr; i++){
    for(j = 0;  j < nA; j++) A.push_back(problems(i,j));
    for(j = nA; j < nc; j++) B.push_back(problems(i,j));
    Aar = arrange(A);
    Bar = arrange(B);
    for(j = 0;  j < Aar.size(); j++){
      tmp = Aar[j];
      As(i,j) = tmp;
      }
    for(j = 0;  j < Bar.size(); j++){
      tmp = Bar[j];
      Bs(i,j) = tmp;
      }
    A.clear();
    B.clear();
    }
  trProblems[0] = As;
  trProblems[1] = Bs;
  return trProblems;
  }








