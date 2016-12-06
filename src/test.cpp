#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> cump_test(std::vector<double> ps){
  int i, n = ps.size();
  double cum = 0;
  std::vector<double> cump;
  for(i = 0; i < n; i++){
    cum += ps[i];
    cump.push_back(cum);
  }
  if(cump.back() > 1.0) cump[i - 1] = 1.0;
  return cump;
}


inline bool incrCompare(const std::pair<double, double>& firstElem, const std::pair<double, double>& secondElem) {
  return firstElem.first < secondElem.first;
  }

inline bool decrCompare(const std::pair<double, double>& firstElem, const std::pair<double, double>& secondElem) {
  return firstElem.first > secondElem.first;
  }

inline std::vector< std::pair<double, double> > mysort(std::vector< std::pair<double, double> > pairs, bool decreasing = true){
  if(decreasing == true)  std::sort(pairs.begin(), pairs.end(), decrCompare);
  if(decreasing == false) std::sort(pairs.begin(), pairs.end(), incrCompare);
  return pairs;
  }

inline std::vector< std::pair<double, int> > mysort2(std::vector< std::pair<double, int> > pairs, bool decreasing = true){
  if(decreasing == true)  std::sort(pairs.begin(), pairs.end(), decrCompare);
  if(decreasing == false) std::sort(pairs.begin(), pairs.end(), incrCompare);
  return pairs;
}

inline std::vector<int> sort_index(std::vector<double> o, bool decreasing = true){
  std::vector< std::pair<double, int> > outcomes;
  std::pair<double, int> outcome;
  int nout = o.size();
  for(int i = 0; i < nout; i++){
    outcome.first  = std::abs(o[i]);
    outcome.second = i;
    outcomes.push_back(outcome);
    }
  outcomes = mysort2(outcomes, decreasing);
  std::vector< std::pair<double, int> >::const_iterator it;
  std::vector<int> order;
  for (it = outcomes.begin(); it != outcomes.end(); ++it) order.push_back(it->second);
  return order;
  }


// [[Rcpp::export]]
std::vector<double> arrange(std::vector<double> opt){
  int i, n = opt.size(), outn = opt.size() / 2;
  std::pair<double, double> event;
  std::vector< std::pair<double, double> > plus, minus;
  std::vector<double> os, ps, all;
  for(i = 0; i < n / 2; i++){
    event.first  = opt[i];
    event.second = opt[i + outn];
    if(opt[i] > 0){
      plus.push_back(event);
    } else {
      minus.push_back(event);
    }
  }
  plus  = mysort(plus,true);
  minus = mysort(minus,false);
  std::vector< std::pair<double, double> >::const_iterator it;
  for (it = minus.begin(); it != minus.end(); ++it){
    os.push_back(it->first);
    ps.push_back(it->second);
  }
  for (it = plus.begin(); it != plus.end(); ++it){
    os.push_back(it->first);
    ps.push_back(it->second);
  }
  all.insert(all.end(), os.begin(), os.end());
  all.insert(all.end(), ps.begin(), ps.end());
  all.push_back(minus.size());
  return all;
}
