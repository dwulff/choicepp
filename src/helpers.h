#ifndef __UTILITIES__
#define __UTILITIES__

inline bool cmp(double a, double b){
  return std::abs(a) > std::abs(b);
  }

inline double rnf(int a = 0, int b = 1){
  double r = double(std::rand()) / RAND_MAX;
  if(a == 0 && b == 1) return r;
  return r * (b - a) + a;
}

inline std::vector<int> sq(int len, int start = 0){
  int i;
  std::vector<int> sequence;
  for(i = start; i < len; i++) sequence.push_back(i);
  return sequence;
}


inline std::vector<double> nrnf(int n, int a = 0, int b = 1, bool norm = true){
  int i;
  double r;
  std::vector<double> rs, rsn;
  for(i = 0; i < n; i++){
    r = double(std::rand()) / RAND_MAX;
    rs.push_back(r * (b - a) + a);
  }
  if(norm){
    double sm = 0;
    for(i = 0; i < n; i++){
      sm += rs[i];
    }
    for(i = 0; i < n; i++){
      rsn.push_back(rs[i] / sm);
    }
    return rsn;
  }
  return rs;
}


inline std::vector<double> cump(std::vector<double> ps){
  int i, n = ps.size();
  double cum = 0;
  std::vector<double> cump;
  for(i = 0; i < n; i++){
    cum += ps[i];
    cump.push_back(cum);
  }
  if(cump.back() > 1) cump[i - 1] = 1.;
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

inline std::vector<double> arrange(std::vector<double> opt){
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

//////////////////////////////////////////////////////////////////////////////
//
//    CHOICE RULE
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////' Exponential choice rule
////'
////' \code{choice_rule} calculates the probability of choosing A using an ex-
////'   ponential choice rule.
////'
////' @param utA numeric specifying the utility of option A
////' @param utB numeric specifying the utility of option B
////' @param phi numeric specifying the choice sensitivity
////'
////' @return a choice probability
////'
////' @export
//// [[Rcpp::export]]
inline double choice_rule(double utA, double utB, double phi){
  return 1 / (1 + exp(phi * (utB - utA)));
  }


#endif // __UTILITIES__
