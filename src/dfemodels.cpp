#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////////////////////
//
//        Natural mean
//
////////////////////////////////////////////////////////////////////////////////////////////////


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
double nm_long(std::vector<double> opt, std::vector<double> out, bool extr = false) {
  std::vector<double> s0, s1;
  int n = opt.size();
  for(int i = 0; i < n; ++i){
    if(opt[i] == 0){
      s0.push_back(out[i]);
      } else {
      s1.push_back(out[i]);
      }
    }
  double m0 = 0, m1 = 0;
  std::vector<double>::const_iterator it;
  for(it = s0.begin(); it != s0.end(); ++it) m0 += *it;
  for(it = s1.begin(); it != s1.end(); ++it) m1 += *it;
  m0 /= s0.size();
  m1 /= s1.size();
  //std::cout << m0 << '\t' << m1 << '\n';
  int choice;
  if(extr == true) return m1 / (m1 + m0);
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

////////////////////////////////////////////////////////////////////////////////////////////////
//
//        Natural mean + recency
//
////////////////////////////////////////////////////////////////////////////////////////////////

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

//' Compute natural mean with recency
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


////////////////////////////////////////////////////////////////////////////////////////////////
//
//        Round-wise
//
////////////////////////////////////////////////////////////////////////////////////////////////

//' Compute round wise
//'
//' @export
// [[Rcpp::export]]
double rw(std::vector<double> opt, std::vector<double> out, bool extr = false) {
  int ind = 1, n = opt.size();
  std::vector<double> ms0;
  std::vector<double> ms1;
  double m0, m1;
  if(opt[0] == 0){
    m0 = out[0];
    } else {
    m1 = out[0];
    }
  for(int i = 1; i < n; ++i){
    if(opt[i - 1] == opt[i]){ //make sure this is not done for n+1
      ind++;
      if(opt[i] == 0){
        m0 += out[i];
        } else {
        m1 += out[i];
        }
      } else {
      if(opt[i-1] == 0){
        ms0.push_back(m0 / ind);
        m1 = out[i];
        } else {
        ms1.push_back(m1 / ind);
        m0 = out[i];
        }
      ind = 1;
      }
    }
  if(opt.back() == 0){
    ms0.push_back(m0 / ind);
    } else {
    ms1.push_back(m1 / ind);
    }


  //for(int i = 0; i < ms0.size(); ++i) std::cout << ms0[i] << '\n';
  //for(int i = 0; i < ms1.size(); ++i) std::cout << ms1[i] << '\n';

  ind = 1;
  std::vector<int> comps;
  if(ms0.size() > ms1.size()){
    int n = ms1.size();
    for(int i = 0; i < n; ++i){
      if(ms0[i] > ms1[i]) comps.push_back(0);
      if(ms0[i] < ms1[i]) comps.push_back(1);
      if(ms0[i + 1] > ms1[i]) comps.push_back(0);
      if(ms0[i + 1] < ms1[i]) comps.push_back(1);
      }
    } else if(ms0.size() < ms1.size()){
    int n = ms0.size();
    for(int i = 0; i < n; ++i){
      if(ms0[i] > ms1[i]) comps.push_back(0);
      if(ms0[i] < ms1[i]) comps.push_back(1);
      if(ms0[i] > ms1[i + 1]) comps.push_back(0);
      if(ms0[i] < ms1[i + 1]) comps.push_back(1);
      }
    } else {
    int n = ms0.size();
    for(int i = 0; i < n; ++i){
      if(ms0[i] > ms1[i]) comps.push_back(0);
      if(ms0[i] < ms1[i]) comps.push_back(1);
      if(i != (n-1)){
        if(ms0[i] > ms1[i + 1]) comps.push_back(0);
        if(ms0[i] < ms1[i + 1]) comps.push_back(1);
        if(ms0[i + 1] > ms1[i]) comps.push_back(0);
        if(ms0[i + 1] < ms1[i]) comps.push_back(1);
        }
      }
    }

  double p = 0;
  n = comps.size();
  if(n == 0){
    return -1;
    } else {
    for(int i = 0; i < n; i++) p += comps[i];
    p /= n;
    }
  if(extr == true) return p;
  if(p == .5) return -1;
  double choice;
  if(p < .5){
    choice = 0;
    } else {
    choice = 1;
    }
  return choice;
  }


////////////////////////////////////////////////////////////////////////////////////////////////
//
//        Selective Integration
//
////////////////////////////////////////////////////////////////////////////////////////////////

//' Compute round wise
//'
//' @export
// [[Rcpp::export]]
double si(std::vector<double> opt, std::vector<double> out, double w, bool extr = false, bool verbose = false) {
  int ind = 1, n = opt.size();
  std::vector<double> ms0, ms1;
  std::vector<double> os0, os1;
  std::vector<int> is0, is1;
  double m0, m1;

  // set first value
  if(opt[0] == 0){
    m0 = out[0];
    os0.push_back(out[0]);
  } else {
    m1 = out[0];
    os1.push_back(out[0]);
  }

  // iterate through outcomes and compute means for option 0 and option 1
  for(int i = 1; i < n; ++i){
    if(opt[i - 1] == opt[i]){ //make sure this is not done for n+1
      ind++;
      if(opt[i] == 0){
        m0 += out[i];
        os0.push_back(out[i]);
      } else {
        m1 += out[i];
        os1.push_back(out[i]);
      }
    } else {
      if(opt[i-1] == 0){
        ms0.push_back(m0 / ind);
        is0.push_back(ind);
        m1 = out[i];
        os1.push_back(out[i]);
      } else {
        ms1.push_back(m1 / ind);
        is1.push_back(ind);
        m0 = out[i];
        os0.push_back(out[i]);
      }
      ind = 1;
    }
  }

  // handle last outcome
  if(opt.back() == 0){
    ms0.push_back(m0 / ind);
    is0.push_back(ind);
  } else {
    ms1.push_back(m1 / ind);
    is1.push_back(ind);
  }

  // make is cumulative
  int ni0 = is0.size();
  int ni1 = is1.size();
  for(int i = 1; i < ni0; ++i) is0[i] += is0[i-1];
  for(int i = 1; i < ni1; ++i) is1[i] += is1[i-1];

  // length of values
  int n0 = ms0.size();
  int n1 = ms1.size();
  int no0 = os0.size();
  int no1 = os1.size();

  // determine comparisons
  std::vector<double> w0(n0, 0.0);
  std::vector<double> w1(n1, 0.0);
  if(ms0.size() > ms1.size()){
    for(int i = 0; i < n1; ++i){
      if(ms0[i] > ms1[i]) w1[i] = w * std::abs(ms0[i] - ms1[i]);
      if(ms0[i] < ms1[i]) w0[i] = w * std::abs(ms0[i] - ms1[i]);
      if(ms0[i + 1] > ms1[i]) w1[i] = w * std::abs(ms0[i + 1] - ms1[i]);
      if(ms0[i + 1] < ms1[i]) w0[i + 1] = w * std::abs(ms0[i + 1] - ms1[i]);
    }
  } else if(ms0.size() < ms1.size()){
    for(int i = 0; i < n0; ++i){
      if(ms0[i] > ms1[i]) w1[i] = w * std::abs(ms0[i] - ms1[i]);
      if(ms0[i] < ms1[i]) w0[i] = w * std::abs(ms0[i] - ms1[i]);
      if(ms0[i] > ms1[i + 1]) w1[i + 1] = w * std::abs(ms0[i] - ms1[i + 1]);
      if(ms0[i] < ms1[i + 1]) w0[i] = w * std::abs(ms0[i] - ms1[i + 1]);
    }
  } else {
    for(int i = 0; i < n0; ++i){
      if(ms0[i] > ms1[i]) w1[i] = w * std::abs(ms0[i] - ms1[i]);
      if(ms0[i] < ms1[i]) w0[i] = w * std::abs(ms0[i] - ms1[i]);
      if(i != (n0-1)){
        if(ms0[i] > ms1[i + 1]) w1[i + 1] = w * std::abs(ms0[i] - ms1[i + 1]);
        if(ms0[i] < ms1[i + 1]) w0[i] = w * std::abs(ms0[i] - ms1[i + 1]);
        if(ms0[i + 1] > ms1[i]) w1[i] = w * std::abs(ms0[i + 1] - ms1[i]);
        if(ms0[i + 1] < ms1[i]) w0[i + 1] = w * std::abs(ms0[i + 1] - ms1[i]);
      }
    }
  }

  // compute utilities
  double u0 = 0, u1 = 0;
  int j;
  j = 0;
  for(int i = 0; i < no0; ++i){
    if(verbose == true) std::cout  << '0' << '\t' << os0[i] << '\t' << w0[j] << '\n';
    u0 += os0[i] - w0[j];
    if(i == (is0[j]-1)) j++;
  }
  j = 0;
  for(int i = 0; i < no1; ++i){
    if(verbose == true) std::cout  << '1' << '\t' << os1[i] << '\t' << w1[j] << '\n';
    u1 += os1[i] - w1[j];
    if(i == (is1[j]-1)) j++;
  }
  u0 /= no0;
  u1 /= no1;

  // choose
  double choice;
  if(extr == true) return u1 / (u1 + u0);
  if(u0 != u1){
    if(u0 > u1){
      choice = 0;
    } else {
      choice = 1;
    }
  } else {
    choice = -1;
  }
  return choice;
}





//' Compute round wise
//'
//' @export
// [[Rcpp::export]]
double si2(std::vector<double> opt, std::vector<double> out, double w, bool extr = false, bool verbose = false) {
  int ind = 1, n = opt.size();
  std::vector<double> ms0, ms1;
  std::vector<double> os0, os1;
  std::vector<int> is0, is1;
  double m0, m1;

  // set first value
  if(opt[0] == 0){
    m0 = out[0];
    os0.push_back(out[0]);
  } else {
    m1 = out[0];
    os1.push_back(out[0]);
  }

  // iterate through outcomes and compute means for option 0 and option 1
  for(int i = 1; i < n; ++i){
    if(opt[i - 1] == opt[i]){ //make sure this is not done for n+1
      ind++;
      if(opt[i] == 0){
        m0 += out[i];
        os0.push_back(out[i]);
      } else {
        m1 += out[i];
        os1.push_back(out[i]);
      }
    } else {
      if(opt[i-1] == 0){
        ms0.push_back(std::abs(m0 / ind));
        is0.push_back(ind);
        m1 = out[i];
        os1.push_back(out[i]);
      } else {
        ms1.push_back(std::abs(m1 / ind));
        is1.push_back(ind);
        m0 = out[i];
        os0.push_back(out[i]);
      }
      ind = 1;
    }
  }

  // handle last outcome
  if(opt.back() == 0){
    ms0.push_back(std::abs(m0 / ind));
    is0.push_back(ind);
  } else {
    ms1.push_back(std::abs(m1 / ind));
    is1.push_back(ind);
  }

  // make is cumulative
  int ni0 = is0.size();
  int ni1 = is1.size();
  for(int i = 1; i < ni0; ++i) is0[i] += is0[i-1];
  for(int i = 1; i < ni1; ++i) is1[i] += is1[i-1];

  // length of values
  int n0 = ms0.size();
  int n1 = ms1.size();
  int no0 = os0.size();
  int no1 = os1.size();

  // determine comparisons
  std::vector<double> w0(n0, 1.0);
  std::vector<double> w1(n1, 1.0);
  if(ms0.size() > ms1.size()){
    for(int i = 0; i < n1; ++i){
      if(ms0[i] > ms1[i]) w1[i] -= w;
      if(ms0[i] < ms1[i]) w0[i] -= w;
      if(ms0[i + 1] > ms1[i]) w1[i] -= w;
      if(ms0[i + 1] < ms1[i]) w0[i + 1] -= w;
    }
  } else if(ms0.size() < ms1.size()){
    for(int i = 0; i < n0; ++i){
      if(ms0[i] > ms1[i]) w1[i] -= w;
      if(ms0[i] < ms1[i]) w0[i] -= w;
      if(ms0[i] > ms1[i + 1]) w1[i + 1] -= w;
      if(ms0[i] < ms1[i + 1]) w0[i] -= w;
    }
  } else {
    for(int i = 0; i < n0; ++i){
      if(ms0[i] > ms1[i]) w1[i] -= w;
      if(ms0[i] < ms1[i]) w0[i] -= w;
      if(i != (n0-1)){
        if(ms0[i] > ms1[i + 1]) w1[i + 1] -= w;
        if(ms0[i] < ms1[i + 1]) w0[i] -= w;
        if(ms0[i + 1] > ms1[i]) w1[i] -= w;
        if(ms0[i + 1] < ms1[i]) w0[i + 1] -= w;
      }
    }
  }

    // compute utilities
    double u0 = 0, u1 = 0;
  int j;
  j = 0;
  for(int i = 0; i < no0; ++i){
    if(verbose == true) std::cout  << '0' << '\t' << os0[i] << '\t' << w0[j] << '\n';
    u0 += os0[i] * w0[j];
    if(i == (is0[j]-1)) j++;
  }
  j = 0;
  for(int i = 0; i < no1; ++i){
    if(verbose == true) std::cout  << '1' << '\t' << os1[i] << '\t' << w1[j] << '\n';
    u1 += os1[i] * w1[j];
    if(i == (is1[j]-1)) j++;
  }
  u0 /= no0;
  u1 /= no1;

  // choose
  double choice;
  if(extr == true) return u1 / (u1 + u0);
  if(u0 != u1){
    if(u0 > u1){
      choice = 0;
    } else {
      choice = 1;
    }
  } else {
    choice = -1;
  }
  return choice;
}


//' Compute round wise
//'
//' @export
// [[Rcpp::export]]
double si3(std::vector<double> opt, std::vector<double> out, double w, bool extr = false, bool verbose = false) {
  int ind = 1, n = opt.size();
  std::vector<double> ms0, ms1;
  std::vector<double> os0, os1;
  std::vector<int> is0, is1;
  double m0, m1;

  // set first value
  if(opt[0] == 0){
    m0 = out[0];
    os0.push_back(out[0]);
  } else {
    m1 = out[0];
    os1.push_back(out[0]);
  }

  // iterate through outcomes and compute means for option 0 and option 1
  for(int i = 1; i < n; ++i){
    if(opt[i - 1] == opt[i]){ //make sure this is not done for n+1
      ind++;
      if(opt[i] == 0){
        m0 += out[i];
        os0.push_back(out[i]);
      } else {
        m1 += out[i];
        os1.push_back(out[i]);
      }
    } else {
      if(opt[i-1] == 0){
        ms0.push_back(std::abs(m0 / ind));
        is0.push_back(ind);
        m1 = out[i];
        os1.push_back(out[i]);
      } else {
        ms1.push_back(std::abs(m1 / ind));
        is1.push_back(ind);
        m0 = out[i];
        os0.push_back(out[i]);
      }
      ind = 1;
    }
  }

  // handle last outcome
  if(opt.back() == 0){
    ms0.push_back(std::abs(m0 / ind));
    is0.push_back(ind);
  } else {
    ms1.push_back(std::abs(m1 / ind));
    is1.push_back(ind);
  }

  // make is cumulative
  int ni0 = is0.size();
  int ni1 = is1.size();
  for(int i = 1; i < ni0; ++i) is0[i] += is0[i-1];
  for(int i = 1; i < ni1; ++i) is1[i] += is1[i-1];

  // length of values
  int n0 = ms0.size();
  int n1 = ms1.size();
  int no0 = os0.size();
  int no1 = os1.size();

  // determine comparisons
  std::vector<double> w0(n0, 1.0);
  std::vector<double> w1(n1, 1.0);
  if(ms0.size() > ms1.size()){
    for(int i = 0; i < n1; ++i){
      if(ms0[i] < ms1[i]) w1[i] += w;
      if(ms0[i] > ms1[i]) w0[i] += w;
      if(ms0[i + 1] < ms1[i]) w1[i] += w;
      if(ms0[i + 1] > ms1[i]) w0[i + 1] += w;
    }
  } else if(ms0.size() < ms1.size()){
    for(int i = 0; i < n0; ++i){
      if(ms0[i] < ms1[i]) w1[i] += w;
      if(ms0[i] > ms1[i]) w0[i] += w;
      if(ms0[i] < ms1[i + 1]) w1[i + 1] += w;
      if(ms0[i] > ms1[i + 1]) w0[i] += w;
    }
  } else {
    for(int i = 0; i < n0; ++i){
      if(ms0[i] < ms1[i]) w1[i] += w;
      if(ms0[i] > ms1[i]) w0[i] += w;
      if(i != (n0-1)){
        if(ms0[i] < ms1[i + 1]) w1[i + 1] += w;
        if(ms0[i] > ms1[i + 1]) w0[i] += w;
        if(ms0[i + 1] < ms1[i]) w1[i] += w;
        if(ms0[i + 1] > ms1[i]) w0[i + 1] += w;
      }
    }
  }

  // compute utilities
  double u0 = 0, u1 = 0;
  int j;
  j = 0;
  for(int i = 0; i < no0; ++i){
    if(verbose == true) std::cout  << '0' << '\t' << os0[i] << '\t' << w0[j] << '\n';
    u0 += (1/(i + 1.0)) * (os0[i] - u0) * w0[j];
    if(i == (is0[j]-1)) j++;
  }
  j = 0;
  for(int i = 0; i < no1; ++i){
    if(verbose == true) std::cout  << '1' << '\t' << os1[i] << '\t' << w1[j] << '\n';
    u1 += (1/(i + 1.0)) * (os1[i] - u1) * w1[j];
    if(i == (is1[j]-1)) j++;
  }

  if(verbose == true) std::cout << u0 << '\t' << u1 << '\n';

  // choose
  double choice;
  if(extr == true) return u1 / (u1 + u0);
  if(u0 != u1){
    if(u0 > u1){
      choice = 0;
    } else {
      choice = 1;
    }
  } else {
    choice = -1;
  }
  return choice;
}







