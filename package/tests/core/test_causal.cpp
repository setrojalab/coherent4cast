#include <librecurcast/causal_target.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
int main() {
  constexpr std::size_t n=8, bins=256;
  const double lambda=.9, censor=.3, e=.4;
  const std::vector<double> times{.1,.4,.7,1.,.1,.4,.7,1.};
  const std::vector<int> events{0,1,0,1,0,1,0,1};
  std::vector<recurcast::CTRecord> rows;
  for(std::size_t i=4;i<n;++i) rows.push_back({i,0,0,0,times[i],events[i]});
  recurcast::CausalContinuousTarget engine(1,bins,1,1,1.,{1},std::vector<std::size_t>(n,0),rows,{censor},{e});
  auto result=engine.fit(std::vector<double>(bins,lambda),0,1e-10);
  const double risk=-std::expm1(-lambda),rmst=risk/lambda;
  assert(std::abs(result.estimate[0]-risk)<1e-11);
  assert(std::abs(result.estimate[1]-rmst)<1e-11);
  double maximum=0;
  for(std::size_t i=0;i<n;++i){
    double a=0,b=0;
    if(i>=4){
      const double t=times[i],death=events[i]==0?1.:0.,inc=std::exp(-lambda*(1-t));
      const double first=std::expm1(censor*t)/censor;
      const double second=std::expm1((lambda+censor)*t)/(lambda+censor);
      a=(death*std::exp(censor*t)*inc-lambda*std::exp(-lambda)*second)/e;
      b=(-death*std::exp(censor*t)*(1-inc)/lambda+first-std::exp(-lambda)*second)/e;
    }
    maximum=std::max(maximum,std::abs(result.eif[2*i]-a));
    maximum=std::max(maximum,std::abs(result.eif[2*i+1]-b));
  }
  assert(maximum<2e-9);
  rows.clear();for(std::size_t i=0;i<n;++i)rows.push_back({i,0,0,0,times[i],events[i]});
  recurcast::ContinuousTarget old(1,bins,1,1,1.,{1},std::vector<std::size_t>(n,0),rows,{censor});
  recurcast::CausalContinuousTarget reduced(1,bins,1,1,1.,{1},std::vector<std::size_t>(n,0),rows,{censor},{1.});
  auto a=old.fit(std::vector<double>(bins,lambda),20,1e-8);
  auto b=reduced.fit(std::vector<double>(bins,lambda),20,1e-8);
  assert(a.estimate==b.estimate && a.eif==b.eif && a.rates==b.rates);
  bool rejected=false;
  try {recurcast::CausalContinuousTarget bad(1,bins,1,1,1.,{1},std::vector<std::size_t>(n,0),rows,{censor},{0.});}
  catch(const std::invalid_argument&){rejected=true;}
  assert(rejected);
  std::cout<<"Causal analytic EIF error="<<maximum<<"; probability-one reduction exact; positivity guard passed\n";
}
