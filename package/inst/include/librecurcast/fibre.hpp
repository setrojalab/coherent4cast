#ifndef LIBRECURCAST_FIBRE_HPP
#define LIBRECURCAST_FIBRE_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace recurcast {
struct FibreEvent { double lower, upper; int mark; bool terminal; };
struct FibrePath { std::vector<int> order; std::vector<double> time; double log_weight; int draw; };
struct FibreResult {
  std::vector<FibrePath> paths;
  std::string method, reason;
  std::size_t permutations=0, attempted=0, rejected=0;
  double volume=0, volume_mcse=0, importance_ess=0;
};

// Integration over unlabelled event histories requires a unique matching of
// same-mark occurrences to overlapping observation windows. Earliest deadline
// among available windows is the canonical matching. It avoids counting the
// same ordered history once for every admissible record assignment.
inline bool fibre_priority(const FibreEvent& a, std::size_t ia,
                           const FibreEvent& b, std::size_t ib) {
  if (a.upper != b.upper) return a.upper < b.upper;
  if (a.lower != b.lower) return a.lower < b.lower;
  return ia < ib;
}

inline FibreResult coarsening_fibre(const std::vector<FibreEvent>& events,
    double followup, const std::vector<double>& breaks,
    const std::vector<double>& nodes, const std::vector<double>& weights,
    std::size_t max_perms=120, std::size_t draws=2000,
    std::uint64_t seed=20261122, std::size_t max_paths=100000) {
  if (!std::isfinite(followup) || followup<=0 || nodes.empty() ||
      nodes.size()!=weights.size() || !max_perms || draws<2 || !max_paths)
    throw std::invalid_argument("Invalid fibre follow-up, quadrature or budget");
  const std::size_t m=events.size();
  if(m>64) throw std::invalid_argument("At most 64 recorded events per observation fibre are supported");
  std::size_t deaths=0;
  for (const auto& e:events) {
    if(!std::isfinite(e.lower)||!std::isfinite(e.upper)||e.lower<0||e.upper<e.lower||e.upper>followup||e.upper<=0||e.mark<0)
      throw std::invalid_argument("Event windows must lie in (0, followup]");
    deaths += e.terminal;
  }
  if(deaths>1) throw std::invalid_argument("Only one terminal event is allowed per subject");
  for(std::size_t j=0;j<nodes.size();++j)
    if(!(nodes[j]>0&&nodes[j]<1&&weights[j]>0&&std::isfinite(weights[j])))
      throw std::invalid_argument("Quadrature nodes must be interior with positive weights");
  FibreResult result;
  if(m==0) {result.paths.push_back({{}, {}, 0.0, 0});result.method="quadrature";
    result.reason="no recorded events";result.permutations=1;result.volume=1;return result;}
  std::vector<bool> used(m,false);
  auto local_upper = [&](std::size_t e, const std::vector<bool>& taken) {
    double up=events[e].upper;
    for(std::size_t j=0;j<m;++j) if(!taken[j]&&j!=e&&
      events[j].mark==events[e].mark&&events[j].terminal==events[e].terminal&&
      fibre_priority(events[j],j,events[e],e)) up=std::min(up,events[j].lower);
    return up;
  };
  std::vector<std::vector<int>> orders;
  std::vector<int> order;
  std::function<void(double)> enumerate = [&](double earliest) {
    if(orders.size()>max_perms)return;
    if(order.size()==m){orders.push_back(order);return;}
    for(std::size_t e=0;e<m;++e) {
      if(used[e] || (events[e].terminal && order.size()+1!=m))continue;
      const double upper=local_upper(e,used),lower=std::max(earliest,events[e].lower);
      const bool point=events[e].upper==events[e].lower;
      if(point ? !(events[e].upper>earliest && events[e].upper<=upper) : !(upper>lower))continue;
      used[e]=true;order.push_back(static_cast<int>(e));enumerate(lower);
      order.pop_back();used[e]=false;if(orders.size()>max_perms)return;
    }
  };
  enumerate(0);
  result.permutations=orders.size();
  if(orders.empty())throw std::invalid_argument("The coarsening fibre is empty: no legal strictly ordered event history");
  bool sampled=orders.size()>max_perms;
  result.reason=sampled?"admissible-order enumeration exceeded max_perms":"all canonical admissible orders enumerated";
  std::vector<double> split_points=breaks;
  for(const auto&e:events){split_points.push_back(e.lower);split_points.push_back(e.upper);}
  std::sort(split_points.begin(),split_points.end());
  split_points.erase(std::unique(split_points.begin(),split_points.end()),split_points.end());
  if(!sampled) {
    bool exhausted=false;
    for(const auto& ord:orders) {
      std::vector<double> lo(m),up(m),time(m);std::fill(used.begin(),used.end(),false);
      for(std::size_t j=0;j<m;++j){const auto e=static_cast<std::size_t>(ord[j]);lo[j]=events[e].lower;up[j]=local_upper(e,used);used[e]=true;}
      // Future upper bounds truncate earlier integration limits. This is an
      // ordered-box integral; the same-day recurrence/death case is a triangle.
      for(std::size_t j=m-1;j>0;--j)up[j-1]=std::min(up[j-1],up[j]);
      std::function<void(std::size_t,double,double)> integrate = [&](std::size_t j,double previous,double logw) {
        if(exhausted)return;
        if(j==m){if(result.paths.size()>=max_paths){exhausted=true;return;}
          result.paths.push_back({ord,time,logw,0});return;}
        const auto&e=events[static_cast<std::size_t>(ord[j])];
        if(e.lower==e.upper){if(e.upper>previous&&e.upper<=up[j]){time[j]=e.upper;integrate(j+1,e.upper,logw);}return;}
        const double a=std::max(lo[j],previous),b=up[j];if(!(b>a))return;
        std::vector<double> cuts{a};for(double v:split_points)if(v>a&&v<b)cuts.push_back(v);cuts.push_back(b);
        for(std::size_t k=1;k<cuts.size();++k) {
          const double width=cuts[k]-cuts[k-1];
          for(std::size_t q=0;q<nodes.size();++q){time[j]=cuts[k-1]+width*nodes[q];
            integrate(j+1,time[j],logw+std::log(width)+std::log(weights[q]));if(exhausted)return;}
        }
      };
      integrate(0,0,0);if(exhausted)break;
    }
    if(exhausted){sampled=true;result.reason="deterministic quadrature exceeded max_paths";result.paths.clear();}
    else if(result.paths.empty())throw std::invalid_argument("No positive-mass legal path remains after window restrictions");
  }
  if(sampled) {
    result.method="importance_sampling";result.attempted=draws;
    if(draws>max_paths)throw std::invalid_argument("Importance draws exceed max_paths; increase the path budget explicitly");
    std::mt19937_64 rng(seed);
    // Use raw engine bits rather than implementation-dependent distributions.
    auto uniform=[&](){return (static_cast<double>(rng()>>11)+0.5)*0x1.0p-53;};
    for(std::size_t draw=0;draw<draws;++draw) {
      std::fill(used.begin(),used.end(),false);order.clear();std::vector<double> time;
      double previous=0,log_inverse_proposal=0;bool failed=false;
      for(std::size_t j=0;j<m;++j) {
        std::vector<int> candidate;std::vector<double> lower,upper;
        for(std::size_t e=0;e<m;++e) {
          if(used[e]||(events[e].terminal&&j+1!=m))continue;
          double a=std::max(previous,events[e].lower),b=local_upper(e,used);
          // Any next event must still leave room before all unassigned upper
          // endpoints. This reduces, but does not conceal, rejected draws.
          for(std::size_t k=0;k<m;++k)if(!used[k]&&k!=e)b=std::min(b,events[k].upper);
          const bool point=events[e].lower==events[e].upper;
          if(point ? !(events[e].upper>previous&&events[e].upper<=b) : !(b>a))continue;
          candidate.push_back(static_cast<int>(e));lower.push_back(a);upper.push_back(b);
        }
        if(candidate.empty()){failed=true;break;}
        std::size_t k=std::min(candidate.size()-1,static_cast<std::size_t>(uniform()*static_cast<double>(candidate.size())));
        const std::size_t e=static_cast<std::size_t>(candidate[k]);
        log_inverse_proposal+=std::log(static_cast<double>(candidate.size()));
        double t=events[e].upper;
        if(events[e].upper!=events[e].lower){const double width=upper[k]-lower[k];t=lower[k]+width*uniform();log_inverse_proposal+=std::log(width);}
        used[e]=true;order.push_back(static_cast<int>(e));time.push_back(t);previous=t;
      }
      if(failed){++result.rejected;continue;}
      result.paths.push_back({order,time,log_inverse_proposal-std::log(static_cast<double>(draws)),static_cast<int>(draw+1)});
    }
    if(result.paths.empty())throw std::runtime_error("No accepted importance paths; increase draws or improve observation windows");
  } else result.method="quadrature";
  double sumsq=0;
  for(const auto&p:result.paths){const double w=std::exp(p.log_weight);result.volume+=w;sumsq+=w*w;}
  if(result.method=="importance_sampling") {
    const double M=static_cast<double>(draws);
    result.volume_mcse=std::sqrt(std::max(0.0,M/(M-1)*(sumsq-result.volume*result.volume/M)));
    result.importance_ess=sumsq>0?result.volume*result.volume/sumsq:0;
  }
  return result;
}
} // namespace recurcast
#endif
