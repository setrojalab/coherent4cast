#ifndef LIBRECURCAST_SPARSE_UNIFORMIZE_HPP
#define LIBRECURCAST_SPARSE_UNIFORMIZE_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace recurcast {
// Off-diagonal generator entries only. Holding rates are the row sums.
// CSR column indices and offsets are zero based; duplicate transitions are legal.
struct CSRGenerator {
 std::size_t n=0;
 std::vector<std::size_t> offsets, columns;
 std::vector<double> rates;
};
struct SparsePropagationResult {
 std::vector<double> probability, reward, reward_error_bound;
 double error_bound=0,roundoff_bound=0;
 std::size_t terms=0,chunks=0,workspace_bytes=0;
};
inline SparsePropagationResult sparse_uniformize(const std::vector<double>& initial,
 const CSRGenerator& q,double duration,const std::vector<double>& rewards,
 std::size_t nr,double tolerance=1e-12,std::size_t max_terms=10000,
 std::size_t max_chunks=100000) {
 const std::size_t n=q.n;
 if(n==0 || initial.size()!=n || q.offsets.size()!=n+1 ||
    q.columns.size()!=q.rates.size() || q.offsets.front()!=0 ||
    q.offsets.back()!=q.rates.size() || rewards.size()!=n*nr ||
    !std::isfinite(duration) || duration<0 || !std::isfinite(tolerance) ||
    tolerance<=0 || tolerance>=1 || max_terms<2 || max_chunks<1)
  throw std::invalid_argument("sparse uniformization: invalid dimensions, duration, tolerance or budget");
 for(std::size_t i=0;i<n;++i)if(q.offsets[i]>q.offsets[i+1] || q.offsets[i+1]>q.rates.size())
  throw std::invalid_argument("sparse uniformization: CSR offsets must be monotone and in range");
 std::vector<long double> exits(n,0),p(initial.begin(),initial.end()),max_reward(nr,0);
 long double mass=0,nu=0;
 for(std::size_t i=0;i<n;++i) {
  if(!std::isfinite(initial[i]) || initial[i]<0 || q.offsets[i]>q.offsets[i+1])
   throw std::invalid_argument("sparse uniformization: invalid probability or CSR offset");
  mass+=initial[i];
  for(std::size_t k=q.offsets[i];k<q.offsets[i+1];++k) {
   if(q.columns[k]>=n || q.columns[k]==i || !std::isfinite(q.rates[k]) || q.rates[k]<0)
    throw std::invalid_argument("sparse uniformization: only finite nonnegative off-diagonal rates are accepted");
   exits[i]+=q.rates[k];
  }
  nu=std::max(nu,exits[i]);
  for(std::size_t j=0;j<nr;++j) {
   if(!std::isfinite(rewards[i*nr+j]))throw std::invalid_argument("sparse uniformization: nonfinite reward");
   max_reward[j]=std::max(max_reward[j],std::abs(static_cast<long double>(rewards[i*nr+j])));
  }
 }
 if(!(mass>0) || mass>1+1e-10L || !std::isfinite(nu))
  throw std::invalid_argument("sparse uniformization: initial mass must be in (0,1], rates must be finite");
 SparsePropagationResult out;out.reward.assign(nr,0);out.reward_error_bound.assign(nr,0);
 out.workspace_bytes=(5*n+5*nr)*sizeof(long double)+(n+2*nr)*sizeof(double);
 if(duration==0 || nu==0) {
  out.probability=initial;
  if(duration>0)for(std::size_t j=0;j<nr;++j) {
   long double total=0;for(std::size_t i=0;i<n;++i)total+=static_cast<long double>(initial[i])*rewards[i*nr+j]*duration;
   if(!std::isfinite(total) || std::abs(total)>std::numeric_limits<double>::max())
    throw std::runtime_error("sparse uniformization: integrated reward overflow");
   out.reward[j]=static_cast<double>(total);
  }
  return out;
 }
 // Chunking avoids dense transition matrices and exp(-nu*t) underflow. The
 // cost scales with nu*t; an explicit budget prevents an unbounded computation.
 const long double total_x=nu*duration;
 const long double wanted=std::max(1.0L,std::ceil(total_x/16.0L));
 if(!std::isfinite(wanted) || wanted>static_cast<long double>(max_chunks))
  throw std::runtime_error("sparse uniformization: max_chunks exceeded; raise budget or rescale the model");
 out.chunks=static_cast<std::size_t>(wanted);
 const long double h=static_cast<long double>(duration)/static_cast<long double>(out.chunks);
 const long double x=nu*h,local_tol=tolerance/(2.0L*static_cast<long double>(out.chunks));
 std::vector<long double> v(n),next(n),answer(n),cumulative(nr),reward_answer(nr),reward_total(nr,0),reward_error(nr,0);
 long double probability_error=0;
 for(std::size_t chunk=0;chunk<out.chunks;++chunk) {
  v=p;std::fill(answer.begin(),answer.end(),0);std::fill(cumulative.begin(),cumulative.end(),0);
  std::fill(reward_answer.begin(),reward_answer.end(),0);
  long double weight=std::exp(-x),tail=1,previous_tail=1;
  for(std::size_t i=0;i<n;++i)answer[i]=weight*v[i];
  bool converged=false;
  for(std::size_t k=1;k<=max_terms;++k) {
   // Integral identity: E[int_0^h r(X_t)dt] = sum_m Pr(N=m)
   // sum_{l=0}^{m-1} p P^l r / nu. It uses only row-vector products.
   for(std::size_t j=0;j<nr;++j) {
    long double pr=0;for(std::size_t i=0;i<n;++i)pr+=v[i]*rewards[i*nr+j];
    cumulative[j]+=pr/nu;
   }
   std::fill(next.begin(),next.end(),0);
   for(std::size_t i=0;i<n;++i) {
    next[i]+=v[i]*(1-exits[i]/nu);
    for(std::size_t e=q.offsets[i];e<q.offsets[i+1];++e)
     next[q.columns[e]]+=v[i]*static_cast<long double>(q.rates[e])/nu;
   }
   v.swap(next);weight*=x/static_cast<long double>(k);
   for(std::size_t i=0;i<n;++i)answer[i]+=weight*v[i];
   for(std::size_t j=0;j<nr;++j)reward_answer[j]+=weight*cumulative[j];
   ++out.terms;
   if(static_cast<long double>(k+2)>x) {
    const long double next_weight=weight*x/static_cast<long double>(k+1);
    tail=next_weight/(1-x/static_cast<long double>(k+2));
    previous_tail=weight+tail;
    if(tail<=local_tol && previous_tail<=local_tol) {converged=true;break;}
   }
  }
  if(!converged)throw std::runtime_error("sparse uniformization: max_terms exceeded");
  long double answer_mass=0;for(long double z:answer)answer_mass+=z;
  if(!(answer_mass>0))throw std::runtime_error("sparse uniformization: lost probability mass");
  for(std::size_t i=0;i<n;++i)p[i]=answer[i]*mass/answer_mass;
  // Normalizing a truncated stochastic series changes its L1 error from
  // at most tail to at most 2*tail. Rewards are deliberately not renormalized.
  for(std::size_t j=0;j<nr;++j) {
   reward_total[j]+=reward_answer[j];
   reward_error[j]+=max_reward[j]*h*(mass*previous_tail+probability_error);
  }
  probability_error+=2*mass*tail;
 }
 out.probability.resize(n);for(std::size_t i=0;i<n;++i)out.probability[i]=static_cast<double>(p[i]);
 for(std::size_t j=0;j<nr;++j) {
  if(!std::isfinite(reward_total[j]) || std::abs(reward_total[j])>std::numeric_limits<double>::max())
   throw std::runtime_error("sparse uniformization: integrated reward overflow");
  out.reward[j]=static_cast<double>(reward_total[j]);out.reward_error_bound[j]=static_cast<double>(reward_error[j]);
 }
 out.error_bound=static_cast<double>(probability_error);
 // This is an operation-count roundoff diagnostic, separate from the rigorous
 // real-arithmetic Poisson truncation bound. It is not an interval-arithmetic certificate.
 out.roundoff_bound=static_cast<double>(std::min(1.0L,
  8*std::numeric_limits<long double>::epsilon()*static_cast<long double>(out.terms)*
  static_cast<long double>(q.rates.size()+n)+4*std::numeric_limits<double>::epsilon()*static_cast<long double>(n)));
 return out;
}
}
#endif
