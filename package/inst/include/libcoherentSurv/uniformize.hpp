#ifndef LIBcoherentSurv_UNIFORMIZE_HPP
#define LIBcoherentSurv_UNIFORMIZE_HPP
#include "gamma_lik.hpp"
namespace coherentSurv {
struct PropagationResult {
 std::vector<double> probability, reward, reward_error_bound;
 double error_bound=0.0,roundoff_bound=0.0;
 std::size_t terms=0,squarings=0;
};
inline std::vector<long double> matrix_product(const std::vector<long double>& a,
 const std::vector<long double>& b,std::size_t rows,std::size_t inner,std::size_t cols) {
 std::vector<long double> c(rows*cols,0.0L);
 for (std::size_t i=0;i<rows;++i) for (std::size_t k=0;k<inner;++k) {
  const long double v=a[i*inner+k];
  if (v!=0.0L) for (std::size_t j=0;j<cols;++j) c[i*cols+j]+=v*b[k*cols+j];
 }
 return c;
}
inline void normalize_transition(std::vector<long double>& matrix,std::size_t n) {
 for(std::size_t i=0;i<n;++i) {
  long double sum=0.0L;for(std::size_t j=0;j<n;++j)sum+=matrix[i*n+j];
  if(!(sum>0.0L) || !std::isfinite(sum))throw std::runtime_error("uniformize: lost transition mass");
  for(std::size_t j=0;j<n;++j)matrix[i*n+j]/=sum;
 }
}
// Ordered segments are composed in R without remixture of latent frailty.
// Error bounds below bound series truncation. Floating-point error is reported separately.
inline PropagationResult uniformize(const std::vector<double>& probability,
 const std::vector<double>& generator,double duration,const std::vector<double>& rewards,
 std::size_t nr,double tolerance=1e-12,std::size_t max_terms=10000U) {
 const std::size_t n=probability.size();
 if (n==0 || n>512U || generator.size()!=n*n || rewards.size()!=n*nr ||
     !std::isfinite(duration) || duration<0.0 || !std::isfinite(tolerance) ||
     tolerance<=0.0 || tolerance>=1.0 || max_terms<2U)
  throw std::invalid_argument("uniformize: invalid dimensions, duration, tolerance or budget");
 long double mass=0.0L,nu=0.0L;std::vector<long double> bound(nr,0.0L);
 for (std::size_t i=0;i<n;++i) {
  if (!std::isfinite(probability[i]) || probability[i]<0.0) throw std::invalid_argument("uniformize: invalid initial probability");
  mass+=probability[i];long double rowsum=0.0L,scale=1.0L;
  for (std::size_t j=0;j<n;++j) {
   const double q=generator[i*n+j];
   if (!std::isfinite(q) || (i!=j && q<0.0) || (i==j && q>0.0)) throw std::invalid_argument("uniformize: invalid generator sign or nonfinite entry");
   rowsum+=q;scale+=std::abs(q);
  }
  if (std::abs(rowsum)>1e-12L*scale) throw std::invalid_argument("uniformize: generator rows must sum to zero");
  nu=std::max(nu,-static_cast<long double>(generator[i*n+i]));
  for (std::size_t j=0;j<nr;++j) {
   if (!std::isfinite(rewards[i*nr+j])) throw std::invalid_argument("uniformize: nonfinite reward");
   bound[j]=std::max(bound[j],std::abs(static_cast<long double>(rewards[i*nr+j])));
  }
 }
 if (mass<=0.0L || mass>1.0L+1e-10L) throw std::invalid_argument("uniformize: initial mass must be in (0,1]");
 PropagationResult out;out.reward.assign(nr,0.0);out.reward_error_bound.assign(nr,0.0);
 out.probability=probability;
 if (duration==0.0) return out;
 if (nu==0.0L) {
  for (std::size_t j=0;j<nr;++j) for (std::size_t i=0;i<n;++i)
   out.reward[j]+=probability[i]*rewards[i*nr+j]*duration;
  for (double v:out.reward) if (!std::isfinite(v)) throw std::overflow_error("uniformize: reward overflow");
  return out;
 }
 long double h=duration,x=nu*h;std::size_t squarings=0;
 while (x>16.0L) {h*=0.5L;x*=0.5L;++squarings;if (squarings>1024U) throw std::overflow_error("uniformize: scaling range exceeded");}
 const long double amplification=std::ldexp(1.0L,static_cast<int>(squarings));
 const long double small_tol=static_cast<long double>(tolerance)/(4.0L*amplification);
 std::vector<long double> P(n*n,0.0L),power(n*n,0.0L),T(n*n,0.0L),R(n*nr),cumulative(n*nr,0.0L),G(n*nr,0.0L);
 for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<n;++j) {
  P[i*n+j]=static_cast<long double>(generator[i*n+j])/nu+(i==j?1.0L:0.0L);
  power[i*n+j]=(i==j?1.0L:0.0L);
 }
 for (std::size_t j=0;j<R.size();++j) R[j]=static_cast<long double>(rewards[j])/nu;
 long double weight=std::exp(-x),tail=1.0L;std::size_t last=0;
 for (std::size_t j=0;j<T.size();++j) T[j]=weight*power[j];
 for (std::size_t k=1;k<=max_terms;++k) {
  if (nr>0U) {
   const auto increment=matrix_product(power,R,n,n,nr);
   for (std::size_t j=0;j<cumulative.size();++j) cumulative[j]+=increment[j];
  }
  power=matrix_product(power,P,n,n,n);weight*=x/static_cast<long double>(k);
  for (std::size_t j=0;j<T.size();++j) T[j]+=weight*power[j];
  for (std::size_t j=0;j<G.size();++j) G[j]+=weight*cumulative[j];
  if (static_cast<long double>(k+2U)>x) {
   const long double next=weight*x/static_cast<long double>(k+1U);
   tail=next/(1.0L-x/static_cast<long double>(k+2U));
   // reward tail <= max |r| * h * P{Poisson(x)>=k}.
   if (tail<small_tol && (weight+tail)<small_tol) {last=k;break;}
  }
 }
 if (last==0U) throw std::runtime_error("uniformize: Poisson truncation budget exceeded");
 // Normalizing the truncated stochastic row costs at most a further tail in L1.
 // It prevents epsilon-scale mass loss from being amplified at extreme rates.
 normalize_transition(T,n);
 const long double probability_error=2.0L*amplification*tail;
 std::vector<long double> reward_error(nr);
 for (std::size_t j=0;j<nr;++j) reward_error[j]=bound[j]*h*(weight+tail);
 long double accumulated_tail=2.0L*tail,elapsed=h;
 for (std::size_t s=0;s<squarings;++s) {
  if (nr>0U) {
   const auto second=matrix_product(T,G,n,n,nr);
   for (std::size_t j=0;j<G.size();++j) G[j]+=second[j];
   for (std::size_t j=0;j<nr;++j) reward_error[j]=2.0L*reward_error[j]+accumulated_tail*elapsed*bound[j];
  }
  T=matrix_product(T,T,n,n,n);normalize_transition(T,n);accumulated_tail*=2.0L;elapsed*=2.0L;
 }
 std::fill(out.probability.begin(),out.probability.end(),0.0);
 for (std::size_t j=0;j<n;++j) {
  long double value=0.0L;for (std::size_t i=0;i<n;++i)value+=probability[i]*T[i*n+j];
  out.probability[j]=finite_double(value,"uniformize: propagated probability overflow");
 }
 for (std::size_t j=0;j<nr;++j) {
  long double value=0.0L;for (std::size_t i=0;i<n;++i)value+=probability[i]*G[i*nr+j];
  out.reward[j]=finite_double(value,"uniformize: integrated reward overflow");
  out.reward_error_bound[j]=finite_double(mass*reward_error[j],"uniformize: reward bound overflow");
 }
 out.error_bound=finite_double(mass*probability_error,"uniformize: error bound overflow");
 out.roundoff_bound=finite_double(std::min(1.0L,amplification*static_cast<long double>(last+squarings+1U)*static_cast<long double>(n)*std::numeric_limits<long double>::epsilon()*16.0L),"uniformize: roundoff bound overflow");
 out.terms=last+1U;out.squarings=squarings;return out;
}
}
#endif
