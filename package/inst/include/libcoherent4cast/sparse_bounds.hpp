#ifndef LIBcoherent4cast_SPARSE_BOUNDS_HPP
#define LIBcoherent4cast_SPARSE_BOUNDS_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
namespace coherent4cast {
struct SparseBoundEndpoint {
 double value=0,certificate_lower=0,certificate_upper=0,dual_gap=0,residual=0;
 double linear_primal=0,linear_dual=0,dual_mass=0,dual_high_mass=0,denominator=0;
 double linearization_point=0,objective_sign=1;
 std::size_t iterations=0;
 std::vector<double> law,recorded,unrecorded,fibre_dual;
};
struct SparseBoundResult {
 SparseBoundEndpoint lower,upper;
 double denominator_lower_bound=0;
 std::size_t support=0,fibres=0,workspace_bytes=0;
};
struct SparseLinearResult {
 std::vector<long double> recorded,unrecorded,fibre_dual;
 long double primal=0,dual=0,dual_mass=0,dual_high_mass=0;
};
// Exact structural reduction of the marginal-odds model, not a general LP.
// a=x+y, b=l*x+u*y; total y=(1-rho-l*rho)/(u-l). Within
// each fibre both x and y maximize the same linear objective c_j.
inline SparseLinearResult sparse_linear_bound(const std::vector<long double>& c,
 const std::vector<std::size_t>& map,const std::vector<double>& lower,
 const std::vector<double>& upper,long double rho,long double l,long double u) {
 const std::size_t j=c.size(),m=lower.size();
 std::vector<long double> best(m,-std::numeric_limits<long double>::infinity()),w(m),y(m,0);
 std::vector<std::size_t> atom(m,0),order(m);std::iota(order.begin(),order.end(),0);
 for(std::size_t k=0;k<j;++k)if(c[k]>best[map[k]]){best[map[k]]=c[k];atom[map[k]]=k;}
 std::stable_sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return best[a]>best[b];});
 long double extra=rho,high=(u==l)?0:(1-rho-l*rho)/(u-l);
 high=std::max(0.0L,std::min(rho,high));
 for(std::size_t h=0;h<m;++h){w[h]=lower[h];extra-=w[h];}
 long double cw=best[order.front()]+1;
 for(std::size_t z=0;z<m && extra>0;++z) {
  const auto h=order[z];const long double add=std::min(extra,static_cast<long double>(upper[h])-w[h]);
  if(add>0){w[h]+=add;extra-=add;cw=best[h];}
 }
 if(extra>1e-10L)throw std::runtime_error("sparse bounds: incompatible fibre masses");
 const long double high_total=high;long double cy=best[order.front()]+1;
 for(std::size_t z=0;z<m && high>0;++z){const auto h=order[z];const long double add=std::min(high,w[h]);if(add>0){y[h]=add;high-=add;cy=best[h];}}
 const long double A=1+l,B=1+u,mu=(B-A)*cy;
 const long double lambda=A*cw+std::max(0.0L,(B-A)*cw-mu);
 SparseLinearResult out;out.recorded.assign(j,0);out.unrecorded.assign(j,0);out.fibre_dual.resize(m);
 out.dual_mass=lambda;out.dual_high_mass=mu;
 out.dual=lambda*rho+mu*high_total;
 for(std::size_t h=0;h<m;++h) {
  const auto k=atom[h];out.recorded[k]=w[h];out.unrecorded[k]=l*(w[h]-y[h])+u*y[h];
  out.primal+=c[k]*(out.recorded[k]+out.unrecorded[k]);
  const long double d=A*best[h]+std::max(0.0L,(B-A)*best[h]-mu)-lambda;
  out.fibre_dual[h]=d;out.dual+=(d>0?upper[h]:lower[h])*d;
 }
 return out;
}
inline SparseBoundResult compatible_bounds_sparse(const std::vector<double>& p,
 const std::vector<double>& numerator,const std::vector<double>& denominator,
 double kappa,double rho,const std::vector<double>& lower,const std::vector<double>& upper,
 const std::vector<std::size_t>& map,double tolerance=1e-10,std::size_t max_iterations=1000) {
 const std::size_t j=numerator.size(),m=p.size();
 if(j==0 || m==0 || denominator.size()!=j || map.size()!=j || lower.size()!=m || upper.size()!=m ||
    !std::isfinite(kappa) || kappa<1 || !std::isfinite(rho) || rho<=0 || rho>=1 ||
    !std::isfinite(tolerance) || tolerance<=0 || max_iterations<1)
  throw std::invalid_argument("sparse bounds: invalid dimensions or control values");
 const long double ratio=(1-static_cast<long double>(rho))/rho,l=ratio/kappa,u=ratio*kappa;
 if(!std::isfinite(u))throw std::invalid_argument("sparse bounds: odds bound overflow");
 std::vector<double> min_den(m,std::numeric_limits<double>::infinity());
 long double pmass=0,lmass=0,umass=0,den_lower=0;
 for(std::size_t k=0;k<j;++k) {
  if(map[k]>=m || !std::isfinite(numerator[k]) || !std::isfinite(denominator[k]) || denominator[k]<0)
   throw std::invalid_argument("sparse bounds: invalid reward, denominator or fibre map");
  min_den[map[k]]=std::min(min_den[map[k]],denominator[k]);
 }
 for(std::size_t h=0;h<m;++h) {
  if(!std::isfinite(p[h]) || p[h]<0 || !std::isfinite(lower[h]) || !std::isfinite(upper[h]) ||
     lower[h]<0 || lower[h]>upper[h] || upper[h]>rho+1e-10 || !std::isfinite(min_den[h]))
   throw std::invalid_argument("sparse bounds: invalid recorded law, limits or empty fibre");
  pmass+=p[h];lmass+=lower[h];umass+=upper[h];den_lower+=lower[h]*min_den[h];
 }
 if(std::abs(pmass-1)>1e-9L || lmass>rho+1e-10L || umass<rho-1e-10L)
  throw std::invalid_argument("sparse bounds: inconsistent recorded probability masses");
 if(den_lower<=0)throw std::invalid_argument("sparse bounds: conditioning denominator lacks a positive recorded lower bound");
 SparseBoundResult out;out.denominator_lower_bound=static_cast<double>(den_lower);
 out.support=j;out.fibres=m;out.workspace_bytes=(8*j+8*m)*sizeof(long double)+(8*j+3*m)*sizeof(double)+(j+2*m)*sizeof(std::size_t);
 auto solve=[&](long double sign) {
  SparseBoundEndpoint endpoint;std::vector<long double> c(j,0);
  auto linear=sparse_linear_bound(c,map,lower,upper,rho,l,u);
  long double r=0,den=0,num=0;
  for(std::size_t k=0;k<j;++k){const long double q=linear.recorded[k]+linear.unrecorded[k];num+=sign*numerator[k]*q;den+=denominator[k]*q;}
  r=num/den;bool converged=false;
  for(std::size_t iteration=1;iteration<=max_iterations;++iteration) {
   for(std::size_t k=0;k<j;++k)c[k]=sign*numerator[k]-r*denominator[k];
   linear=sparse_linear_bound(c,map,lower,upper,rho,l,u);
   num=0;den=0;for(std::size_t k=0;k<j;++k){const long double q=linear.recorded[k]+linear.unrecorded[k];num+=sign*numerator[k]*q;den+=denominator[k]*q;}
   if(den<den_lower*(1-1e-9L))throw std::runtime_error("sparse bounds: denominator certificate failed");
   const long double next=num/den;
   const long double upper_r=std::max(next,r+std::max(0.0L,linear.dual)/den_lower);
   const long double gap=std::max(0.0L,upper_r-next);
   endpoint.iterations=iteration;
   if(gap<=tolerance*(1+std::abs(next))) {
    endpoint.value=static_cast<double>(sign*next);
    endpoint.certificate_lower=static_cast<double>(sign>0?next:-upper_r);
    endpoint.certificate_upper=static_cast<double>(sign>0?upper_r:-next);
    endpoint.dual_gap=static_cast<double>(gap);converged=true;break;
   }
   r=next;
  }
  if(!converged)throw std::runtime_error("sparse bounds: fractional iteration limit exceeded");
  endpoint.denominator=static_cast<double>(den);
  endpoint.linearization_point=static_cast<double>(r);endpoint.objective_sign=static_cast<double>(sign);
  endpoint.linear_primal=static_cast<double>(linear.primal);endpoint.linear_dual=static_cast<double>(linear.dual);
  endpoint.dual_mass=static_cast<double>(linear.dual_mass);endpoint.dual_high_mass=static_cast<double>(linear.dual_high_mass);
  endpoint.law.resize(j);endpoint.recorded.resize(j);endpoint.unrecorded.resize(j);endpoint.fibre_dual.resize(m);
  std::vector<long double> fibre_mass(m,0);long double aa=0,bb=0,residual=0;
  for(std::size_t k=0;k<j;++k) {
   const long double a=linear.recorded[k],b=linear.unrecorded[k];aa+=a;bb+=b;fibre_mass[map[k]]+=a;
   residual=std::max(residual,std::max(l*a-b,b-u*a));
   endpoint.recorded[k]=static_cast<double>(a);endpoint.unrecorded[k]=static_cast<double>(b);endpoint.law[k]=static_cast<double>(a+b);
  }
  residual=std::max(residual,std::max(std::abs(aa-rho),std::abs(bb-(1-rho))));
  for(std::size_t h=0;h<m;++h){residual=std::max(residual,std::max(static_cast<long double>(lower[h])-fibre_mass[h],fibre_mass[h]-upper[h]));endpoint.fibre_dual[h]=static_cast<double>(linear.fibre_dual[h]);}
  endpoint.residual=static_cast<double>(residual);
  if(residual>1e-8L || linear.dual+1e-10L<linear.primal)throw std::runtime_error("sparse bounds: primal-dual verification failed");
  return endpoint;
 };
 out.lower=solve(-1);out.upper=solve(1);return out;
}
}
#endif
