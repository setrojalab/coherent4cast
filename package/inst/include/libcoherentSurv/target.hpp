#ifndef LIBcoherentSurv_TARGET_HPP
#define LIBcoherentSurv_TARGET_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace coherentSurv {
using TargetMatrix = std::vector<std::vector<double>>;
inline std::vector<double> target_solve(TargetMatrix a, std::vector<double> b) {
  const std::size_t n=b.size();
  for(std::size_t k=0;k<n;++k){
    std::size_t p=k;for(std::size_t i=k+1;i<n;++i)if(std::abs(a[i][k])>std::abs(a[p][k]))p=i;
    if(std::abs(a[p][k])<1e-14)throw std::runtime_error("targeting information is singular; reduce dependent features");
    std::swap(a[p],a[k]);std::swap(b[p],b[k]);const double v=a[k][k];
    for(std::size_t j=k;j<n;++j){a[k][j]/=v;}b[k]/=v;
    for(std::size_t i=0;i<n;++i)if(i!=k){const double z=a[i][k];for(std::size_t j=k;j<n;++j)a[i][j]-=z*a[k][j];b[i]-=z*b[k];}
  }return b;
}
struct TargetEvaluation {double loss=0;TargetMatrix q,mu,H;std::vector<double> gradient;};
inline TargetEvaluation target_evaluate(const TargetMatrix& counts,const TargetMatrix& initial,const TargetMatrix& f,const std::vector<double>& recording,const std::vector<double>& epsilon,double n){
  const std::size_t g=counts.size(),j=f.size(),d=epsilon.size();TargetEvaluation out;
  out.q.assign(g,std::vector<double>(j));out.mu.assign(g,std::vector<double>(d));out.H.assign(d,std::vector<double>(d));out.gradient.assign(d,0);
  for(std::size_t a=0;a<g;++a){
    const double c=1/recording[a];double mx=-std::numeric_limits<double>::infinity(),nr=0;
    for(std::size_t b=0;b<j;++b){double z=std::log(initial[a][b]);for(std::size_t k=0;k<d;++k)z+=c*f[b][k]*epsilon[k];out.q[a][b]=z;mx=std::max(mx,z);nr+=counts[a][b];}
    double den=0;for(double& z:out.q[a]){z=std::exp(z-mx);den+=z;}
    for(std::size_t b=0;b<j;++b){out.q[a][b]/=den;out.loss-=counts[a][b]*std::log(std::max(out.q[a][b],1e-300))/n;for(std::size_t k=0;k<d;++k)out.mu[a][k]+=out.q[a][b]*f[b][k];}
    for(std::size_t k=0;k<d;++k){double observed=0;for(std::size_t b=0;b<j;++b)observed+=counts[a][b]*f[b][k];out.gradient[k]+=c*(nr*out.mu[a][k]-observed)/n;
      for(std::size_t l=0;l<d;++l)for(std::size_t b=0;b<j;++b)out.H[k][l]+=nr*c*c*out.q[a][b]*(f[b][k]-out.mu[a][k])*(f[b][l]-out.mu[a][l])/n;
    }
  }return out;
}
struct TargetResult {TargetMatrix q,mu;std::vector<double> epsilon,gradient;int iterations=0;bool converged=false;};
inline TargetResult target_fit(const TargetMatrix& counts,const TargetMatrix& initial,const TargetMatrix& features,const std::vector<double>& recording,double n,int maxit=60,double tolerance=0){
  if(counts.empty()||features.empty()||features[0].empty()||!std::isfinite(n)||n<=1||maxit<1)throw std::invalid_argument("invalid targeting dimensions, sample size or iteration limit");
  const std::size_t g=counts.size(),j=features.size(),d=features[0].size();
  if(initial.size()!=g||recording.size()!=g)throw std::invalid_argument("targeting group dimension mismatch");
  for(const auto& row:features){if(row.size()!=d)throw std::invalid_argument("ragged feature matrix");for(double x:row)if(!std::isfinite(x))throw std::invalid_argument("nonfinite target feature");}
  double total_recorded=0;
  for(std::size_t a=0;a<g;++a){
    if(counts[a].size()!=j||initial[a].size()!=j||!std::isfinite(recording[a])||recording[a]<=0||recording[a]>1)throw std::invalid_argument("invalid targeting group or recording probability");
    double mass=0;for(std::size_t b=0;b<j;++b){if(!std::isfinite(initial[a][b])||initial[a][b]<=0||!std::isfinite(counts[a][b])||counts[a][b]<0)throw std::invalid_argument("initial law must be strictly positive and counts nonnegative");mass+=initial[a][b];total_recorded+=counts[a][b];}if(std::abs(mass-1)>1e-9)throw std::invalid_argument("initial path law does not sum to one");
  }
  if(total_recorded==0)throw std::invalid_argument("no recorded paths available for targeting");
  if(tolerance<=0){tolerance=1e-6/std::sqrt(n);}TargetResult result;result.epsilon.assign(d,0);
  TargetEvaluation ev=target_evaluate(counts,initial,features,recording,result.epsilon,n);
  for(int it=0;it<maxit;++it){result.iterations=it;double score=0;for(double x:ev.gradient)score=std::max(score,std::abs(x));if(score<tolerance){result.converged=true;break;}
    std::vector<double> step=target_solve(ev.H,ev.gradient);double dot=0;for(std::size_t k=0;k<d;++k)dot+=ev.gradient[k]*step[k];double scale=1;bool accepted=false;
    while(scale>1e-10){std::vector<double> next=result.epsilon;for(std::size_t k=0;k<d;++k)next[k]=std::max(-30.0,std::min(30.0,next[k]-scale*step[k]));TargetEvaluation candidate=target_evaluate(counts,initial,features,recording,next,n);if(candidate.loss<=ev.loss-1e-4*scale*dot){result.epsilon=next;ev=candidate;accepted=true;break;}scale*=.5;}
    if(!accepted)break;
  }
  result.q=ev.q;result.mu=ev.mu;result.gradient=ev.gradient;double score=0;for(double x:ev.gradient)score=std::max(score,std::abs(x));result.converged=score<tolerance;return result;
}
}
#endif
