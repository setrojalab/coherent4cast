#ifndef LIBRECURCAST_CONTINUOUS_TARGET_HPP
#define LIBRECURCAST_CONTINUOUS_TARGET_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace recurcast {
struct CTRecord {std::size_t subject,group,state; double start,stop; int event;};
struct CTResult {
 std::vector<double> estimate,initial,eif,covariance,se,rates,values,polynomial,score_trace;
 double score=0,interpolation_error=0; int iterations=0; bool converged=false;
};
// Symmetric Moore-Penrose solve for the small target-coordinate information.
inline std::vector<double> ct_psolve(std::vector<double> a,const std::vector<double>& b){
 const std::size_t d=b.size();std::vector<double> v(d*d,0),out(d,0);for(std::size_t i=0;i<d;++i)v[i*d+i]=1;
 for(std::size_t it=0;it<100*d*d;++it){std::size_t p=0,q=0;double largest=0,diag=0;
  for(std::size_t i=0;i<d;++i){diag=std::max(diag,std::abs(a[i*d+i]));for(std::size_t j=i+1;j<d;++j)if(std::abs(a[i*d+j])>largest){largest=std::abs(a[i*d+j]);p=i;q=j;}}
  if(largest<1e-14*std::max(diag,1e-30))break;
  const double angle=.5*std::atan2(2*a[p*d+q],a[q*d+q]-a[p*d+p]),c=std::cos(angle),s=std::sin(angle);
  for(std::size_t i=0;i<d;++i){const double x=a[i*d+p],y=a[i*d+q];a[i*d+p]=c*x-s*y;a[i*d+q]=s*x+c*y;}
  for(std::size_t i=0;i<d;++i){const double x=a[p*d+i],y=a[q*d+i];a[p*d+i]=c*x-s*y;a[q*d+i]=s*x+c*y;}
  for(std::size_t i=0;i<d;++i){const double x=v[i*d+p],y=v[i*d+q];v[i*d+p]=c*x-s*y;v[i*d+q]=s*x+c*y;}
 }
 double top=0;for(std::size_t k=0;k<d;++k)top=std::max(top,a[k*d+k]);
 for(std::size_t k=0;k<d;++k)if(a[k*d+k]>top*1e-11&&a[k*d+k]>0){double w=0;for(std::size_t j=0;j<d;++j)w+=v[j*d+k]*b[j];w/=a[k*d+k];for(std::size_t j=0;j<d;++j)out[j]+=v[j*d+k]*w;}
 return out;
}
class ContinuousTarget {
 std::size_t groups_,bins_,states_,marks_,dim_,n_;double tau_,h_;
 std::vector<std::size_t> next_,subject_group_;
 std::vector<CTRecord> rows_;
 std::vector<double> censor_,exposure_,value_,poly_,integrand_,prefix_;
 std::size_t li(std::size_t g,std::size_t b,std::size_t q,std::size_t k) const{return ((g*bins_+b)*states_+q)*marks_+k;}
 std::size_t vi(std::size_t g,std::size_t b,std::size_t q,std::size_t j) const{return ((g*(bins_+1)+b)*(states_+1)+q)*dim_+j;}
 std::size_t pi(std::size_t g,std::size_t b,std::size_t q,std::size_t j,std::size_t power=0) const{return (((g*bins_+b)*(states_+1)+q)*dim_+j)*4+power;}
 std::vector<double> action(const std::vector<double>& rates,std::size_t g,std::size_t b,const std::vector<double>& v) const{
  std::vector<double> out((states_+1)*dim_,0);
  for(std::size_t q=0;q<states_;++q){for(std::size_t k=0;k<marks_;++k){const double rate=rates[li(g,b,q,k)];for(std::size_t j=0;j<dim_;++j)out[q*dim_+j]+=rate*(v[next_[q*marks_+k]*dim_+j]-v[q*dim_+j]);if(k+1<marks_)out[q*dim_+k]+=rate;}out[q*dim_+dim_-1]+=1;}
  return out;
 }
 std::vector<double> propagate(const std::vector<double>& rates,std::size_t g,std::size_t b,std::vector<double> v,double duration) const{
  double nu=0;for(std::size_t q=0;q<states_;++q){double sum=0;for(std::size_t k=0;k<marks_;++k)sum+=rates[li(g,b,q,k)];nu=std::max(nu,sum);}
  if(nu==0){for(std::size_t q=0;q<states_;++q)v[q*dim_+dim_-1]+=duration;return v;}
  const double parts_d=std::max(1.,std::ceil(nu*duration/12.));if(parts_d>10000)throw std::runtime_error("continuous targeting propagation budget exceeded");const auto parts=static_cast<std::size_t>(parts_d);const double step=duration/static_cast<double>(parts);
  for(std::size_t part=0;part<parts;++part){std::vector<double> term=v,answer(v.size());double weight=std::exp(-nu*step),total=weight;for(std::size_t j=0;j<v.size();++j)answer[j]=weight*term[j];bool done=false;
   for(int m=1;m<300;++m){const auto a=action(rates,g,b,term);for(std::size_t j=0;j<v.size();++j)term[j]+=a[j]/nu;weight*=nu*step/static_cast<double>(m);total+=weight;for(std::size_t j=0;j<v.size();++j)answer[j]+=weight*term[j];if(1-total<2e-15){done=true;break;}}
   if(!done){throw std::runtime_error("continuous targeting uniformization did not converge");}v=answer;
  }return v;
 }
 std::vector<double> moments(double u,double a) const{
  std::vector<double> out(4,0);double w=1;
  for(int m=0;m<40;++m){for(int k=0;k<4;++k)out[static_cast<std::size_t>(k)]+=w*std::pow(u,k+m+1)/static_cast<double>(k+m+1);w*=a/static_cast<double>(m+1);}return out;
 }
 std::vector<double> value_at(std::size_t g,std::size_t q,double time) const{
  const std::size_t b=std::min(static_cast<std::size_t>(time/h_),bins_-1);const double u=time/h_-static_cast<double>(b);std::vector<double> out(dim_);
  for(std::size_t j=0;j<dim_;++j){const std::size_t p=pi(g,b,q,j);out[j]=poly_[p]+u*(poly_[p+1]+u*(poly_[p+2]+u*poly_[p+3]));}return out;
 }
 std::vector<double> cumulative(std::size_t g,std::size_t q,double time) const{
  const std::size_t b=std::min(static_cast<std::size_t>(time/h_),bins_-1);const double u=time/h_-static_cast<double>(b);const auto moment=moments(u,censor_[g]*h_);std::vector<double> out(dim_);
  for(std::size_t j=0;j<dim_;++j){out[j]=prefix_[vi(g,b,q,j)];double sum=0;for(std::size_t k=0;k<4;++k)sum+=integrand_[pi(g,b,q,j,k)]*moment[k];out[j]+=h_*std::exp(censor_[g]*static_cast<double>(b)*h_)*sum;}return out;
 }
 void construct(const std::vector<double>& rates){
  value_.assign(groups_*(bins_+1)*(states_+1)*dim_,0);poly_.assign(groups_*bins_*(states_+1)*dim_*4,0);integrand_=poly_;prefix_=value_;
  for(std::size_t g=0;g<groups_;++g){for(std::size_t b=0;b<=bins_;++b)value_[vi(g,b,states_,dim_-2)]=1;
   for(std::size_t bb=bins_;bb>0;--bb){const std::size_t b=bb-1;std::vector<double> right((states_+1)*dim_);for(std::size_t q=0;q<=states_;++q)for(std::size_t j=0;j<dim_;++j)right[q*dim_+j]=value_[vi(g,b+1,q,j)];const auto left=propagate(rates,g,b,right,h_),a0=action(rates,g,b,left),a1=action(rates,g,b,right);
    for(std::size_t q=0;q<=states_;++q)for(std::size_t j=0;j<dim_;++j){const std::size_t u=q*dim_+j,p=pi(g,b,q,j);const double start=left[u],end=right[u],d0=-h_*a0[u],d1=-h_*a1[u];value_[vi(g,b,q,j)]=start;poly_[p]=start;poly_[p+1]=d0;poly_[p+2]=3*(end-start)-2*d0-d1;poly_[p+3]=2*(start-end)+d0+d1;}
    for(std::size_t q=0;q<states_;++q)for(std::size_t k=0;k<marks_;++k){const double rate=rates[li(g,b,q,k)];for(std::size_t j=0;j<dim_;++j)for(std::size_t p=0;p<4;++p)integrand_[pi(g,b,q,j,p)]+=rate*(poly_[pi(g,b,next_[q*marks_+k],j,p)]-poly_[pi(g,b,q,j,p)]);if(k+1<marks_)integrand_[pi(g,b,q,k)]+=rate;}
   }
   const auto moment=moments(1,censor_[g]*h_);
   for(std::size_t b=0;b<bins_;++b)for(std::size_t q=0;q<=states_;++q)for(std::size_t j=0;j<dim_;++j){double sum=0;for(std::size_t p=0;p<4;++p)sum+=integrand_[pi(g,b,q,j,p)]*moment[p];prefix_[vi(g,b+1,q,j)]=prefix_[vi(g,b,q,j)]+h_*std::exp(censor_[g]*static_cast<double>(b)*h_)*sum;}
  }
 }
 void evaluate(const std::vector<double>& rates,CTResult& out,std::vector<double>& information,std::vector<double>& clever){
  construct(rates);out.estimate.assign(dim_,0);out.eif.assign(n_*dim_,0);
  for(std::size_t i=0;i<n_;++i)for(std::size_t j=0;j<dim_;++j){const double z=value_[vi(subject_group_[i],0,0,j)];out.eif[i*dim_+j]=z;out.estimate[j]+=z/static_cast<double>(n_);}
  for(std::size_t i=0;i<n_;++i)for(std::size_t j=0;j<dim_;++j)out.eif[i*dim_+j]-=out.estimate[j];
  for(const auto& row:rows_){const auto a=cumulative(row.group,row.state,row.start),b=cumulative(row.group,row.state,row.stop);for(std::size_t j=0;j<dim_;++j)out.eif[row.subject*dim_+j]-=b[j]-a[j];
   if(row.event>=0&&static_cast<std::size_t>(row.event)<marks_){const std::size_t k=static_cast<std::size_t>(row.event);const auto before=value_at(row.group,row.state,row.stop),after=value_at(row.group,next_[row.state*marks_+k],row.stop);for(std::size_t j=0;j<dim_;++j)out.eif[row.subject*dim_+j]+=std::exp(censor_[row.group]*row.stop)*(after[j]-before[j]+((j==k&&k+1<marks_)?1.:0.));}
  }
  information.assign(dim_*dim_,0);clever.assign(rates.size()*dim_,0);
  for(std::size_t g=0;g<groups_;++g)for(std::size_t b=0;b<bins_;++b){const double t=(static_cast<double>(b)+.5)*h_;
   for(std::size_t q=0;q<states_;++q){const auto before=value_at(g,q,t);for(std::size_t k=0;k<marks_;++k){const auto after=value_at(g,next_[q*marks_+k],t);const std::size_t index=li(g,b,q,k);for(std::size_t j=0;j<dim_;++j)clever[index*dim_+j]=std::exp(censor_[g]*t)*(after[j]-before[j]+((j==k&&k+1<marks_)?1.:0.));const double w=exposure_[(g*bins_+b)*states_+q]*rates[index]/static_cast<double>(n_);for(std::size_t j=0;j<dim_;++j)for(std::size_t k2=0;k2<dim_;++k2)information[j*dim_+k2]+=w*clever[index*dim_+j]*clever[index*dim_+k2];}}
  }
 }
public:
 ContinuousTarget(std::size_t groups,std::size_t bins,std::size_t states,std::size_t marks,double tau,std::vector<std::size_t> next,std::vector<std::size_t> subject_group,std::vector<CTRecord> rows,std::vector<double> censor):groups_(groups),bins_(bins),states_(states),marks_(marks),dim_(marks+1),n_(subject_group.size()),tau_(tau),h_(tau/static_cast<double>(bins)),next_(std::move(next)),subject_group_(std::move(subject_group)),rows_(std::move(rows)),censor_(std::move(censor)){
  if(!groups_||!bins_||!states_||!marks_||marks_>7||n_<2||!std::isfinite(tau_)||tau_<=0||next_.size()!=states_*marks_||censor_.size()!=groups_)throw std::invalid_argument("invalid continuous targeting dimensions");
  for(auto q:next_){if(q>states_)throw std::invalid_argument("transition outside continuous state space");}for(auto g:subject_group_)if(g>=groups_)throw std::invalid_argument("subject group outside continuous nuisance folds");
  for(double g:censor_)if(!std::isfinite(g)||g<0||g*tau_>25)throw std::invalid_argument("censor survival outside numerical positivity range");
  exposure_.assign(groups_*bins_*states_,0);
  for(const auto& r:rows_){if(r.subject>=n_||r.group!=subject_group_[r.subject]||r.state>=states_||!std::isfinite(r.start)||!std::isfinite(r.stop)||r.start<0||r.stop<=r.start||r.stop>tau_+1e-12||r.event < -1||r.event>static_cast<int>(marks_))throw std::invalid_argument("invalid continuous observation row");const std::size_t first=std::min(static_cast<std::size_t>(r.start/h_),bins_-1),last=std::min(static_cast<std::size_t>(r.stop/h_),bins_-1);for(std::size_t b=first;b<=last;++b)exposure_[(r.group*bins_+b)*states_+r.state]+=std::max(0.,std::min(r.stop,static_cast<double>(b+1)*h_)-std::max(r.start,static_cast<double>(b)*h_));}
 }
 CTResult fit(std::vector<double> rates,int maxit=60,double tolerance=0){
  if(rates.size()!=groups_*bins_*states_*marks_||maxit<0){throw std::invalid_argument("initial continuous intensity array has incorrect dimensions");}for(double r:rates)if(!std::isfinite(r)||r<0||r>1e4)throw std::invalid_argument("continuous intensities must be finite and within the declared numerical range");if(tolerance<=0)tolerance=1e-7/std::sqrt(static_cast<double>(n_));
  CTResult out;std::vector<double> information,clever,score(dim_);
  for(int iteration=0;iteration<=maxit;++iteration){evaluate(rates,out,information,clever);if(iteration==0)out.initial=out.estimate;std::fill(score.begin(),score.end(),0);for(std::size_t i=0;i<n_;++i)for(std::size_t j=0;j<dim_;++j)score[j]+=out.eif[i*dim_+j]/static_cast<double>(n_);out.score=0;for(double s:score)out.score=std::max(out.score,std::abs(s));out.score_trace.push_back(out.score);out.iterations=iteration;if(out.score<tolerance){out.converged=true;break;}if(iteration==maxit)break;
   const auto step=ct_psolve(information,score);std::vector<double> delta(rates.size(),0);double maxdelta=1e-12;for(std::size_t i=0;i<rates.size();++i){for(std::size_t j=0;j<dim_;++j)delta[i]+=clever[i*dim_+j]*step[j];maxdelta=std::max(maxdelta,std::abs(delta[i]));}const double scale=std::min(1.,.35/maxdelta);for(std::size_t i=0;i<rates.size();++i){rates[i]*=std::exp(scale*delta[i]);if(!std::isfinite(rates[i])||rates[i]>1e4)throw std::runtime_error("continuous intensity update exceeded numerical range");}
  }
  out.covariance.assign(dim_*dim_,0);out.se.resize(dim_);std::vector<double> mean(dim_,0);for(std::size_t i=0;i<n_;++i)for(std::size_t j=0;j<dim_;++j)mean[j]+=out.eif[i*dim_+j]/static_cast<double>(n_);
  for(std::size_t i=0;i<n_;++i){for(std::size_t j=0;j<dim_;++j)for(std::size_t k=0;k<dim_;++k)out.covariance[j*dim_+k]+=(out.eif[i*dim_+j]-mean[j])*(out.eif[i*dim_+k]-mean[k])/(static_cast<double>(n_)*static_cast<double>(n_-1));}for(std::size_t j=0;j<dim_;++j)out.se[j]=std::sqrt(std::max(0.,out.covariance[j*dim_+j]));
  out.interpolation_error=0;for(std::size_t g=0;g<groups_;++g)for(std::size_t b=0;b<bins_;b+=std::max<std::size_t>(1,bins_/4)){std::vector<double> right((states_+1)*dim_);for(std::size_t q=0;q<=states_;++q)for(std::size_t j=0;j<dim_;++j)right[q*dim_+j]=value_[vi(g,b+1,q,j)];const auto exact=propagate(rates,g,b,right,.63*h_);for(std::size_t q=0;q<states_;++q){const auto interpolated=value_at(g,q,(static_cast<double>(b)+.37)*h_);for(std::size_t j=0;j<dim_;++j)out.interpolation_error=std::max(out.interpolation_error,std::abs(interpolated[j]-exact[q*dim_+j]));}}
  out.rates=std::move(rates);out.values=value_;out.polynomial=poly_;return out;
 }
};
}
#endif
