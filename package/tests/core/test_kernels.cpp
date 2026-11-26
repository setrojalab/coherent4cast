#include <libcoherentSurv/gamma_lik.hpp>
#include <libcoherentSurv/coarsen.hpp>
#include <libcoherentSurv/uniformize.hpp>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
void require(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
double maxdiff(const std::vector<double>& a,const std::vector<double>& b) {
 require(a.size()==b.size(),"test dimension mismatch");double out=0.0;
 for(std::size_t j=0;j<a.size();++j) {out=std::max(out,std::abs(a[j]-b[j]));}
 return out;
}
}
int main() {
 try {
 const std::vector<std::size_t> id={0,0,1,1,2,3};
 const std::vector<double> dt={.4,.8,.5,1.,.7,1.2},X={1,-1,1,-1,1,0,1,0,1,.5,1,1},y={0,1,0,1,0,1};
 std::vector<double> theta={-.7,.25,std::log(.6)};
 auto evaluate=[&](const std::vector<double>& t){return coherentSurv::gamma_evaluate(t,id,dt,X,y,4,2);};
 const auto fit=evaluate(theta);std::vector<double> score(3,0.0),numeric(3),hessian(9);
 for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<3;++j)score[j]+=fit.scores[i*3+j];
 const double step=0.0002;
 for(std::size_t j=0;j<3;++j) {
  auto a=theta,b=theta,c=theta,d=theta;a[j]-=2*step;b[j]-=step;c[j]+=step;d[j]+=2*step;
  const auto fa=evaluate(a),fb=evaluate(b),fc=evaluate(c),fd=evaluate(d);
  numeric[j]=(fa.loglik-8*fb.loglik+8*fc.loglik-fd.loglik)/(12*step);
  for(std::size_t k=0;k<3;++k) {double v=0.0;for(std::size_t i=0;i<4;++i)v+=fa.scores[i*3+k]-8*fb.scores[i*3+k]+8*fc.scores[i*3+k]-fd.scores[i*3+k];hessian[k*3+j]=-v/(12*step);}
 }
 const double grad_error=maxdiff(score,numeric),hessian_error=maxdiff(fit.information,hessian);
 require(grad_error<1e-9,"analytic gamma score disagrees with finite difference");
 require(hessian_error<1e-9,"analytic observed information disagrees with finite difference");
 const auto pois=coherentSurv::gamma_evaluate({-.7,.25},id,dt,X,y,4,2,0.0);
 const auto near=coherentSurv::gamma_evaluate({-.7,.25},id,dt,X,y,4,2,1e-12);
 require(std::abs(near.loglik-pois.loglik)<1e-10,"Poisson limit is unstable");
 const std::vector<std::size_t> owner={0,0,1,1};const std::vector<double> lw={std::log(.2),std::log(.3),std::log(.5),std::log(.7)};
 const auto mix=coherentSurv::coarsen_paths(fit.subject_loglik,fit.scores,fit.information_individual,owner,lw,2,3);
 require(std::abs(mix.posterior[0]+mix.posterior[1]-1.0)<1e-14,"coarsening posterior not normalized");
 double coarsened_error=0.0;
 for(std::size_t j=0;j<3;++j) {
  auto a=theta,b=theta,c=theta,d=theta;a[j]-=2*step;b[j]-=step;c[j]+=step;d[j]+=2*step;
  auto mixed=[&](const std::vector<double>& t){const auto f=evaluate(t);return coherentSurv::coarsen_paths(f.subject_loglik,f.scores,f.information_individual,owner,lw,2,3).likelihood;};
  const auto fa=mixed(a),fb=mixed(b),fc=mixed(c),fd=mixed(d);
  for(std::size_t k=0;k<3;++k) {double v=0;for(std::size_t i=0;i<2;++i)v+=fa.scores[i*3+k]-8*fb.scores[i*3+k]+8*fc.scores[i*3+k]-fd.scores[i*3+k];coarsened_error=std::max(coarsened_error,std::abs(-v/(12*step)-mix.likelihood.information[k*3+j]));}
 }
 require(coarsened_error<1e-9,"Louis information disagrees with derivative");
 const std::vector<double> Q={-.5,.3,.2,0,-.4,.4,0,0,0};
 const std::vector<double> reward={.3,1.,0,1.,0,0};
 const auto forecast=coherentSurv::uniformize({1,0,0},Q,2,reward,2);
 const double p0=std::exp(-1.),p1=3*(std::exp(-.8)-std::exp(-1.));
 const double probability_error=maxdiff(forecast.probability,{p0,p1,1-p0-p1});
 const double alive0=-std::expm1(-1.)/.5,alive1=3*(-std::expm1(-.8)/.4-alive0);
 const double reward_error=maxdiff(forecast.reward,{.3*alive0,alive0+alive1});
 require(probability_error<1e-12 && reward_error<1e-12,"uniformization disagrees with closed form");
 const auto half=coherentSurv::uniformize({1,0,0},Q,1,reward,2);
 const auto second=coherentSurv::uniformize(half.probability,Q,1,reward,2);
 const double tower=maxdiff(second.probability,forecast.probability);require(tower<1e-12,"ordered semigroup mismatch");
 const double big=std::exp(30.);
 const auto extreme=coherentSurv::uniformize({1,0},{-big,big,0,0},1,{1,0},1);
 for(double v:extreme.probability)require(std::isfinite(v),"nonfinite stress probability");
 require(std::isfinite(extreme.reward[0]),"nonfinite stress reward");
 const auto zero=coherentSurv::uniformize({1},{0},3,{2},1);require(zero.reward[0]==6.,"zero generator reward");
 bool invalid=false;try {coherentSurv::uniformize({1,0},{-.5,.4,0,0},1,{},0);}catch(const std::invalid_argument&){invalid=true;}
 require(invalid,"invalid generator was accepted");
 std::cout<<std::setprecision(17)<<"{\n  \"gradient_error\": "<<grad_error<<",\n  \"information_error\": "<<hessian_error<<",\n  \"Louis_information_error\": "<<coarsened_error<<",\n  \"uniformization_error\": "<<probability_error<<",\n  \"reward_error\": "<<reward_error<<",\n  \"tower_error\": "<<tower<<",\n  \"loglik\": "<<fit.loglik<<",\n  \"stress_terminal_probability\": "<<extreme.probability[1]<<"\n}\n";
 return 0;
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
