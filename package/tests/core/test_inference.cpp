#include <libcoherentSurv/target.hpp>
#include <libcoherentSurv/lp_bounds.hpp>
#include <libcoherentSurv/conformal.hpp>
#include <libcoherentSurv/pivot.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
int main(){
 const std::vector<double> p{.5,1.0/3,1.0/6},c{0,0,1},d{0,1,1};
 for(double k:{1.0,1.5,2.0,3.0}){std::vector<double> r=p;for(double& v:r){v*=.7;}auto b=coherentSurv::compatible_bounds(p,c,d,k,.7,r,r);std::cout<<k<<" "<<b.lower<<" "<<b.upper<<" "<<b.gap<<"\n";require(b.gap<1e-12,"LP dual gap");require(b.lower<=b.upper,"LP order");if(k==1){require(std::abs(b.lower-1.0/3)<1e-12,"MAR bound");}}
 auto co=coherentSurv::compatible_bounds({.5,.5},{0,1,0},{1,1,1},1,.7,{.35,.35},{.35,.35},{0,1,1});require(std::abs(co.lower)<1e-12&&std::abs(co.upper-.5)<1e-12,"fibre bounds");
 const auto ranks=coherentSurv::weighted_rank({.1,.5,.9},{.7},{1,1,1},{2,2,2},{2});require(std::abs(ranks[0]-2.0/3)<1e-12,"conformal corner");
 const auto extreme=coherentSurv::weighted_rank({.1,.5,.9},{.7},{1e308,5e307,1e308},{1e308,5e307,1e308},{1e308});require(std::abs(extreme[0]-4.0/7)<1e-15,"finite extreme weights retain scale invariance");
 const auto inactive=coherentSurv::weighted_rank({.1,.5,.9},{.7},{1e-100,1e-100,1e-100},{1e308,1e308,1e-100},{1e-100});require(std::abs(inactive[0]-.5)<1e-15,"inactive upper weights do not distort normalization");
 const double rate=coherentSurv::invert_pivot({1,2,3},0,6);require(std::abs(rate-1)<1e-12,"exponential pivot");
 auto t=coherentSurv::target_fit({{1,1},{0,2}},{{.5,.5},{.5,.5}},{{0},{1}},{.8,.8},4);
 require(t.converged,"target convergence");require(std::abs(t.q[0][1]-.75)<1e-6,"target solution");
 std::cout<<"inference core checks passed\n";
}
