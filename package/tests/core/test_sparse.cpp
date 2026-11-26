#include <libcoherentSurv/sparse_uniformize.hpp>
#include <libcoherentSurv/sparse_bounds.hpp>
#include <libcoherentSurv/uniformize.hpp>
#include <libcoherentSurv/lp_bounds.hpp>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <random>
int main(int argc,char** argv) {
 using namespace coherentSurv;std::cout<<std::setprecision(17);
 CSRGenerator q{3,{0,2,3,3},{1,2,2},{0.3,0.2,0.4}};
 const std::vector<double> dense{-0.5,0.3,0.2,0,-0.4,0.4,0,0,0},p{1,0,0},reward{1,0.3,1,0.5,0,0};
 const auto s=sparse_uniformize(p,q,2.7,reward,2,1e-13);
 const auto d=uniformize(p,dense,2.7,reward,2,1e-13);
 double pe=0,re=0;
 for(std::size_t i=0;i<3;++i)pe=std::max(pe,std::abs(s.probability[i]-d.probability[i]));
 for(std::size_t i=0;i<2;++i)re=std::max(re,std::abs(s.reward[i]-d.reward[i]));
 assert(pe<1e-12 && re<1e-12);
 const auto s1=sparse_uniformize(p,q,40,reward,2,1e-13);
 const auto d1=uniformize(p,dense,40,reward,2,1e-13);
 for(std::size_t i=0;i<3;++i)assert(std::abs(s1.probability[i]-d1.probability[i])<1e-12);
 for(std::size_t i=0;i<2;++i)assert(std::abs(s1.reward[i]-d1.reward[i])<1e-11);
 std::mt19937 rng(20261122);std::uniform_real_distribution<double> unif(0.1,1);
 double max_parity=0,max_gap=0;std::size_t dense_failed=0,dense_passed=0;
 std::ofstream fixtures;if(argc>1){fixtures.open(argv[1]);fixtures<<std::setprecision(17)<<"rep,atom,fibre,p,lo,hi,num,den,rho,kappa,lower,upper,gap\n";}
 for(std::size_t rep=0;rep<80;++rep) {
  const std::size_t j=6,m=3;const double rho=0.35+0.4*unif(rng),kappa=rep%4==0?1:1+3*unif(rng);
  std::vector<double> law(m),a(j),b(j),lo(m),hi(m);std::vector<std::size_t> map{0,0,1,1,2,2};
  double mass=0;for(double& x:law){x=unif(rng);mass+=x;}for(double& x:law)x/=mass;
  for(std::size_t k=0;k<j;++k){a[k]=2*unif(rng)-1;b[k]=unif(rng);}
  for(std::size_t h=0;h<m;++h){lo[h]=rho*law[h]*(rep%2==0?1:0.6);hi[h]=std::min(rho,rho*law[h]*(rep%2==0?1:1.6));}
  const auto fast=compatible_bounds_sparse(law,a,b,kappa,rho,lo,hi,map,1e-12);
  for(std::size_t k=0;k<j;++k)if(fixtures)fixtures<<rep<<","<<k<<","<<map[k]<<","<<law[map[k]]<<","<<lo[map[k]]<<","<<hi[map[k]]<<","<<a[k]<<","<<b[k]<<","<<rho<<","<<kappa<<","<<fast.lower.value<<","<<fast.upper.value<<","<<std::max(fast.lower.dual_gap,fast.upper.dual_gap)<<"\n";
  BoundResult oracle;
  try {oracle=compatible_bounds(law,a,b,kappa,rho,lo,hi,map);}
  catch(const std::exception&){++dense_failed;continue;}
  ++dense_passed;
  max_parity=std::max(max_parity,std::max(std::abs(fast.lower.value-oracle.lower),std::abs(fast.upper.value-oracle.upper)));
  max_gap=std::max(max_gap,std::max(fast.lower.dual_gap,fast.upper.dual_gap));
  if(max_parity>1e-8){std::cerr<<"Failed replicate "<<rep<<" "<<max_parity<<"\n";return 1;}
 }
 // Repeated positive/negative values, exact ties and 50,000 distinct atoms.
 const std::size_t j=50000,m=1000;std::vector<double> law(m,1.0/static_cast<double>(m)),a(j),b(j),lo(m),hi(m);
 std::vector<std::size_t> map(j);
 for(std::size_t k=0;k<j;++k){map[k]=k%m;b[k]=0.2+0.8*static_cast<double>((k*104729)%50023)/50023;a[k]=b[k]*(0.1+0.8*static_cast<double>((k*7919)%50021)/50021);}
 for(std::size_t h=0;h<m;++h){lo[h]=0.6*law[h]*0.8;hi[h]=0.6*law[h]*1.2;}
 auto begin=std::chrono::steady_clock::now();const auto big=compatible_bounds_sparse(law,a,b,3,0.6,lo,hi,map,1e-11);
 const double bounds_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
 assert(big.lower.dual_gap<1e-9 && big.upper.dual_gap<1e-9);
 CSRGenerator chain;chain.n=j;chain.offsets.push_back(0);
 for(std::size_t k=0;k<j;++k){if(k+1<j){chain.columns.push_back(k+1);chain.rates.push_back(0.1);}chain.offsets.push_back(chain.rates.size());}
 std::vector<double> initial(j,0),rw(j,1);initial[0]=1;
 begin=std::chrono::steady_clock::now();const auto large=sparse_uniformize(initial,chain,10,rw,1,1e-12);
 const double propagation_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
 double total=0;for(double z:large.probability)total+=z;
 assert(std::abs(total-1)<1e-12);assert(std::abs(large.probability[0]-std::exp(-1.0))<1e-12);
 assert(std::abs(large.reward[0]-10)<large.reward_error_bound[0]+1e-12);
 std::cout<<"{\n  \"small_probability_error\": "<<pe<<",\n  \"small_reward_error\": "<<re
 <<",\n  \"random_dense_lp_parity_max\": "<<max_parity<<",\n  \"random_fractional_certificate_gap_max\": "<<max_gap
 <<",\n  \"dense_oracle_passed\": "<<dense_passed<<",\n  \"dense_oracle_failed\": "<<dense_failed
 <<",\n  \"support_atoms\": "<<j<<",\n  \"coarsening_fibres\": "<<m
 <<",\n  \"bounds_seconds\": "<<bounds_seconds<<",\n  \"bounds_workspace_bytes_estimate\": "<<big.workspace_bytes
 <<",\n  \"bounds_lower\": "<<big.lower.value<<",\n  \"bounds_upper\": "<<big.upper.value
 <<",\n  \"bounds_gap\": "<<std::max(big.lower.dual_gap,big.upper.dual_gap)
 <<",\n  \"csr_states\": "<<j<<",\n  \"csr_transitions\": "<<chain.rates.size()
 <<",\n  \"csr_seconds\": "<<propagation_seconds<<",\n  \"csr_workspace_bytes_estimate\": "<<large.workspace_bytes
 <<",\n  \"csr_terms\": "<<large.terms<<",\n  \"csr_probability_error_bound\": "<<large.error_bound
 <<",\n  \"csr_reward_error_bound\": "<<large.reward_error_bound[0]<<"\n}\n";
}
