#include <Rcpp.h>
#include <libcoherent4cast/target.hpp>
#include <libcoherent4cast/lp_bounds.hpp>
#include <libcoherent4cast/conformal.hpp>
#include <libcoherent4cast/pivot.hpp>
namespace {
coherent4cast::TargetMatrix to_target_matrix(const Rcpp::NumericMatrix& a){coherent4cast::TargetMatrix z(static_cast<std::size_t>(a.nrow()),std::vector<double>(static_cast<std::size_t>(a.ncol())));for(int i=0;i<a.nrow();++i)for(int j=0;j<a.ncol();++j)z[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]=a(i,j);return z;}
Rcpp::NumericMatrix from_target_matrix(const coherent4cast::TargetMatrix& a){const int n=static_cast<int>(a.size()),p=n?static_cast<int>(a[0].size()):0;Rcpp::NumericMatrix z(n,p);for(int i=0;i<n;++i)for(int j=0;j<p;++j)z(i,j)=a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];return z;}
}
// [[Rcpp::export]]
Rcpp::List rc_cpp_target(Rcpp::NumericMatrix counts,Rcpp::NumericMatrix initial,Rcpp::NumericMatrix features,Rcpp::NumericVector recording,double n,int maxit,double tolerance){const auto a=coherent4cast::target_fit(to_target_matrix(counts),to_target_matrix(initial),to_target_matrix(features),Rcpp::as<std::vector<double>>(recording),n,maxit,tolerance);return Rcpp::List::create(Rcpp::Named("q")=from_target_matrix(a.q),Rcpp::Named("mu")=from_target_matrix(a.mu),Rcpp::Named("epsilon")=a.epsilon,Rcpp::Named("gradient")=a.gradient,Rcpp::Named("iterations")=a.iterations,Rcpp::Named("converged")=a.converged);}
// [[Rcpp::export]]
Rcpp::List rc_cpp_bounds(Rcpp::NumericVector p,Rcpp::NumericVector numerator,Rcpp::NumericVector denominator,double kappa,double rho,Rcpp::NumericVector lower,Rcpp::NumericVector upper,Rcpp::IntegerVector observation_map){std::vector<std::size_t> map;for(int x:observation_map){if(x<0)Rcpp::stop("observation map is negative");map.push_back(static_cast<std::size_t>(x));}const auto a=coherent4cast::compatible_bounds(Rcpp::as<std::vector<double>>(p),Rcpp::as<std::vector<double>>(numerator),Rcpp::as<std::vector<double>>(denominator),kappa,rho,Rcpp::as<std::vector<double>>(lower),Rcpp::as<std::vector<double>>(upper),map);return Rcpp::List::create(Rcpp::Named("lower")=a.lower,Rcpp::Named("upper")=a.upper,Rcpp::Named("primal_dual_gap")=a.gap,Rcpp::Named("primal_residual")=a.residual,Rcpp::Named("lower_law")=a.lower_law,Rcpp::Named("upper_law")=a.upper_law);}
// [[Rcpp::export]]
Rcpp::NumericVector rc_cpp_conformal(Rcpp::NumericVector scores,Rcpp::NumericVector candidates,Rcpp::NumericVector lower,Rcpp::NumericVector upper,Rcpp::NumericVector candidate_upper){return Rcpp::wrap(coherent4cast::weighted_rank(Rcpp::as<std::vector<double>>(scores),Rcpp::as<std::vector<double>>(candidates),Rcpp::as<std::vector<double>>(lower),Rcpp::as<std::vector<double>>(upper),Rcpp::as<std::vector<double>>(candidate_upper)));}
// [[Rcpp::export]]
Rcpp::NumericVector rc_cpp_pivot(Rcpp::NumericVector times,double dispersion,Rcpp::NumericVector quantiles){Rcpp::NumericVector out(quantiles.size());const auto t=Rcpp::as<std::vector<double>>(times);for(R_xlen_t i=0;i<quantiles.size();++i)out[i]=coherent4cast::invert_pivot(t,dispersion,quantiles[i]);return out;}
// [[Rcpp::export]]
Rcpp::NumericVector rc_cpp_cf(Rcpp::NumericVector z,double n,double skewness,double kurtosis,int order){Rcpp::NumericVector out(z.size());for(R_xlen_t i=0;i<z.size();++i)out[i]=coherent4cast::cf_quantile(z[i],n,skewness,kurtosis,order);return out;}
