#include <Rcpp.h>
#include <libcoherent4cast/gamma_lik.hpp>
#include <libcoherent4cast/coarsen.hpp>
#include <libcoherent4cast/uniformize.hpp>
using namespace Rcpp;
namespace {
std::vector<double> row_major(NumericMatrix x) {
 std::vector<double> out(static_cast<std::size_t>(x.nrow()*x.ncol()));
 for(int i=0;i<x.nrow();++i)for(int j=0;j<x.ncol();++j)out[static_cast<std::size_t>(i*x.ncol()+j)]=x(i,j);
 return out;
}
std::vector<std::size_t> indices(IntegerVector x) {
 std::vector<std::size_t> out;out.reserve(static_cast<std::size_t>(x.size()));
 for(int value:x){if(value==NA_INTEGER || value<1)stop("unit indices must be positive integers");out.push_back(static_cast<std::size_t>(value-1));}return out;
}
NumericMatrix to_matrix(const std::vector<double>& values,std::size_t nr,std::size_t nc) {
 NumericMatrix out(static_cast<int>(nr),static_cast<int>(nc));
 for(std::size_t i=0;i<nr;++i)for(std::size_t j=0;j<nc;++j) {out(static_cast<int>(i),static_cast<int>(j))=values[i*nc+j];}
 return out;
}
List likelihood_list(const coherent4cast::GammaResult& x) {
 const std::size_t d=x.dimension,n=x.n;NumericVector g(static_cast<int>(d));
 for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<d;++j)g[static_cast<int>(j)]-=x.scores[i*d+j]/static_cast<double>(n);
 List out=List::create(_["loglik"]=x.loglik,_["objective"]=-x.loglik/static_cast<double>(n),_["gradient"]=g,
 _["subject_loglik"]=x.subject_loglik,_["scores"]=to_matrix(x.scores,n,d),_["counts"]=x.counts,_["cumulative_hazard"]=x.cumulative_hazard);
 if(!x.information.empty()){
  out["information"]=to_matrix(x.information,d,d);NumericVector a(static_cast<int>(n*d*d));a.attr("dim")=IntegerVector::create(static_cast<int>(d),static_cast<int>(d),static_cast<int>(n));
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<d;++j)for(std::size_t k=0;k<d;++k) {a[static_cast<int>(j+d*k+d*d*i)]=x.information_individual[(i*d+j)*d+k];}
  out["information_individual"]=a;
 }return out;
}
}
// [[Rcpp::export]]
List cpp_rc_loglik(NumericVector theta,IntegerVector unit,NumericVector exposure,NumericMatrix design,NumericVector event,double phi_fixed=-1.0,bool information=true) {
 const auto id=indices(unit);if(id.empty())stop("no observations");const std::size_t n=*std::max_element(id.begin(),id.end())+1U;
 return likelihood_list(coherent4cast::gamma_evaluate(as<std::vector<double>>(theta),id,as<std::vector<double>>(exposure),row_major(design),as<std::vector<double>>(event),n,static_cast<std::size_t>(design.ncol()),phi_fixed,information));
}
// [[Rcpp::export]]
List cpp_rc_coarsen(NumericVector loglik,NumericMatrix scores,NumericVector information,IntegerVector owner,NumericVector log_weights) {
 const auto id=indices(owner);if(id.empty())stop("no paths");const std::size_t n=*std::max_element(id.begin(),id.end())+1U,d=static_cast<std::size_t>(scores.ncol());std::vector<double> info(static_cast<std::size_t>(information.size()));
 if(!info.empty()){if(info.size()!=id.size()*d*d)stop("incorrect path information dimensions");for(std::size_t i=0;i<id.size();++i)for(std::size_t j=0;j<d;++j)for(std::size_t k=0;k<d;++k)info[(i*d+j)*d+k]=information[static_cast<int>(j+d*k+d*d*i)];}
 const auto result=coherent4cast::coarsen_paths(as<std::vector<double>>(loglik),row_major(scores),info,id,as<std::vector<double>>(log_weights),n,d);List out=likelihood_list(result.likelihood);out["posterior"]=result.posterior;return out;
}
// [[Rcpp::export]]
List cpp_rc_propagate(NumericVector probability,NumericMatrix generator,double duration,NumericMatrix reward,double tolerance=1e-12,int max_terms=10000) {
 if(max_terms<2)stop("max_terms must be at least two");
 const auto result=coherent4cast::uniformize(as<std::vector<double>>(probability),row_major(generator),duration,row_major(reward),static_cast<std::size_t>(reward.ncol()),tolerance,static_cast<std::size_t>(max_terms));
 return List::create(_["probability"]=result.probability,_["reward"]=result.reward,_["error_bound"]=result.error_bound,_["reward_error_bound"]=result.reward_error_bound,_["roundoff_bound"]=result.roundoff_bound,_["terms"]=static_cast<double>(result.terms),_["squarings"]=static_cast<double>(result.squarings));
}
