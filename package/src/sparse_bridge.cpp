#include <Rcpp.h>
#include <librecurcast/sparse_uniformize.hpp>
#include <librecurcast/sparse_bounds.hpp>
using namespace Rcpp;
namespace {
std::vector<std::size_t> sparse_indices(IntegerVector x) {
 std::vector<std::size_t> y; y.reserve(static_cast<std::size_t>(x.size()));
 for(int v:x){if(v==NA_INTEGER || v<0)stop("CSR indices must be nonnegative integers");y.push_back(static_cast<std::size_t>(v));}return y;
}
List bound_endpoint(const recurcast::SparseBoundEndpoint& e) {
 return List::create(_["value"]=e.value,_["certificate_lower"]=e.certificate_lower,
 _["certificate_upper"]=e.certificate_upper,_["dual_gap"]=e.dual_gap,_["residual"]=e.residual,
 _["linear_primal"]=e.linear_primal,_["linear_dual"]=e.linear_dual,_["dual_mass"]=e.dual_mass,
 _["dual_high_mass"]=e.dual_high_mass,_["denominator"]=e.denominator,
 _["linearization_point"]=e.linearization_point,_["objective_sign"]=e.objective_sign,
 _["iterations"]=static_cast<double>(e.iterations),_["law"]=e.law,_["recorded"]=e.recorded,
 _["unrecorded"]=e.unrecorded,_["fibre_dual"]=e.fibre_dual);
}
}
// [[Rcpp::export]]
List cpp_rc_propagate_sparse(NumericVector probability,IntegerVector offsets,IntegerVector columns,
 NumericVector rates,double duration,NumericMatrix reward,double tolerance=1e-12,
 int max_terms=10000,int max_chunks=100000) {
 if(max_terms<2 || max_chunks<1)stop("invalid sparse propagation budget");
 recurcast::CSRGenerator q;q.n=static_cast<std::size_t>(probability.size());q.offsets=sparse_indices(offsets);
 q.columns=sparse_indices(columns);q.rates=as<std::vector<double>>(rates);
 std::vector<double> r(static_cast<std::size_t>(reward.nrow()*reward.ncol()));
 for(int i=0;i<reward.nrow();++i)for(int j=0;j<reward.ncol();++j)r[static_cast<std::size_t>(i*reward.ncol()+j)]=reward(i,j);
 const auto z=recurcast::sparse_uniformize(as<std::vector<double>>(probability),q,duration,r,
 static_cast<std::size_t>(reward.ncol()),tolerance,static_cast<std::size_t>(max_terms),static_cast<std::size_t>(max_chunks));
 return List::create(_["probability"]=z.probability,_["reward"]=z.reward,
 _["error_bound"]=z.error_bound,_["reward_error_bound"]=z.reward_error_bound,
 _["roundoff_bound"]=z.roundoff_bound,_["terms"]=static_cast<double>(z.terms),
 _["chunks"]=static_cast<double>(z.chunks),_["workspace_bytes"]=static_cast<double>(z.workspace_bytes));
}
// [[Rcpp::export]]
List cpp_rc_bounds_sparse(NumericVector p,NumericVector numerator,NumericVector denominator,
 double kappa,double rho,NumericVector lower,NumericVector upper,IntegerVector observation_map,
 double tolerance=1e-10,int max_iterations=1000) {
 if(max_iterations<1)stop("max_iterations must be positive");
 const auto z=recurcast::compatible_bounds_sparse(as<std::vector<double>>(p),as<std::vector<double>>(numerator),
 as<std::vector<double>>(denominator),kappa,rho,as<std::vector<double>>(lower),as<std::vector<double>>(upper),
 sparse_indices(observation_map),tolerance,static_cast<std::size_t>(max_iterations));
 return List::create(_["lower"]=z.lower.certificate_lower,_["upper"]=z.upper.certificate_upper,
 _["lower_feasible"]=z.lower.value,_["upper_feasible"]=z.upper.value,
 _["gap"]=std::max(z.lower.dual_gap,z.upper.dual_gap),_["residual"]=std::max(z.lower.residual,z.upper.residual),
 _["lower_law"]=z.lower.law,_["upper_law"]=z.upper.law,
 _["lower_certificate"]=bound_endpoint(z.lower),_["upper_certificate"]=bound_endpoint(z.upper),
 _["denominator_lower_bound"]=z.denominator_lower_bound,_["workspace_bytes"]=static_cast<double>(z.workspace_bytes));
}
