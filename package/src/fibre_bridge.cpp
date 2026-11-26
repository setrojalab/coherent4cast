#include <Rcpp.h>
#include <libcoherentSurv/fibre.hpp>
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::export]]
List cpp_cs_coarsen_fibre(NumericVector lower, NumericVector upper,
  IntegerVector mark, LogicalVector terminal, double followup,
  NumericVector breaks, NumericVector nodes, NumericVector weights,
  int max_perms, int draws, double seed, int max_paths) {
  if(lower.size()!=upper.size()||lower.size()!=mark.size()||lower.size()!=terminal.size()||
    max_perms<1||draws<2||max_paths<1||!R_finite(seed)||seed<0)stop("Invalid fibre inputs");
  std::vector<coherentSurv::FibreEvent> events;
  for(R_xlen_t j=0;j<lower.size();++j){if(mark[j]==NA_INTEGER||terminal[j]==NA_LOGICAL)stop("Missing event mark or terminal indicator");
    events.push_back({lower[j],upper[j],mark[j],static_cast<bool>(terminal[j])});}
  const auto out=coherentSurv::coarsening_fibre(events,followup,as<std::vector<double>>(breaks),
    as<std::vector<double>>(nodes),as<std::vector<double>>(weights),
    static_cast<std::size_t>(max_perms),static_cast<std::size_t>(draws),
    static_cast<std::uint64_t>(seed),static_cast<std::size_t>(max_paths));
  const int n=static_cast<int>(out.paths.size()),m=static_cast<int>(events.size());
  NumericMatrix times(n,m);IntegerMatrix orders(n,m);NumericVector log_weights(n);IntegerVector draw(n);
  for(int i=0;i<n;++i){log_weights[i]=out.paths[static_cast<std::size_t>(i)].log_weight;draw[i]=out.paths[static_cast<std::size_t>(i)].draw;
    for(int j=0;j<m;++j){times(i,j)=out.paths[static_cast<std::size_t>(i)].time[static_cast<std::size_t>(j)];orders(i,j)=out.paths[static_cast<std::size_t>(i)].order[static_cast<std::size_t>(j)]+1;}}
  return List::create(_["times"]=times,_["order"]=orders,_["log_weights"]=log_weights,_["draw"]=draw,
    _["method"]=out.method,_["reason"]=out.reason,_["permutations"]=static_cast<double>(out.permutations),
    _["attempted"]=static_cast<double>(out.attempted),_["rejected"]=static_cast<double>(out.rejected),
    _["volume"]=out.volume,_["volume_mcse"]=out.volume_mcse,_["importance_ess"]=out.importance_ess);
}
