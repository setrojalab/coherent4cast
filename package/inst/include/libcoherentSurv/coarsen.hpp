#ifndef LIBcoherentSurv_COARSEN_HPP
#define LIBcoherentSurv_COARSEN_HPP
#include "gamma_lik.hpp"
namespace coherentSurv {
struct CoarsenedResult {GammaResult likelihood;std::vector<double> posterior;};
inline CoarsenedResult coarsen_paths(const std::vector<double>& loglik,
 const std::vector<double>& scores,const std::vector<double>& path_information,
 const std::vector<std::size_t>& owner,const std::vector<double>& log_weights,
 std::size_t n,std::size_t d) {
 const std::size_t paths=owner.size();const bool info=!path_information.empty();
 if (paths==0 || n==0 || d==0 || loglik.size()!=paths || log_weights.size()!=paths || scores.size()!=paths*d ||
     (info && path_information.size()!=paths*d*d)) throw std::invalid_argument("coarsen: inconsistent dimensions");
 CoarsenedResult result;auto& out=result.likelihood;out.n=n;out.dimension=d;
 out.subject_loglik.assign(n,0.0);out.scores.assign(n*d,0.0);result.posterior.resize(paths);
 std::vector<long double> maxima(n,-std::numeric_limits<long double>::infinity()),mass(n,0.0L);
 std::vector<std::size_t> number(n,0U);
 for (std::size_t h=0;h<paths;++h) {
  if (owner[h]>=n || !std::isfinite(loglik[h]) || !std::isfinite(log_weights[h])) throw std::invalid_argument("coarsen: invalid path or quadrature log weight");
  const long double z=static_cast<long double>(loglik[h])+log_weights[h];
  maxima[owner[h]]=std::max(maxima[owner[h]],z);++number[owner[h]];
 }
 if (std::find(number.begin(),number.end(),0U)!=number.end()) throw std::invalid_argument("coarsen: an observed unit has no legal path");
 for (std::size_t h=0;h<paths;++h) mass[owner[h]]+=std::exp(static_cast<long double>(loglik[h])+log_weights[h]-maxima[owner[h]]);
 for (std::size_t i=0;i<n;++i) {out.subject_loglik[i]=finite_double(maxima[i]+std::log(mass[i]),"coarsen: log likelihood overflow");out.loglik+=out.subject_loglik[i];}
 if (info) {out.information.assign(d*d,0.0);out.information_individual.assign(n*d*d,0.0);}
 for (std::size_t h=0;h<paths;++h) {
  const std::size_t i=owner[h];const double w=static_cast<double>(std::exp(static_cast<long double>(loglik[h])+log_weights[h]-maxima[i])/mass[i]);result.posterior[h]=w;
  for (std::size_t j=0;j<d;++j) {
   if (!std::isfinite(scores[h*d+j])) throw std::invalid_argument("coarsen: nonfinite path score");
   out.scores[i*d+j]+=w*scores[h*d+j];
   if (info) for (std::size_t k=0;k<d;++k) {
    const double v=path_information[(h*d+j)*d+k];
    if (!std::isfinite(v)) throw std::invalid_argument("coarsen: nonfinite path information");
    out.information_individual[(i*d+j)*d+k]+=w*(v-scores[h*d+j]*scores[h*d+k]);
   }
  }
 }
 if (info) for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<d;++j) for (std::size_t k=0;k<d;++k) {
  const std::size_t cell=(i*d+j)*d+k;
  out.information_individual[cell]+=out.scores[i*d+j]*out.scores[i*d+k];
  out.information[j*d+k]+=out.information_individual[cell];
 }
 for(double value:out.scores)if(!std::isfinite(value))throw std::overflow_error("coarsen: score overflow");
 for(double value:out.information_individual)if(!std::isfinite(value))throw std::overflow_error("coarsen: individual information overflow");
 for(double value:out.information)if(!std::isfinite(value))throw std::overflow_error("coarsen: summed information overflow");
 if(!std::isfinite(out.loglik))throw std::overflow_error("coarsen: summed log likelihood overflow");
 return result;
}
}
#endif
