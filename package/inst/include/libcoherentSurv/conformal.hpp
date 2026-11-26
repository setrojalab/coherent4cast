#ifndef LIBcoherentSurv_CONFORMAL_HPP
#define LIBcoherentSurv_CONFORMAL_HPP
#include <cmath>
#include <stdexcept>
#include <vector>
namespace coherentSurv {
inline std::vector<double> weighted_rank(const std::vector<double>& scores,const std::vector<double>& candidates,const std::vector<double>& lower,const std::vector<double>& upper,const std::vector<double>& candidate_upper){
 if(scores.empty()||lower.size()!=scores.size()||upper.size()!=scores.size()||candidate_upper.size()!=candidates.size())throw std::invalid_argument("weighted rank dimension mismatch");
 for(std::size_t i=0;i<scores.size();++i)if(!std::isfinite(scores[i])||!std::isfinite(lower[i])||lower[i]<0||std::isnan(upper[i])||upper[i]<lower[i])throw std::invalid_argument("invalid conformal score or weight interval");
 std::vector<double> out(candidates.size());
 for(std::size_t k=0;k<candidates.size();++k){
  if(!std::isfinite(candidates[k])||std::isnan(candidate_upper[k])||candidate_upper[k]<=0)throw std::invalid_argument("candidate scores must be finite and candidate weights positive");
  // Normalize only masses participating in this candidate's worst-case rank.
  // This preserves scale invariance without overflow, including irrelevant
  // infinite upper bounds on observations that are in the rank's head.
  double scale=candidate_upper[k];
  for(std::size_t i=0;i<scores.size();++i){const double w=scores[i]>=candidates[k]?upper[i]:lower[i];if(w>scale)scale=w;}
  if(std::isinf(scale)){out[k]=1;continue;}
  double tail=candidate_upper[k]/scale,head=0;
  for(std::size_t i=0;i<scores.size();++i){if(scores[i]>=candidates[k])tail+=upper[i]/scale;else head+=lower[i]/scale;}
  if(!std::isfinite(tail+head)||tail+head<=0)throw std::runtime_error("invalid normalized conformal weight mass");
  out[k]=tail/(tail+head);
 }return out;
}
}
#endif
