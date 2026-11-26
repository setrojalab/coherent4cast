#ifndef LIBcoherentSurv_HISTORY_HPP
#define LIBcoherentSurv_HISTORY_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace coherentSurv {
struct HistoryResult {std::vector<double> counts,recency,gap;};
inline HistoryResult compile_history(const std::vector<int>& unit,
 const std::vector<double>& start,const std::vector<double>& stop,
 const std::vector<int>& mark,std::size_t types,double cap,double decay) {
 const std::size_t m=unit.size();
 if(start.size()!=m||stop.size()!=m||mark.size()!=m||cap<0||!std::isfinite(cap)||!std::isfinite(decay))
  throw std::invalid_argument("history: inconsistent arguments");
 HistoryResult out;out.counts.resize(m*types);out.recency.resize(m*types);out.gap.resize(m*types);
 std::vector<double> count(types),recent(types),last(types,std::numeric_limits<double>::quiet_NaN());
 double previous=0;int current=0;
 for(std::size_t i=0;i<m;++i){
  if(unit[i]<1||!std::isfinite(start[i])||!std::isfinite(stop[i])||start[i]<0||stop[i]<start[i]||mark[i]<0||static_cast<std::size_t>(mark[i])>types)
   throw std::invalid_argument("history: invalid row");
  if(unit[i]!=current){if(unit[i]<current)throw std::invalid_argument("history: unsorted units");current=unit[i];std::fill(count.begin(),count.end(),0);std::fill(recent.begin(),recent.end(),0);std::fill(last.begin(),last.end(),std::numeric_limits<double>::quiet_NaN());previous=start[i];}
  if(start[i]<previous-1e-12)throw std::invalid_argument("history: overlapping intervals");
  for(std::size_t k=0;k<types;++k){
   if(decay>=0)recent[k]*=std::exp(-decay*(start[i]-previous));
   out.counts[i*types+k]=std::min(count[k],cap);out.recency[i*types+k]=recent[k];out.gap[i*types+k]=start[i]-last[k];
   if(decay>=0)recent[k]*=std::exp(-decay*(stop[i]-start[i]));
  }
  if(mark[i]>0){const std::size_t k=static_cast<std::size_t>(mark[i]-1);count[k]+=1;last[k]=stop[i];recent[k]+=1;}
  previous=stop[i];
 }
 return out;
}
}
#endif
