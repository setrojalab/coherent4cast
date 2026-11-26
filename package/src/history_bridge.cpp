#include <Rcpp.h>
#include "../inst/include/libcoherentSurv/history.hpp"
// [[Rcpp::export]]
Rcpp::List cs_cpp_history(Rcpp::IntegerVector unit,Rcpp::NumericVector start,
 Rcpp::NumericVector stop,Rcpp::IntegerVector mark,int types,double cap,double decay){
 if(types<0)Rcpp::stop("history: negative type count");
 auto z=coherentSurv::compile_history(Rcpp::as<std::vector<int>>(unit),Rcpp::as<std::vector<double>>(start),
 Rcpp::as<std::vector<double>>(stop),Rcpp::as<std::vector<int>>(mark),static_cast<std::size_t>(types),cap,decay);
 Rcpp::NumericMatrix counts(unit.size(),types),recent(unit.size(),types),gap(unit.size(),types);
 for(R_xlen_t i=0;i<unit.size();++i)for(int k=0;k<types;++k){const std::size_t pos=static_cast<std::size_t>(i)*static_cast<std::size_t>(types)+static_cast<std::size_t>(k);counts(i,k)=z.counts[pos];recent(i,k)=z.recency[pos];gap(i,k)=std::isnan(z.gap[pos])?NA_REAL:z.gap[pos];}
 return Rcpp::List::create(Rcpp::_ ["counts"]=counts,Rcpp::_ ["recency"]=recent,Rcpp::_ ["gap"]=gap);
}
