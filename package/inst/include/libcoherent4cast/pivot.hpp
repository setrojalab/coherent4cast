#ifndef LIBcoherent4cast_PIVOT_HPP
#define LIBcoherent4cast_PIVOT_HPP
#include <cmath>
#include <stdexcept>
#include <vector>
namespace coherent4cast {
inline double softplus(double x){return x>0?x+std::log1p(std::exp(-x)):std::log1p(std::exp(x));}
inline double pivot_sum(const std::vector<double>& times,double phi,double log_lambda){
 double value=0;for(double t:times){const double z=log_lambda+std::log(t);value+=phi==0?std::exp(z):softplus(std::log(phi)+z)/phi;}return value;
}
inline double invert_pivot(const std::vector<double>& times,double phi,double quantile){
 if(times.empty()||!std::isfinite(phi)||phi<0||!std::isfinite(quantile)||quantile<=0)throw std::invalid_argument("invalid pivot inputs");
 for(double t:times)if(!std::isfinite(t)||t<=0)throw std::invalid_argument("pivot requires positive uncensored terminal times");
 double lower=-700,upper=700;
 if(pivot_sum(times,phi,lower)>quantile||pivot_sum(times,phi,upper)<quantile)throw std::runtime_error("pivot root outside representable positive rate range");
 for(int it=0;it<200;++it){const double middle=(lower+upper)/2;if(pivot_sum(times,phi,middle)>quantile)upper=middle;else lower=middle;if(upper-lower<1e-13)break;}
 return std::exp((lower+upper)/2);
}
inline double cf_quantile(double z,double n,double skewness,double kurtosis,int order){
 if(!std::isfinite(z)||!std::isfinite(n)||n<=0||!std::isfinite(skewness)||!std::isfinite(kurtosis)||order<0||order>2)throw std::invalid_argument("invalid Cornish-Fisher arguments");
 double q=z;if(order>0)q+=skewness*(z*z-1)/(6*std::sqrt(n));if(order>1)q+=kurtosis*(z*z*z-3*z)/(24*n)-skewness*skewness*(2*z*z*z-5*z)/(36*n);return q;
}
}
#endif
