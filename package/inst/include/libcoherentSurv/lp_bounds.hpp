#ifndef LIBcoherentSurv_LP_BOUNDS_HPP
#define LIBcoherentSurv_LP_BOUNDS_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
namespace coherentSurv {
using LPMatrix=std::vector<std::vector<double>>;
struct LPResult {std::vector<double> solution;double objective=0,dual_objective=0,gap=0,residual=0;int iterations=0;};
inline std::vector<double> lp_solve_system(LPMatrix a,std::vector<double> b){
 const std::size_t n=b.size();
 for(std::size_t k=0;k<n;++k){std::size_t p=k;for(std::size_t i=k+1;i<n;++i)if(std::abs(a[i][k])>std::abs(a[p][k]))p=i;
  if(std::abs(a[p][k])<1e-13){throw std::runtime_error("singular LP basis");}std::swap(a[p],a[k]);std::swap(b[p],b[k]);const double v=a[k][k];
  for(std::size_t j=k;j<n;++j){a[k][j]/=v;}b[k]/=v;
  for(std::size_t i=0;i<n;++i)if(i!=k){const double z=a[i][k];for(std::size_t j=k;j<n;++j){a[i][j]-=z*a[k][j];}b[i]-=z*b[k];}
 }return b;
}
// Two-phase primal simplex. The supplied form is max c'x, A x <= b, x >= 0.
// Bland's index rule and explicit phase-I feasibility prevent cycling at ties.
inline LPResult simplex_max(const LPMatrix& A,const std::vector<double>& b,const std::vector<double>& c){
 const std::size_t m=A.size(),n=c.size();if(m==0||n==0||b.size()!=m)throw std::invalid_argument("empty or incompatible LP");
 std::size_t na=0;for(double x:b){if(!std::isfinite(x))throw std::invalid_argument("nonfinite LP right hand side");if(x<0)++na;}
 const std::size_t total=n+m+na;LPMatrix initial(m,std::vector<double>(total));LPMatrix tab(m+1,std::vector<double>(total+1));std::vector<std::size_t> basis(m);std::vector<double> rhs(m);std::size_t next=n+m;
 for(std::size_t i=0;i<m;++i){if(A[i].size()!=n)throw std::invalid_argument("ragged LP matrix");const double sign=b[i]<0?-1:1;rhs[i]=sign*b[i];
  for(std::size_t j=0;j<n;++j){if(!std::isfinite(A[i][j]))throw std::invalid_argument("nonfinite LP coefficient");initial[i][j]=sign*A[i][j];}
  initial[i][n+i]=sign;
  if(b[i]<0){basis[i]=next;initial[i][next]=1;++next;}else basis[i]=n+i;
  for(std::size_t j=0;j<total;++j){tab[i][j]=initial[i][j];}tab[i][total]=rhs[i];
 }
 int iterations=0;
 auto objective=[&](const std::vector<double>& cost){std::fill(tab[m].begin(),tab[m].end(),0);for(std::size_t j=0;j<total;++j){tab[m][j]=-cost[j];}for(std::size_t i=0;i<m;++i){const double cb=cost[basis[i]];for(std::size_t j=0;j<=total;++j){tab[m][j]+=cb*tab[i][j];}}};
 auto pivot=[&](std::size_t row,std::size_t col){const double value=tab[row][col];for(std::size_t j=0;j<=total;++j){tab[row][j]/=value;}for(std::size_t i=0;i<=m;++i)if(i!=row){const double value2=tab[i][col];for(std::size_t j=0;j<=total;++j){tab[i][j]-=value2*tab[row][j];}}basis[row]=col;};
 auto optimize=[&](std::size_t permitted){for(;;){if(++iterations>100000)throw std::runtime_error("LP iteration limit exceeded");std::size_t col=permitted;for(std::size_t j=0;j<permitted;++j)if(tab[m][j]<-1e-12){col=j;break;}if(col==permitted)return;
   std::size_t row=m;double ratio=std::numeric_limits<double>::infinity();for(std::size_t i=0;i<m;++i)if(tab[i][col]>1e-13){const double z=tab[i][total]/tab[i][col];if(z<ratio-1e-12||(std::abs(z-ratio)<=1e-12&&(row==m||basis[i]<basis[row]))){ratio=z;row=i;}}
   if(row==m){throw std::runtime_error("unbounded LP");}pivot(row,col);
  }};
 std::vector<double> cost(total,0);for(std::size_t j=n+m;j<total;++j){cost[j]=-1;}objective(cost);optimize(total);
 if(tab[m][total]<-1e-9)throw std::runtime_error("infeasible compatibility region");
 for(std::size_t i=0;i<m;++i)if(basis[i]>=n+m){std::size_t col=n+m;for(std::size_t j=0;j<n+m;++j)if(std::abs(tab[i][j])>1e-12){col=j;break;}if(col==n+m)throw std::runtime_error("redundant artificial row in LP; simplify constraints");pivot(i,col);}
 std::fill(cost.begin(),cost.end(),0);for(std::size_t j=0;j<n;++j){if(!std::isfinite(c[j]))throw std::invalid_argument("nonfinite LP objective");cost[j]=c[j];}objective(cost);optimize(n+m);
 LPResult result;result.solution.assign(n,0);result.iterations=iterations;
 for(std::size_t i=0;i<m;++i)if(basis[i]<n){result.solution[basis[i]]=std::max(0.0,tab[i][total]);}
 for(std::size_t j=0;j<n;++j){result.objective+=c[j]*result.solution[j];}
 for(std::size_t i=0;i<m;++i){double value=-b[i];for(std::size_t j=0;j<n;++j){value+=A[i][j]*result.solution[j];}result.residual=std::max(result.residual,value);}
 LPMatrix bt(m,std::vector<double>(m));std::vector<double> cb(m);
 for(std::size_t i=0;i<m;++i){cb[i]=cost[basis[i]];for(std::size_t j=0;j<m;++j){bt[i][j]=initial[j][basis[i]];}}
 const auto y=lp_solve_system(bt,cb);for(std::size_t i=0;i<m;++i){result.dual_objective+=rhs[i]*y[i];}
 result.gap=std::abs(result.objective-result.dual_objective);
 if(result.residual>1e-8||result.gap>1e-8){throw std::runtime_error("LP primal-dual verification failed");}return result;
}
struct BoundResult {double lower=0,upper=0,gap=0,residual=0;std::vector<double> lower_law,upper_law;};
inline BoundResult compatible_bounds(const std::vector<double>& p,const std::vector<double>& numerator,const std::vector<double>& denominator,double kappa,double rho,const std::vector<double>& recorded_lower,const std::vector<double>& recorded_upper,const std::vector<std::size_t>& observation_map={}){
 const std::size_t j=numerator.size(),m=p.size();
 if(j<1||j>64||m<1||denominator.size()!=j||recorded_lower.size()!=m||recorded_upper.size()!=m||!std::isfinite(kappa)||kappa<1||!std::isfinite(rho)||rho<=0||rho>=1)throw std::invalid_argument("invalid finite-support bound inputs");
 std::vector<std::size_t> map=observation_map;
 if(map.empty()){if(m!=j)throw std::invalid_argument("observation map is required for coarsened support");for(std::size_t k=0;k<j;++k){map.push_back(k);}}
 if(map.size()!=j)throw std::invalid_argument("observation map length differs from complete support");
 const std::size_t n=2*j+1;LPMatrix A;std::vector<double>b;double known_den=0,mass=0;
 auto le=[&](std::vector<double> row,double rhs){A.push_back(row);b.push_back(rhs);};
 auto eq=[&](std::vector<double> row,double rhs){le(row,rhs);for(double& x:row){x=-x;}le(row,-rhs);};
 const double ratio=(1-rho)/rho;
 std::vector<double> min_den(m,std::numeric_limits<double>::infinity());
 for(std::size_t k=0;k<j;++k){
  if(map[k]>=m||!std::isfinite(numerator[k])||!std::isfinite(denominator[k])||denominator[k]<0)throw std::invalid_argument("invalid reward or observation map");
  min_den[map[k]]=std::min(min_den[map[k]],denominator[k]);
  std::vector<double> row(n,0);row[j+k]=1;row[k]=-ratio*kappa;le(row,0);
  row.assign(n,0);row[k]=ratio/kappa;row[j+k]=-1;le(row,0);
 }
 for(std::size_t h=0;h<m;++h){
  if(!std::isfinite(p[h])||p[h]<0||!std::isfinite(recorded_lower[h])||!std::isfinite(recorded_upper[h])||recorded_lower[h]<0||recorded_lower[h]>recorded_upper[h]||recorded_upper[h]>rho+1e-10||!std::isfinite(min_den[h]))throw std::invalid_argument("invalid recorded probabilities, limits, or empty observation fibre");
  mass+=p[h];known_den+=recorded_lower[h]*min_den[h];
  std::vector<double> row(n,0);for(std::size_t k=0;k<j;++k){if(map[k]==h)row[k]=1;}row.back()=-recorded_upper[h];le(row,0);
  row.assign(n,0);for(std::size_t k=0;k<j;++k){if(map[k]==h)row[k]=-1;}row.back()=recorded_lower[h];le(row,0);
 }
 if(std::abs(mass-1)>1e-9)throw std::invalid_argument("recorded law must sum to one");
 if(known_den<=0)throw std::invalid_argument("conditioning denominator lacks a positive recorded lower bound");
 std::vector<double> row(n,0),cost(n,0);for(std::size_t k=0;k<j;++k){row[k]=1;}row.back()=-rho;eq(row,0);
 row.assign(n,0);for(std::size_t k=0;k<j;++k){row[j+k]=1;}row.back()=-(1-rho);eq(row,0);
 row.assign(n,0);for(std::size_t k=0;k<j;++k){row[k]=row[j+k]=denominator[k];cost[k]=cost[j+k]=numerator[k];}eq(row,1);
 const auto upper=simplex_max(A,b,cost);for(double& x:cost){x=-x;}const auto lower=simplex_max(A,b,cost);
 BoundResult result;result.lower=-lower.objective;result.upper=upper.objective;result.gap=std::max(lower.gap,upper.gap);result.residual=std::max(lower.residual,upper.residual);result.lower_law.resize(j);result.upper_law.resize(j);
 if(lower.solution.back()<=0||upper.solution.back()<=0)throw std::runtime_error("invalid Charnes-Cooper scale");
 for(std::size_t k=0;k<j;++k){result.lower_law[k]=(lower.solution[k]+lower.solution[j+k])/lower.solution.back();result.upper_law[k]=(upper.solution[k]+upper.solution[j+k])/upper.solution.back();}
 return result;
}
}
#endif
