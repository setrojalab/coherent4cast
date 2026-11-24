#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
// General mean-one discrete quadrature likelihood. Louis' identity is applied
// within each subject; derivatives of nodes and weights are supplied analytically.
// [[Rcpp::export]]
List cpp_rc_mixture_lik(NumericVector beta, IntegerVector unit, NumericVector exposure,
 NumericMatrix X, NumericVector event, NumericVector h, NumericVector g,
 NumericMatrix dh, NumericMatrix dg, NumericVector d2h, NumericVector d2g, bool information) {
 int m=X.nrow(),p=X.ncol(),q=dh.ncol(),L=h.size(),n=max(unit),P=p+q;
 if(beta.size()!=p || unit.size()!=m || exposure.size()!=m || event.size()!=m ||
    g.size()!=L || dh.nrow()!=L || dg.nrow()!=L || dg.ncol()!=q ||
    d2h.size()!=L*q*q || d2g.size()!=L*q*q) stop("Incompatible mixture dimensions");
 NumericVector A(n),D(n),E(n),A1(n*p),A2(information?n*p*p:0),EX(n*p);
 for(int r=0;r<m;++r){int i=unit[r]-1; if(i<0||i>=n||!std::isfinite(exposure[r])||exposure[r]<0||(event[r]!=0&&event[r]!=1))stop("Invalid exposure/event/unit");
   double eta=0; for(int a=0;a<p;++a){if(!std::isfinite(X(r,a)))stop("Nonfinite design");eta+=X(r,a)*beta[a];}
   double z=exposure[r]*std::exp(eta);if(!std::isfinite(z))stop("Likelihood overflow");
   A[i]+=z;D[i]+=event[r];E[i]+=event[r]*eta;
   for(int a=0;a<p;++a){A1[i+n*a]+=z*X(r,a);EX[i+n*a]+=event[r]*X(r,a);
     if(information)for(int b=0;b<p;++b)A2[i+n*(a+p*b)]+=z*X(r,a)*X(r,b);}
 }
 NumericVector ll(n),ind(information?P*P*n:0),suminfo(information?P*P:0);
 NumericMatrix scores(n,P),post(n,L);
 for(int i=0;i<n;++i){
   std::vector<double> v(L),s(L*P),H(information?L*P*P:0); double mx=R_NegInf;
   for(int l=0;l<L;++l){double z=std::exp(h[l]),za=z*A[i];v[l]=g[l]+D[i]*h[l]-za;mx=std::max(mx,v[l]);
     for(int a=0;a<p;++a){s[l+L*a]=EX[i+n*a]-z*A1[i+n*a];
       if(information){for(int b=0;b<p;++b)H[l+L*(a+P*b)]=-z*A2[i+n*(a+p*b)];
         for(int b=0;b<q;++b){double x=-z*A1[i+n*a]*dh(l,b);H[l+L*(a+P*(p+b))]=x;H[l+L*(p+b+P*a)]=x;}}}
     for(int a=0;a<q;++a){s[l+L*(p+a)]=dg(l,a)+(D[i]-za)*dh(l,a);
       if(information)for(int b=0;b<q;++b)H[l+L*(p+a+P*(p+b))]=d2g[l+L*(a+q*b)]+(D[i]-za)*d2h[l+L*(a+q*b)]-za*dh(l,a)*dh(l,b);}
   }
   if(!std::isfinite(mx)) {
     stop("Nonfinite quadrature likelihood");
   }
   double den=0;
   for(int l=0;l<L;++l)den+=std::exp(v[l]-mx);
   ll[i]=E[i]+mx+std::log(den);
   for(int l=0;l<L;++l){post(i,l)=std::exp(v[l]-mx)/den;for(int a=0;a<P;++a)scores(i,a)+=post(i,l)*s[l+L*a];}
   if(information)for(int a=0;a<P;++a)for(int b=0;b<P;++b){double info=scores(i,a)*scores(i,b);
     for(int l=0;l<L;++l)info-=post(i,l)*(H[l+L*(a+P*b)]+s[l+L*a]*s[l+L*b]);
     ind[a+P*(b+P*i)]=info;suminfo[a+P*b]+=info;}
 }
 NumericVector gradient(P);for(int a=0;a<P;++a)for(int i=0;i<n;++i)gradient[a]-=scores(i,a)/n;
 if(information){ind.attr("dim")=IntegerVector::create(P,P,n);suminfo.attr("dim")=IntegerVector::create(P,P);}
 double total_ll=Rcpp::sum(ll);
 return List::create(_["subject_loglik"]=ll,_["scores"]=scores,_["information_individual"]=ind,
 _["information"]=suminfo,_["posterior_weights"]=post,_["counts"]=D,_["cumulative_hazard"]=A,
 _["loglik"]=total_ll,_["objective"]=-total_ll/n,_["gradient"]=gradient);
}
