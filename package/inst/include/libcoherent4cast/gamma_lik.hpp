#ifndef LIBcoherent4cast_GAMMA_LIK_HPP
#define LIBcoherent4cast_GAMMA_LIK_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace coherent4cast {
struct GammaResult {
  std::size_t n=0, dimension=0;
  double loglik=0.0;
  std::vector<double> subject_loglik, scores, information, information_individual;
  std::vector<double> counts, cumulative_hazard;
};

inline double finite_double(long double value, const char* message) {
  const double answer=static_cast<double>(value);
  if (!std::isfinite(answer)) throw std::overflow_error(message);
  return answer;
}

// log(1+x)-x/(1+x), evaluated without cancellation near the Poisson limit.
inline long double gamma_curvature(long double x) {
  if (x<0.001L) {
    long double answer=0.0L,power=x*x;
    for (int k=2;k<=12;++k) {
      const long double sign=(k%2==0)?1.0L:-1.0L;
      answer+=sign*static_cast<long double>(k-1)/static_cast<long double>(k)*power;
      power*=x;
    }
    return answer;
  }
  return std::log1p(x)-x/(1.0L+x);
}

// Design is row-major. Unit indices are zero-based, contiguous and nonempty.
// phi_fixed < 0 estimates log(phi), phi_fixed == 0 is the Poisson boundary.
inline GammaResult gamma_evaluate(const std::vector<double>& theta,
  const std::vector<std::size_t>& unit,const std::vector<double>& exposure,
  const std::vector<double>& design,const std::vector<double>& event,
  std::size_t n,std::size_t p,double phi_fixed=-1.0,bool information=true) {
  const std::size_t m=unit.size();
  const bool estimate=phi_fixed<0.0;
  const std::size_t d=p+(estimate?1U:0U);
  if (m==0 || n==0 || p==0 || exposure.size()!=m || event.size()!=m ||
      design.size()!=m*p || theta.size()!=d || !std::isfinite(phi_fixed))
    throw std::invalid_argument("gamma_lik: inconsistent dimensions or dispersion");
  for (double x:theta) if (!std::isfinite(x)) throw std::invalid_argument("gamma_lik: nonfinite parameter");
  const long double phi=estimate?std::exp(static_cast<long double>(theta[p])):static_cast<long double>(phi_fixed);
  if (!(phi>=0.0L) || !std::isfinite(phi) || (estimate && phi==0.0L))
    throw std::invalid_argument("gamma_lik: dispersion outside numerical range");
  std::vector<long double> A(n,0.0L), counts(n,0.0L), evt(n,0.0L), A1(n*p,0.0L), es(n*p,0.0L);
  std::vector<long double> A2(information?n*p*p:0U,0.0L);
  std::vector<bool> seen(n,false);
  for (std::size_t row=0;row<m;++row) {
    const std::size_t i=unit[row];
    if (i>=n || !std::isfinite(exposure[row]) || exposure[row]<0.0 ||
        !std::isfinite(event[row]) || (event[row]!=0.0 && event[row]!=1.0))
      throw std::invalid_argument("gamma_lik: invalid unit, exposure or binary event");
    seen[i]=true;
    long double eta=0.0L;
    for (std::size_t j=0;j<p;++j) {
      const double x=design[row*p+j];
      if (!std::isfinite(x)) throw std::invalid_argument("gamma_lik: nonfinite design");
      eta+=static_cast<long double>(x)*theta[j];
    }
    const long double v=exposure[row]==0.0?0.0L:std::exp(eta+std::log(static_cast<long double>(exposure[row])));
    if (!std::isfinite(v)) throw std::overflow_error("gamma_lik: cumulative hazard overflow");
    A[i]+=v;counts[i]+=event[row];evt[i]+=event[row]*eta;
    for (std::size_t j=0;j<p;++j) {
      const long double x=design[row*p+j];
      A1[i*p+j]+=v*x;es[i*p+j]+=event[row]*x;
      if (information) for (std::size_t k=0;k<p;++k) A2[(i*p+j)*p+k]+=v*x*design[row*p+k];
    }
  }
  if (std::find(seen.begin(),seen.end(),false)!=seen.end()) throw std::invalid_argument("gamma_lik: missing unit index");
  GammaResult out;out.n=n;out.dimension=d;
  out.subject_loglik.resize(n);out.scores.resize(n*d);out.counts.resize(n);out.cumulative_hazard.resize(n);
  if (information) {out.information.assign(d*d,0.0);out.information_individual.assign(n*d*d,0.0);}
  long double ll_total=0.0L;
  for (std::size_t i=0;i<n;++i) {
    const long double ni=counts[i];
    long double ll=evt[i]-A[i], w=1.0L, denom=1.0L;
    long double score_phi=0.0L, info_phi=0.0L;
    if (phi>0.0L) {
      const long double a=1.0L/phi,q=phi*A[i];denom=1.0L+q;
      ll=evt[i]-(a+ni)*std::log1p(q);w=(1.0L+ni*phi)/denom;
      long double sj=0.0L,hj=0.0L;
      const std::size_t jumps=static_cast<std::size_t>(ni);
      for (std::size_t j=0;j<jumps;++j) {
        const long double jp=static_cast<long double>(j)*phi;
        ll+=std::log1p(jp);sj+=jp/(1.0L+jp);hj+=jp/((1.0L+jp)*(1.0L+jp));
      }
      if (estimate) {
        const long double g=gamma_curvature(q);
        score_phi=sj+a*g-ni*q/denom;
        info_phi=-hj+a*g-a*q*q/(denom*denom)+ni*q/(denom*denom);
      }
    }
    out.subject_loglik[i]=finite_double(ll,"gamma_lik: log likelihood overflow");ll_total+=ll;
    out.counts[i]=static_cast<double>(ni);out.cumulative_hazard[i]=finite_double(A[i],"gamma_lik: hazard outside double range");
    for (std::size_t j=0;j<p;++j) out.scores[i*d+j]=finite_double(es[i*p+j]-w*A1[i*p+j],"gamma_lik: score overflow");
    if (estimate) out.scores[i*d+p]=finite_double(score_phi,"gamma_lik: dispersion score overflow");
    if (information) {
      for (std::size_t j=0;j<p;++j) for (std::size_t k=0;k<p;++k) {
        const long double h=w*A2[(i*p+j)*p+k]-(phi>0.0L?w*phi/denom*A1[i*p+j]*A1[i*p+k]:0.0L);
        out.information_individual[(i*d+j)*d+k]=finite_double(h,"gamma_lik: information overflow");
      }
      if (estimate) {
        for (std::size_t j=0;j<p;++j) {
          const double h=finite_double(phi*(ni-A[i])/(denom*denom)*A1[i*p+j],"gamma_lik: cross information overflow");
          out.information_individual[(i*d+j)*d+p]=h;out.information_individual[(i*d+p)*d+j]=h;
        }
        out.information_individual[(i*d+p)*d+p]=finite_double(info_phi,"gamma_lik: dispersion information overflow");
      }
      for (std::size_t j=0;j<d*d;++j) out.information[j]+=out.information_individual[i*d*d+j];
    }
  }
  out.loglik=finite_double(ll_total,"gamma_lik: total log likelihood overflow");
  for(double value:out.information)if(!std::isfinite(value))throw std::overflow_error("gamma_lik: summed information overflow");
  return out;
}
}
#endif
