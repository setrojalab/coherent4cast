\name{cs_causal}
\alias{cs_causal}
\title{Baseline Treatment Effects on Recurrent and Terminal Event Outcomes}
\usage{
cs_causal(formula, data, treatment, horizon, censor, time_unit,
  unit = "id", start = "start", stop = "stop",
  event_type = "event_type", terminal = "terminal",
  observed = NULL, event_types = NULL, folds = NULL,
  propensity = NULL, spec = list(), seed = 1L, estimand = NULL)
}
\arguments{
\item{formula}{One-sided formula naming baseline confounder strata. Continuous predictors must be discretized into scientifically declared strata or represented through supplied out-of-fold nuisance arrays.}
\item{data}{Long-format event histories, grouped by independent subject and ordered in time.}
\item{treatment}{Name of the baseline binary 0/1 treatment column, fixed within subject.}
\item{horizon}{Positive follow-up duration measured from the common baseline.}
\item{censor}{Name of a 0/1 right-censor indicator. Censoring terminates a subject history and cannot coincide with an exact event.}
\item{time_unit}{Declared time unit, for example \code{"years"}.}
\item{unit,start,stop,event_type,terminal,observed,event_types}{Trajectory columns and recurrent mark labels, as in \code{\link{cs_as_trajectory}}. Delayed entry, observation gaps, coarsened times and simultaneous distinct jumps are not supported by this causal entry point. Zero recurrent marks are allowed for a terminal-only analysis.}
\item{folds}{Optional subject fold labels in first-appearance order, or a column name constant within subject. By default folds are balanced within baseline stratum and treatment arm. Every training fold/stratum needs both arms.}
\item{propensity}{Probability of treatment 1: \code{NULL} for cross-fitted stratum proportions with a 0.5 pseudocount; one known probability for a randomized design; or an externally estimated fold-by-stratum matrix. Supplied estimates must exclude validation outcomes. Matrix rows follow sorted fold labels, columns follow sorted stratum keys.}
\item{spec}{Named list. Options: \code{strata} (formula variables), \code{nfolds} (5), \code{grid} (64), \code{history_cap} (2), \code{working} (\code{"loglinear"} or \code{"pooled"}), \code{nonrecurrent} (labels), \code{complete_followup} (FALSE), \code{min_propensity} (0.01), \code{min_censor_survival} (0.001), \code{maxit} (100), and \code{tolerance} (\code{1e-7/sqrt(n)}). \code{initial_rates} and \code{censor_rates} may be lists named \code{"0"} and \code{"1"}; each rate array is marks by living states by time grid by fold-stratum, each censor vector has one nonnegative constant hazard per fold-stratum. Groups are ordered first by fold, then stratum. Marks are recurrent labels followed by the terminal mark. Supplied arrays must be fitted outside the corresponding validation fold.}
\item{seed}{Seed used for fold construction; the caller's random state is restored.}
\item{estimand}{Optional analyst declaration retained in the audit and summary. A named list may specify population, treatment, endpoint, intercurrent events and population summary. It documents the intended question; it does not change the implemented baseline binary-treatment estimand, validate identification or establish regulatory compliance.}
}
\description{Estimates arm-specific mean recurrent counts, terminal risk and restricted mean time alive, standardized over the same baseline population. Arm-specific intensity targeting uses inverse treatment probabilities and inverse censor survival. Treatment contrasts retain the full joint influence-function covariance.}
\details{
This implements the static binary baseline-treatment extension in Supplement S12 of the accompanying methodological manuscript. It is distinct from prognostic \code{predict.coherentSurv} output.

Let \eqn{v_j^w(0,X)} be the arm-specific backward value for outcome \eqn{j}, \eqn{B_{jk}^w} its jump increment, \eqn{e_w(X)} the treatment probability, and \eqn{G_w(t)} censor survival. The gradient is
\deqn{D_j^w = v_j^w(0,X)-\psi_j(w) + I(W=w)e_w(X)^{-1}\sum_k\int B_{jk}^w/G_w\{dN_k-\lambda_k^w dt\}.}
Here \eqn{W} is baseline treatment, \eqn{X} baseline covariates, \eqn{N_k} the type-\eqn{k} counting process, \eqn{\lambda_k^w} its arm-specific intensity, and integrals stop at the observed censor or terminal time. Exponential rate updates use \eqn{B_{jk}^w/(e_wG_w)}. Every subject contributes baseline standardization to both arms; only the observed arm contributes its event martingale.

Causal identification requires consistency, no interference, conditional treatment exchangeability, positivity, and conditionally independent censoring. The finite retained count state must be sufficient. These assumptions cannot be established by software checks. The default nuisance models are stratum-specific log-linear count-state event rates and constant censor hazards, estimated from training folds. Rates are piecewise constant after targeting. External arrays support alternative out-of-fold event estimates on the same declared finite state.

Double robustness concerns a correct event-intensity model or a correct product of treatment probability and censor survival. It does not protect against unmeasured treatment confounding. Efficiency and Wald/multiplier inference are asymptotic, requiring nuisance convergence, a sufficiently small product remainder, and regular score solutions. \code{interpolation_error} is a sampled numerical diagnostic, not a certified uniform error bound; refine \code{grid} to assess stability. Nonregular coordinates retain point results but do not receive regular confidence intervals.

The implementation supports at most 30 baseline strata, six recurrent marks and 128 living states. Time-varying treatment regimes, baseline treatment-confounding sensitivity bounds, arbitrary censor-history regression and automatic coarsened causal likelihoods are not implemented by this function.
}
\value{An object of class \code{coherentSurv_causal} containing arm laws, standardized estimates, differences, individual arm and contrast EIFs, joint covariance, nuisance arrays, fold provenance, overlap diagnostics, numerical diagnostics, assumptions and convergence flags. \code{predict} and \code{tidy} return contrast tables; \code{coef} and \code{vcov} return differences and their covariance. \code{\link{cs_audit_report}} exports its recorded provenance.}
\references{Rytgaard, H. C., Gerds, T. A. and van der Laan, M. J. (2022). Continuous-time targeted minimum loss-based estimation of intervention-specific mean outcomes. \emph{The Annals of Statistics}, 50, 2469--2491. \doi{10.1214/21-AOS2114}.}
\seealso{\code{\link{cs_causal_contrast}}, \code{\link{cs_audit_report}}}
\examples{
d <- cs_simulate(80, horizon = 1, rates = c(.9, .5, .4),
                 dispersion = 0, censor_rate = 0, seed = 71)
d$treatment <- as.integer(d$id \%\% 2)
d$censored <- 0L
fit <- cs_causal(~1, d, "treatment", 1, "censored", "years",
  propensity = .5, spec = list(complete_followup = TRUE, working = "pooled"))
predict(fit, interval = "none")
fit$diagnostics
}

