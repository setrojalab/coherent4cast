\name{predict.recurcast_targeted}
\alias{predict.recurcast_targeted}
\title{Forecasts and finite-vector inference from a targeted common law}
\usage{
\method{predict}{recurcast_targeted}(object, newdata = NULL, landmark = 0,
  horizon = NULL, type = "all", interval = "auto", level = .95,
  draws = 1000L, seed = 1L, ...)
}
\arguments{
\item{object}{A finite-path or continuous-time targeted recurcast fit.}
\item{newdata}{Optional continuous histories using the fitted schema.}
\item{landmark}{Nonnegative assessment time. Zero selects baseline population inference when newdata is NULL and the horizon is the fitted target horizon.}
\item{horizon}{Forecast duration; defaults to the prespecified fitted horizon. Continuous forecasts must end no later than that horizon.}
\item{type}{A fitted target name, burden, risk, survival, rmst, or all.}
\item{interval}{auto, none, multiplier, logit_gaussian, or projection. auto selects transformed Gaussian multiplier intervals for the fitted population target.}
\item{level}{Simultaneous confidence level for the fitted finite target vector.}
\item{draws}{Number of independent Gaussian subject multiplier draws, at least 50.}
\item{seed}{Local random seed; the caller's RNG state is restored.}
\item{...}{Reserved; unsupported arguments are rejected.}
}
\description{Returns the fitted population vector with automatic EIF-based uncertainty, or conditional dynamic forecasts under the same targeted law.}
\details{Multiplier intervals use the empirical subject EIF covariance and a max-t critical value. logit_gaussian uses a logit transform for bounded targets and a log transform for nonnegative unbounded burdens before back transformation. These are regular asymptotic approximations. They are not finite-sample confidence projections, do not prove sparse-cell coverage, and are not Bayesian posterior intervals. Failed score solutions, zero observed event marks, zero variances and transformation boundaries retain explicit inference-status flags. The automatic boundary fallback is a separate observed-reward coarsening outer set. For each known-bounded reward, observed lower and upper values receive simultaneous Hoeffding-Bonferroni mean bounds. Death uses observed death and observed death plus early censoring; restricted survival uses observed time and the horizon for censored subjects. Independent subjects and valid prespecified bounds suffice; CAR and estimated censor weights are not required. Uncapped nonnegative burdens retain [0, infinity]. The explicit projection option requests these conservative outer sets even for a regular fit. They are not sharp LP regions or Gaussian confidence intervals.

The data-dependent automatic Gaussian/projection switch does not inherit an unconditional or selection-conditional finite-sample guarantee. Request interval="projection" throughout when invoking the unconditional observed-envelope theorem. The returned inference_scope states whether the construction was selected automatically.

For conditional continuous-history forecasts, the fold-stratum mixture weights are updated by Bayes' rule using the observed history and the targeted intensities. Future rewards are propagated under each retained component, then mixed. The population EIF does not estimate uncertainty for an individual conditional history. Such calls return point forecasts with an explicit conditional-inference status; explicitly requesting population multiplier intervals for them raises an error.}
\value{A recurcast_prediction data frame containing targets, estimates, standard errors and limits when defined, and inference-status labels. Conditional forecasts retain posterior fold weights as an attribute.}
\seealso{\code{targeting_spec}, \code{rc_band}}
