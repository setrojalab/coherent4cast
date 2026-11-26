\name{cs_frailty_diagnose}
\alias{cs_frailty_diagnose}
\title{Frailty boundary score and family comparison}
\description{Evaluate a no-frailty boundary score test and optionally refit alternative frailty families.}
\usage{
cs_frailty_diagnose(object, alternatives = FALSE) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{alternatives}{Refit other frailty families and return a descriptive AIC comparison; no nonnested chi-square likelihood-ratio test is applied.}
}
\details{The gamma variance score at zero is ((d-A)^2-d)/2 per complete path. For coarsened data its conditional expectation is computed across compatible paths. Unit-weight fits use the sample nuisance-score projection. Nonunit case weights use the derivative of the weighted variance score and weighted observed information for nuisance adjustment, retaining Louis' correction under coarsening. Weighted calibration requires declared fixed or exogenous weights and a valid weighted-score central limit theorem; arbitrary outcome-dependent selection weights are not validated. The null fit must have full-rank positive information, no coefficient-bound contact, a small score and optimizer convergence; otherwise its p value is withheld. A one-sided positive-part squared score uses the asymptotic half-point-mass, half-chi-square(1) reference only under a regular correctly specified no-frailty null and an identifiable variance direction. This is not a test that proves the frailty family correct. All alternative fits retain the same case weights. Nonnested families receive descriptive AIC comparisons; with arbitrary survey or importance weights the weighted criterion does not inherit the standard AIC risk interpretation. No universal chi-square likelihood-ratio reference is used.}
\value{A list with null_variance_test, null regularity and nuisance-sensitivity diagnostics, comparison table, optional alternative fits and calibration scope.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
