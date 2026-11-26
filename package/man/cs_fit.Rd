\name{cs_fit}
\alias{cs_fit}
\title{Fit a coherent recurrent and terminal event law}
\description{Estimate a common law with shared frailty or an explicitly specified targeted estimation experiment.}
\usage{
cs_fit(formula, data, unit = "id", start = "start", 
    stop = "stop", event_type = "event_type", terminal = "terminal", 
    observed = NULL, event_types = NULL, time_unit, 
    breaks = c(0, Inf), observation = c("exact", "coarsened"), 
    targeting = "none", control = cs_control(), 
    weights = NULL, subset = NULL, na.action = stats::na.fail, 
    labels = NULL, nonrecurrent = character(), coarsening = NULL, 
    intervention = NULL, estimand = NULL, seed = 1L, 
    frailty = c("gamma", "lognormal", "discrete", "none"), 
    frailty_points = 2L, frailty_diagnostics = TRUE, 
    targeting_spec = NULL, coarsening_control = list()) 


}
\arguments{
\item{formula}{A one-sided formula for predictable measured covariates. Interval outcomes are prohibited as covariates.}
\item{data}{A data frame of complete counting-process intervals; one row is one observed interval.}
\item{unit}{Name of the unit identifier column, or a selected unit identifier for plotting.}
\item{start}{Name of the interval start column.}
\item{stop}{Name of the interval end column.}
\item{event_type}{Name of the nonterminal event-label column. Zero or 'none' denotes no event.}
\item{terminal}{Name of the binary absorbing-event column.}
\item{observed}{Optional name of a binary monitoring indicator column. Unobserved intervals cannot be fitted as observed zeros.}
\item{event_types}{Declared nonterminal event labels, including types with no observed events.}
\item{time_unit}{Explicit nonempty unit of follow-up time, such as years or operating hours.}
\item{breaks}{Increasing piecewise-baseline boundaries beginning at zero and ending at Inf.}
\item{observation}{Exact complete histories or an explicitly enumerated coarsened-path experiment.}
\item{targeting}{Likelihood fitting ('none'), finite-path MAR targeting ('finite_path'), or continuous-time CAR intensity targeting ('continuous'). See targeting_spec for the distinct observation models and required specification.}
\item{control}{Settings returned by cs_control.}
\item{weights}{Optional positive case weights, one per independent unit or observed coarsened subject.}
\item{subset}{Optional expression selecting whole units. Partial-history selection is rejected.}
\item{na.action}{na.fail by default. na.omit removes whole affected units and records their identifiers.}
\item{labels}{Optional user-defined target labels stored with the fit.}
\item{nonrecurrent}{Event labels that may occur at most once; their post-event eligibility is zero.}
\item{coarsening}{For enumerated paths, a list with owner (path-to-observed-subject indices) and log_weights (quadrature integration weights).}
\item{intervention}{Reserved extension argument. A non-NULL value raises an explicit unsupported-causal-analysis error.}
\item{estimand}{Optional analyst-supplied declaration retained in the fit audit, summary and audit export. A named list may describe population, treatment, endpoint, intercurrent events and population summary. It does not select a different model, establish identification or certify ICH E9(R1) compliance.}
\item{seed}{Integer seed recorded in provenance or used by a simulation.}
\item{frailty}{Mean-one gamma, lognormal or estimated discrete frailty, or 'none' to fix the shared variance to zero through the existing gamma-kernel limit. A positive phi_fixed conflicts with 'none'.}
\item{frailty_points}{Number of estimated mass points for discrete frailty, from two to eight; three deterministic optimization starts are compared.}
\item{frailty_diagnostics}{Compute a no-frailty boundary score diagnostic and the fitted-model AIC, stored for summary().}
\item{targeting_spec}{A list defining the target experiment, folds, observation process and horizon; see targeting_spec.}
\item{coarsening_control}{Named arguments passed to cs_coarsen_fibre when legal histories are generated automatically. Declare the actual event-time windows or rounding unit.}
}
\details{Likelihood baselines are piecewise constant, with type-specific covariate and capped-history coefficients. Gamma integration is analytic; lognormal integration uses quadrature; discrete frailty estimates mean-one support and weights with three starts. Scores and Louis information are analytic for the numerical likelihood. Coarsened histories are generated automatically when coarsening is NULL; quadrature and importance-sampling error are reported separately. Targeting selects a separate finite-path MAR or continuous-time CAR sufficient-state model, as documented in targeting_spec. Fits at numerical bounds or with deficient information retain visible nonregular diagnostics. Prognostic targets are not causal treatment effects.}
\value{A coherentSurv or coherentSurv_targeted S3 object containing the common law, estimates, covariance, convergence diagnostics, scientific assumptions and provenance.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
\examples{
data(sim_perinatal)
fit <- cs_fit(~x,sim_perinatal,time_unit='months',event_types=c('1','2'))
print(fit)
glance(fit)
}

