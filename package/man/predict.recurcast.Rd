\name{predict.recurcast}
\alias{predict.recurcast}
\title{Landmark forecasts from a common law}
\description{Forecast terminal risk, survival, recurrent burden and restricted mean time alive.}
\usage{
\method{predict}{recurcast}(object, newdata = NULL, landmark, horizon, 
    type = c("risk", "burden", "survival", "rmst", 
        "all"), quadrature = object$control$quadrature, 
    covariate_evolution = c("freeze", "locf", "joint_submodel"), 
    newdata_future = NULL, covariate_model = NULL, 
    covariate_draws = 100L, covariate_grid = 16L, seed = 1L, 
    covariate_history = NULL, newdata_observation = c("exact", 
        "coarsened"), coarsening_control = list(), 
    ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{newdata}{Optional histories following the fitted data schema. Only information available by the landmark is used.}
\item{landmark}{A single nonnegative assessment time; terminal-event-free and observed eligibility is enforced.}
\item{horizon}{One or more positive forecast durations; a single duration for validation or simulation.}
\item{type}{Requested risk, survival probability, restricted mean time alive (rmst), recurrent burden, or all outputs.}
\item{quadrature}{Starting node count for frailty integration: generalized Gauss-Laguerre for gamma forecasts and Gaussian quadrature for lognormal integration. Lognormal fitting can refine this count.}
\item{covariate_evolution}{Freeze landmark covariates, carry forward explicitly supplied scheduled updates (locf), or integrate a Gaussian longitudinal joint_submodel for future covariates.}
\item{newdata_future}{A deterministic scenario table with unit column, time and every formula covariate. Updates must be known or stipulated at forecast origin, not leaked future measurements.}
\item{covariate_model}{An object returned by rc_covariate_model, trained on an appropriate estimation sample.}
\item{covariate_draws}{Number of conditional random-coefficient trajectory draws (at least two); covariate_mcse is reported separately.}
\item{covariate_grid}{Number of piecewise-constant future covariate segments; assess discretization by increasing this grid.}
\item{seed}{Integer seed recorded in provenance or used by a simulation.}
\item{covariate_history}{Explicit measured unit, time and variable table for joint_submodel prediction. Only measurements at or before landmark are used; event interval rows are not repeated biomarker measurements.}
\item{newdata_observation}{Whether the supplied pre-landmark newdata is exact or coarsened. Coarsened histories are expanded automatically and their posterior mixture is retained.}
\item{coarsening_control}{Named arguments passed to rc_coarsen_fibre when legal histories are generated automatically. Declare the actual event-time windows or rounding unit.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{The frailty posterior is conditioned on the observed history and retained across all ordered CSR propagation segments. Future covariates can be frozen, carried forward through prespecified scenario updates, or integrated over Gaussian random-coefficient draws conditional on past measurements. Monte Carlo standard errors and uniformization bounds are distinct from quadrature and covariate-grid errors; increase nodes, draws and grid resolution to check stability. Pre-landmark coarsened newdata is automatically expanded and posterior history mixtures are retained. Full training fibres cannot be reused as landmark histories because they may contain future information.}
\value{A recurcast_prediction data frame with unit, landmark, horizon, target, estimate, time unit, covariate Monte Carlo error and numerical diagnostics. Excluded ineligible units are recorded in an attribute.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
