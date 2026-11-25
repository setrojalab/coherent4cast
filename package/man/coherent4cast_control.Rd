\name{coherent4cast_control}
\alias{coherent4cast_control}
\alias{recurcast_control}
\title{Estimation and propagation controls}
\description{Declare numerical and sufficient-state settings.}
\usage{
coherent4cast_control(maxit = 300L, tolerance = 1e-08, history_cap = 2L, 
    phi_fixed = NULL, start_phi = 0.5, quadrature = 30L, 
    propagation_tolerance = 1e-12, max_states = 50000L, 
    coefficient_bounds = c(-12, 8), log_phi_bounds = c(-5, 
        2.5), quadrature_tolerance = 1e-07, quadrature_max = 1024L) 


recurcast_control(maxit = 300L, tolerance = 1e-08, history_cap = 2L, 
    phi_fixed = NULL, start_phi = 0.5, quadrature = 30L, 
    propagation_tolerance = 1e-12, max_states = 50000L, 
    coefficient_bounds = c(-12, 8), log_phi_bounds = c(-5, 
        2.5), quadrature_tolerance = 1e-07, quadrature_max = 1024L) 
}
\arguments{
\item{maxit}{Maximum optimizer iterations.}
\item{tolerance}{Optimization convergence tolerance.}
\item{history_cap}{Count-feature cap used for estimation and finite-state forecasting.}
\item{phi_fixed}{A fixed nonnegative gamma or lognormal frailty variance. NULL estimates log dispersion; zero selects the no-frailty boundary. Estimated discrete support requires NULL.}
\item{start_phi}{Positive optimizer starting dispersion, not a clinical assumption.}
\item{quadrature}{Starting node count for frailty integration: generalized Gauss-Laguerre for gamma forecasts and Gaussian quadrature for lognormal integration. Lognormal fitting can refine this count.}
\item{propagation_tolerance}{Poisson-series truncation tolerance for each propagation.}
\item{max_states}{Maximum number of retained finite states, including absorption.}
\item{coefficient_bounds}{Lower and upper numerical bounds for log-rate coefficients.}
\item{log_phi_bounds}{Lower and upper numerical bounds for log dispersion.}
\item{quadrature_tolerance}{Lognormal fitting doubles Gaussian quadrature nodes until maximum subject log-likelihood and mean-score changes satisfy this tolerance. This is a numerical stability diagnostic, not a rigorous integration error bound.}
\item{quadrature_max}{Maximum lognormal quadrature nodes, up to 1024. Unresolved refinement is reported as nonconverged rather than silently accepted.}
}
\details{Defaults are numerical starting values and tolerances. State truncation and history capping change the working intensity family and require scientific justification.}
\value{A named list of control settings.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
\examples{
coherent4cast_control(phi_fixed=0,history_cap=1)
}

\note{The legacy function name \code{recurcast_control} remains an alias for migration from the former package name. New code should use \code{coherent4cast_control}.}
