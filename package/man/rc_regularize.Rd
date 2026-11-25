\name{rc_regularize}
\alias{rc_regularize}
\title{Observed-information log-determinant stabilization}
\description{Refit an exact-history shared-gamma model with one-half the log determinant of observed information added to its log likelihood.}
\usage{
rc_regularize(object, maxit = 120L, tolerance = 1e-08) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{maxit}{Maximum optimizer iterations.}
\item{tolerance}{Optimization convergence tolerance.}
}
\details{Requires estimated positive dispersion and unweighted exact histories. Optimization is restricted to the positive-information domain. The log-determinant score is analytic; penalized curvature uses central differences. Sandwich variance uses centered independent-subject scores. This is local regularization, not automatic first-order Firth bias cancellation, a Bayesian posterior or a finite-sample coverage guarantee. The returned fit can be used with predict and rc_band.}
\value{An updated coherent4cast object with penalized convergence diagnostics, sandwich covariance and regularization provenance.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
