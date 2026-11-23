\name{rc_band}
\alias{rc_band}
\title{Frequentist Gaussian simulation of nonlinear forecast bands}
\usage{
rc_band(theta, covariance, forecast, draws = 1000L, alpha = 0.05,
  seed = 1L, scale = c("logit", "identity"))
}
\arguments{
\item{theta}{Fitted parameter vector on an unconstrained scale; use log dispersion where appropriate.}
\item{covariance}{Estimated parameter covariance, usually a subject-level sandwich covariance.}
\item{forecast}{Function mapping each parameter vector to a fixed-length numeric forecast vector.}
\item{draws}{Number of Gaussian parameter draws, at least two.}
\item{alpha}{Simultaneous error probability.}
\item{seed}{R random seed; the caller random-number state is restored.}
\item{scale}{Logit for probability forecasts strictly inside (0,1), or identity for other coordinates.}
}
\description{Propagates parameter covariance through a user-supplied nonlinear forecast and forms a simultaneous studentized maximum-statistic band.}
\details{This is a frequentist sandwich-Gaussian approximation, not a Bayesian Laplace posterior or an exact confidence projection. The input estimator and covariance determine validity. Invalid draws are reported rather than silently removed. Logit bands require strictly interior probabilities; structural zeros and ones need a different inferential construction. The function uses the supplied nonlinear forecasting operator on every draw.}
\value{Point forecasts, lower and upper limits, transformed standard deviations, simulated critical value and provenance.}
\examples{
rc_band(0, matrix(.04), function(b) stats::plogis(b + c(-1, 0, 1)),
        draws = 200, seed = 19)
}
