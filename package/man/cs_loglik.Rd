\name{cs_loglik}
\alias{cs_loglik}
\title{Exact Shared Gamma Likelihood and Analytic Derivatives}
\usage{cs_loglik(theta, unit, exposure, design, event, phi_fixed = NULL,
  information = TRUE)}
\arguments{
\item{theta}{Numeric coefficients followed by log frailty variance when estimated.}
\item{unit}{Contiguous positive integer subject indices for expanded rows.}
\item{exposure}{Nonnegative exposure lengths.}
\item{design}{Numeric expanded-row design matrix.}
\item{event}{Binary marked-event indicators.}
\item{phi_fixed}{Fixed nonnegative frailty variance or \code{NULL} to estimate it.}
\item{information}{Whether to compute summed and individual observed information.}}
\description{Evaluates the shared-gamma exact-history likelihood, analytic subject scores and observed information. Independent sampling units are subjects, not expanded exposure rows.}
\details{This implements equation (3) of the accompanying theory manuscript. Integer event-count recurrences evaluate gamma-function differences without subtracting large log-gamma values. A fixed variance of zero is the Poisson boundary. Events require eligibility in the supplied design. Information matrices describe the working finite-dimensional likelihood and are not a general semiparametric efficiency bound.}
\value{A list with \code{loglik}, negative-mean \code{objective} and \code{gradient}, \code{subject_loglik}, \code{scores}, \code{counts}, \code{cumulative_hazard}, and optionally \code{information} and \code{information_individual}.}
\examples{
cs_loglik(c(log(0.3), log(0.5)), 1:3, rep(1, 3),
          matrix(1, 3, 1), c(1, 0, 1))
}
