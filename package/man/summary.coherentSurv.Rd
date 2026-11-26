\name{summary.coherentSurv}
\alias{summary.coherentSurv}
\alias{print.coherentSurv}
\alias{print.summary.coherentSurv}
\alias{coef.coherentSurv}
\alias{vcov.coherentSurv}
\alias{logLik.coherentSurv}
\title{Summarize a fitted coherent event-history law}
\description{Display coefficients and regular sandwich intervals.}
\usage{
\method{summary}{coherentSurv}(object, ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Intervals are asymptotic and remain conditional on identifiability and
regularity. If \code{object$convergence$valid} is not true, regular Wald limits
are withheld as \code{NA} and \code{inference_status} is
\code{"unavailable_nonregular_fit"}. Estimates and numerical standard errors
remain available for diagnosis; those standard errors do not certify valid
inference. The status is retained in the coefficient table and in
\code{tidy()} output. Regular fits receive \code{"regular_asymptotic"}.}
\value{A summary.coherentSurv list.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
