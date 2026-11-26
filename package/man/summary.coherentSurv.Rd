\name{summary.coherentSurv}
\alias{summary.coherentSurv}
\alias{print.coherentSurv}
\alias{print.summary.coherentSurv}
\alias{coef.coherentSurv}
\alias{vcov.coherentSurv}
\alias{logLik.coherentSurv}
\title{Summarize a shared gamma fit}
\description{Display coefficients and regular sandwich intervals.}
\usage{
\method{summary}{coherentSurv}(object, ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Intervals are asymptotic and remain conditional on identifiability and regularity. Nonregular diagnostics are retained.}
\value{A summary.coherentSurv list.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
