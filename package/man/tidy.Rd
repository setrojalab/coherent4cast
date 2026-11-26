\name{tidy}
\alias{tidy}
\alias{tidy.coherentSurv}
\title{coherentSurv tidy interface}
\description{S3 or convenience interface for tidy output.}
\usage{
tidy(x, ...) 
}
\arguments{
\item{x}{A fitted model or prediction object.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Arguments are dispatched to the applicable model or finite-support procedure. A fitted continuous-time likelihood is not silently treated as a finite-support sensitivity experiment. For ordinary coherentSurv fits, the coefficient table includes \code{inference_status}; regular Wald limits are withheld as \code{NA} when the fit has not met convergence and regularity criteria. See \code{summary.coherentSurv}.}
\value{A model summary, forecast, sensitivity region or calibration result, according to the method.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
