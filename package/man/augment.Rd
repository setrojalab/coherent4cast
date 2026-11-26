\name{augment}
\alias{augment}
\alias{augment.coherentSurv}
\title{coherentSurv augment interface}
\description{S3 or convenience interface for augment output.}
\usage{
augment(x, ...) 
}
\arguments{
\item{x}{A fitted model or prediction object.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Arguments are dispatched to the applicable model or finite-support procedure. A fitted continuous-time likelihood is not silently treated as a finite-support sensitivity experiment.}
\value{A model summary, forecast, sensitivity region or calibration result, according to the method.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
