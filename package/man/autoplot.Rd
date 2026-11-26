\name{autoplot}
\alias{autoplot}
\alias{autoplot.coherentSurv_prediction}
\alias{plot.coherentSurv_prediction}
\title{coherentSurv autoplot interface}
\description{S3 or convenience interface for autoplot output.}
\usage{
autoplot(object, ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Arguments are dispatched to the applicable model or finite-support procedure. A fitted continuous-time likelihood is not silently treated as a finite-support sensitivity experiment.}
\value{A model summary, forecast, sensitivity region or calibration result, according to the method.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
