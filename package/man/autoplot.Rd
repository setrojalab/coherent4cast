\name{autoplot}
\alias{autoplot}
\alias{autoplot.coherent4cast_prediction}
\alias{plot.coherent4cast_prediction}
\title{coherent4cast autoplot interface}
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
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
