\name{tidy}
\alias{tidy}
\alias{tidy.coherent4cast}
\title{coherent4cast tidy interface}
\description{S3 or convenience interface for tidy output.}
\usage{
tidy(x, ...) 
}
\arguments{
\item{x}{A fitted model or prediction object.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Arguments are dispatched to the applicable model or finite-support procedure. A fitted continuous-time likelihood is not silently treated as a finite-support sensitivity experiment.}
\value{A model summary, forecast, sensitivity region or calibration result, according to the method.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
