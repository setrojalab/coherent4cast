\name{summary.coherent4cast}
\alias{summary.coherent4cast}
\alias{print.coherent4cast}
\alias{print.summary.coherent4cast}
\alias{coef.coherent4cast}
\alias{vcov.coherent4cast}
\alias{logLik.coherent4cast}
\title{Summarize a shared gamma fit}
\description{Display coefficients and regular sandwich intervals.}
\usage{
\method{summary}{coherent4cast}(object, ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Intervals are asymptotic and remain conditional on identifiability and regularity. Nonregular diagnostics are retained.}
\value{A summary.coherent4cast list.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
