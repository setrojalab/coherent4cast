\name{summary.recurcast}
\alias{summary.recurcast}
\alias{print.recurcast}
\alias{print.summary.recurcast}
\alias{coef.recurcast}
\alias{vcov.recurcast}
\alias{logLik.recurcast}
\title{Summarize a shared gamma fit}
\description{Display coefficients and regular sandwich intervals.}
\usage{
\method{summary}{recurcast}(object, ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Intervals are asymptotic and remain conditional on identifiability and regularity. Nonregular diagnostics are retained.}
\value{A summary.recurcast list.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
