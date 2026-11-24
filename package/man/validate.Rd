\name{validate}
\alias{validate}
\alias{validate.default}
\alias{validate.recurcast}
\title{Validate dynamic predictions}
\description{Calculate censoring-adjusted discrimination, calibration and threshold statistics.}
\usage{
validate(object, ...) 
}
\arguments{
\item{object}{A fitted model, prediction object, or numeric risk vector, as appropriate.}
\item{...}{Additional arguments for the relevant S3 method or plot.}
}
\details{Validation is at a fixed horizon among independently sampled eligible units. Predictions should come from held-out subjects. This routine does not claim optimism correction for training-set predictions or a causal net benefit. Undefined metrics are labelled rather than silently replaced.}
\value{A list with metric estimates, a weighted two-by-two table, eligibility counts, censoring weights and assumptions.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
\examples{
validate(c(.1,.2,.7,.8),c(1,3,2,4),c(1,1,1,0),horizon=2.5)
}
