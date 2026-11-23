\name{validate.recurcast_competing}
\alias{validate.recurcast_competing}
\title{Validate a cumulative-incidence forecast}
\description{Evaluates a selected cause's cumulative-incidence forecasts with
inverse-censoring weights while retaining earlier competing events as
observed noncases.}
\usage{
\method{validate}{recurcast_competing}(object, time, status,
  horizon, cause, predicted, censor_survival = NULL,
  threshold = 0.2, strata = NULL, ...)
}
\arguments{
\item{object}{A fitted \code{recurcast_competing} object.}
\item{time}{Aligned nonnegative follow-up times in the fitted time unit.}
\item{status}{Aligned cause or censoring labels from the fitted mapping.}
\item{horizon}{One positive finite assessment horizon on the same time scale.}
\item{cause}{One fitted absorbing cause defining the binary event-by-horizon target.}
\item{predicted}{Aligned predicted cumulative incidence for that cause and horizon.}
\item{censor_survival}{\code{NULL} for a marginal censoring Kaplan-Meier estimate,
or a list with vectors \code{event} containing \eqn{G(T_i-)} and \code{horizon}
containing \eqn{G(h)}, both of length equal to \code{time}. Here \eqn{G} is the
censoring-survival function, \eqn{T_i} is unit \eqn{i}'s observed time, and
\eqn{h} is the assessment horizon. Conditional probabilities must be correctly
specified and aligned.}
\item{threshold}{Probability threshold for classification and net benefit.}
\item{strata}{Optional aligned labels for stratified reporting.}
\item{...}{Unsupported extra arguments are rejected.}
}
\details{A selected-cause failure at or before the horizon is a case. An
earlier failure from another absorbing cause is an observed noncase with
weight based on censoring survival immediately before its failure time.
An event-free observation beyond the horizon is a noncase weighted at the
horizon. Censoring before an outcome can be established does not create a
known noncase. The method delegates metrics to \code{validate.default()}
after this mapping, including Brier score, weighted AUROC, calibration and
threshold-based summaries. Undefined denominators retain the underlying
diagnostics. Validation needs held-out or otherwise appropriately assessed
forecasts; in-sample use does not establish external performance.

The default weights assume marginal independent right censoring. Supplying
conditional weights replaces that assumption with the validity and support
conditions of the supplied censoring model. Competing outcomes themselves
are never treated as independent censoring for the cumulative-incidence target.}
\value{The validation list, or a list by stratum, annotated with the target
cause, endpoint definition and censoring assumptions.}
\seealso{\code{\link{rc_competing}}, \code{\link{validate}}}
\keyword{survival}
