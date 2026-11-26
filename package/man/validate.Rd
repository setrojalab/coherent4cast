\name{validate}
\alias{validate}
\alias{validate.default}
\alias{validate.coherentSurv}
\title{Validate dynamic predictions}
\description{Calculate censoring-adjusted discrimination, calibration and threshold statistics.}
\usage{
validate(object, ...) 
\method{validate}{default}(object, time, status, horizon,
  censor_survival = NULL, threshold = 0.2, strata = NULL, ...)
\method{validate}{coherentSurv}(object, time, status, horizon, predicted,
  censor_survival = NULL, ...)
}
\arguments{
\item{object}{A numeric vector of terminal-event risks for the default method,
or a fitted model for the model method. A long prediction data frame must first
be restricted to one risk target and horizon and aligned to the evaluation units.}
\item{time}{Nonnegative follow-up durations measured from the assessment
landmark, one per eligible independent evaluation unit. Subtract the landmark
from absolute endpoint times before calling this function.}
\item{status}{Binary indicator of an observed terminal event: one for an event,
zero for right censoring.}
\item{horizon}{One positive forecast duration after the same assessment
landmark, in the same units as \code{time}.}
\item{predicted}{Aligned terminal-event risks for the fitted-model method.
These must be supplied explicitly; the method does not generate predictions.}
\item{censor_survival}{Optional list with vectors \code{case = G(time-)} and
\code{horizon = G(horizon)}, each aligned to the evaluation units. Needed
probabilities must be positive. The default is a marginal Kaplan--Meier
censoring estimate, requiring marginal independent censoring. Supplied
conditional probabilities require an appropriate censoring model.}
\item{threshold}{One decision threshold strictly between zero and one.}
\item{strata}{Optional nonmissing group labels for separate validation
calculations. Without supplied censoring probabilities, the censoring estimate
is recalculated within each group.}
\item{...}{For the model method, \code{threshold} and \code{strata} are
forwarded to the default method. Other arguments are rejected; in particular,
\code{landmark} is not an accepted argument.}
}
\details{Validation is at a fixed horizon among independently sampled eligible
units. The caller must select units that are observed and terminal-event-free
at the landmark, retain one row per unit, align the risk estimates and outcomes,
and rebase follow-up to that landmark. Predictions should come from held-out
subjects. This routine does not claim optimism correction for training-set
predictions or a causal net benefit. Undefined metrics are labelled rather than
silently replaced. A calibration slope with complete or quasi-complete
separation is returned as \code{NA} with \code{calibration_status = "separated"};
constant predictions give \code{"constant_predictions"}. Calibration-in-the-large
and other metrics may remain defined. Nonconverged calibration coefficients
are also withheld. Numerical convergence alone does not make a separated
logistic slope finite.}
\value{A list with metric estimates, a weighted two-by-two table, eligibility counts, censoring weights and assumptions.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
\examples{
validate(c(.1,.2,.7,.8),c(1,3,2,4),c(1,1,1,0),horizon=2.5)
# Evaluate two-year risks for units eligible at year five.
validate(c(.8,.2),time=c(6,10)-5,status=c(1,0),horizon=2)
}
