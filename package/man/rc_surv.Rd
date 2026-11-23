\name{rc_surv}
\alias{rc_surv}
\alias{predict.recurcast_surv}
\title{No-frailty right-censored survival interface}
\description{Maps ordinary right-censored survival records to the existing
piecewise-exponential likelihood with frailty variance fixed at zero.}
\usage{
rc_surv(formula, data, time = "time", status = "status",
  unit = NULL, start = NULL, time_unit,
  breaks = c(0, Inf), control = recurcast_control(), ...)
\method{predict}{recurcast_surv}(object, newdata = NULL,
  landmark = 0, horizon, newdata_future = NULL,
  covariate_history = NULL, ...)
}
\arguments{
\item{formula}{A one-sided covariate formula or
\code{Surv(time, status) ~ predictors} / \code{Surv(start, stop, status) ~ predictors}.
The response arguments must be column names; a two-sided formula determines
the response-column mapping. Outcome labels and times cannot be predictors.}
\item{data}{A data frame with exact observed follow-up and covariates.}
\item{time, status}{Names of follow-up and event-status columns for a one-sided
formula. Status is 0 for right censoring and 1 for the absorbing event.}
\item{unit}{Optional independent-unit identifier column. Required for start-stop
records; otherwise rows are independent units.}
\item{start}{Optional interval-start column for start-stop input.}
\item{time_unit}{Explicit unit for every time and horizon, such as \code{"days"}.}
\item{breaks}{Baseline-hazard cut points from zero to infinity.}
\item{control}{Controls from \code{recurcast_control()}. A positive fixed frailty
variance conflicts with this no-frailty interface.}
\item{object}{An object returned by \code{rc_surv()}.}
\item{newdata}{Follow-up records using the original column mapping, or one
covariate-only profile per unit declared event-free at the landmark.}
\item{landmark}{Nonnegative assessment time measured from the time origin.}
\item{horizon}{Positive forecast durations after the landmark.}
\item{newdata_future}{Optional declared future covariate updates; see
\code{predict.recurcast()}. The unit column uses the original mapping, or
\code{unit} when no identifier column was supplied.}
\item{covariate_history}{Optional observed covariate history with the same
unit-column mapping.}
\item{...}{Additional named fit arguments passed to \code{recurcast()}, or
prediction arguments passed to \code{predict.recurcast()}. Schema overrides,
frailty, targeting, coarsening and observation-model overrides are rejected
by the fitting wrapper.}
}
\details{This is an interface to the unchanged zero-variance gamma limit,
not a new likelihood or a Cox partial-likelihood fit. Independent units,
the specified proportional-hazards structure, and a valid censoring
observation model remain assumptions. Observation gaps and delayed entry
are subject to the underlying engine's explicit eligibility checks.
Covariate-only prediction profiles declare eligibility rather than document
observed follow-up. Numerical propagation bounds exclude estimation and
model error. The estimates are prognostic, not treatment effects.}
\value{\code{rc_surv()} returns a \code{recurcast_surv} object inheriting
\code{recurcast}, including coefficients, covariance, convergence information
and audit metadata. Its prediction method returns the ordinary RecurCast
forecast data frame with assumptions and numerical error diagnostics.}
\examples{
d <- data.frame(time = c(1, 2, 3, 4, 2, 4),
                status = c(1, 0, 1, 0, 1, 0))
fit <- rc_surv(~ 1, d, time_unit = "years")
coef(fit)
predict(fit, landmark = 0, horizon = 1, type = "survival")
}
\seealso{\code{\link{recurcast}}, \code{\link{rc_competing}},
\code{\link{rc_audit_report}}}
\keyword{survival}
