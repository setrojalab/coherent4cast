\name{rc_competing}
\alias{rc_competing}
\alias{predict.coherent4cast_competing}
\alias{print.coherent4cast_competing}
\alias{coef.coherent4cast_competing}
\alias{vcov.coherent4cast_competing}
\alias{logLik.coherent4cast_competing}
\alias{tidy.coherent4cast_competing}
\alias{summary.coherent4cast_competing}
\alias{glance.coherent4cast_competing}
\title{Cause-specific hazards and coherent competing-risk forecasts}
\description{Fits separate no-frailty cause-specific piecewise-exponential
hazards and integrates all fitted causes jointly to obtain cumulative
incidence and event-free survival.}
\usage{
rc_competing(formula, data, time = "time", status = "status",
  unit = NULL, start = NULL, causes = NULL, censor_code = 0,
  time_unit, breaks = c(0, Inf),
  control = coherent4cast_control(), ...)
\method{predict}{coherent4cast_competing}(object, newdata = NULL,
  landmark = 0, horizon,
  type = c("cif", "survival", "risk", "rmst", "all"),
  covariate_evolution = c("freeze", "locf"),
  newdata_future = NULL, ...)
\method{print}{coherent4cast_competing}(x, ...)
\method{coef}{coherent4cast_competing}(object, ...)
\method{vcov}{coherent4cast_competing}(object, ...)
\method{logLik}{coherent4cast_competing}(object, ...)
\method{tidy}{coherent4cast_competing}(x, ...)
\method{summary}{coherent4cast_competing}(object, ...)
\method{glance}{coherent4cast_competing}(x, ...)
}
\arguments{
\item{formula}{A one-sided covariate formula or a supported \code{Surv}
formula with column names as response arguments. See \code{rc_surv()}.}
\item{data}{Exact observed first-event histories and predictors. Causes are
mutually exclusive absorbing outcomes.}
\item{time, status}{Follow-up and cause-label column names for a one-sided formula.}
\item{unit}{Optional independent-unit column; required for start-stop records.}
\item{start}{Optional interval-start column. Histories must begin at time zero
and be contiguous; this interface rejects delayed entry and observation gaps.}
\item{causes}{Distinct nonzero absorbing-cause labels, excluding the censoring
code. When \code{NULL}, labels are inferred from observed noncensoring statuses.
Each declared cause must have an observed event. Label \code{"0"} is reserved.}
\item{censor_code}{One explicit censoring label, also required to match a
factor's censoring level when the status column is a factor.}
\item{time_unit}{Common explicit unit for all times and horizons.}
\item{breaks}{Baseline-hazard cut points from zero to infinity.}
\item{control}{Controls from \code{coherent4cast_control()}. Frailty variance is
fixed at zero by this interface.}
\item{object, x}{A fitted \code{coherent4cast_competing} object.}
\item{newdata}{Records using the fitted mapping or eligible covariate-only
profiles. A profile declares event-free eligibility at the landmark.}
\item{landmark}{Nonnegative assessment time.}
\item{horizon}{Positive finite durations after the landmark.}
\item{type}{Cause-specific cumulative incidence, event-free survival,
any-cause risk, restricted mean event-free time, or all targets.}
\item{covariate_evolution}{Freeze the assessment-time covariates, or apply
declared last-observation-carried-forward updates.}
\item{newdata_future}{Declared future covariate updates with the original
unit-column mapping, or \code{unit} when rows supplied the identifiers.}
\item{...}{Named fit arguments forwarded to \code{coherent4cast()}; fixed schema,
observation and no-frailty settings cannot be overridden. The prediction
method rejects unsupported arguments. Other methods accept generic arguments.}
}
\details{Each cause-specific fit treats other absorbing outcomes as exits
from its risk set. Prediction then integrates all cause-specific rates in one
absorbing-state generator. Consequently, cumulative incidences sum to
any-cause risk, and event-free survival plus any-cause risk sums to one up to
the reported numerical error. This is not one minus separate cause-specific
survival curves.

The joint sandwich covariance stacks subject scores across causes and uses
block-diagonal observed information; cross-cause covariance is retained.
Coefficient intervals are regular asymptotic approximations. No automatic
cumulative-incidence confidence bands are returned. Independent units,
valid censoring assumptions and the specified hazard model remain required.
The log likelihood is the sum of fitted cause-specific contributions. With
arbitrary case weights, its AIC-style summary is a weighted descriptive
criterion and need not have ordinary likelihood calibration. The forecasts
are prognostic; competing-risk mapping does not identify intervention effects.}
\value{A \code{coherent4cast_competing} object containing cause-specific fits,
joint coefficients and covariance, schema, convergence and audit information.
Prediction returns a long forecast data frame including targets,
assumptions and numerical propagation error bounds.}
\examples{
d <- data.frame(time = c(1, 2, 3, 4, 2, 4),
                status = c("a", "b", "0", "a", "b", "0"))
fit <- rc_competing(~ 1, d, causes = c("a", "b"),
                    time_unit = "years")
coef(fit)
predict(fit, landmark = 0, horizon = 1, type = "all")
}
\seealso{\code{\link{rc_surv}}, \code{\link{validate.coherent4cast_competing}},
\code{\link{rc_audit_report}}}
\keyword{survival}
