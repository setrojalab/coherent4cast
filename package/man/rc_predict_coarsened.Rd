\name{rc_predict_coarsened}
\alias{rc_predict_coarsened}
\title{Automatic forecasting from coarsened histories at a landmark}
\usage{rc_predict_coarsened(object, history, landmark, horizon,
  type=c("risk","burden","survival","rmst","all"),
  quadrature=object$control$quadrature,coarsening_control=list(),...)}
\arguments{
\item{object}{A fitted coherent4cast model, whether trained on exact or coarsened histories.}
\item{history}{Only records available through the landmark, in the counting-row or event-window format accepted by rc_coarsen_fibre. Every subject must be terminal-event-free and observed through the landmark.}
\item{landmark}{Positive finite assessment time, in the fitted time unit.}
\item{horizon}{Positive forecast horizons in that time unit.}
\item{type}{Requested prognostic summaries.}
\item{quadrature}{Frailty quadrature nodes for forecasting.}
\item{coarsening_control}{Named lower/upper/followup column choices, format, rounding_unit, nodes, max_perms, draws, seed and max_paths. It cannot override the fitted marks, time units or baseline breaks.}
\item{...}{Prediction controls forwarded to rc_predict_paths: covariate_evolution, newdata_future, covariate_model, covariate_history, covariate_draws, covariate_grid and seed. Future schedules and measured histories use the fitted original unit column and a numeric column named time. They are replicated to each candidate path. The joint submodel requires an actual measured covariate_history; event rows are never treated as repeated measurements.}
}
\description{Generates a history-only observation fibre automatically, updates its path mixture and returns dynamic forecasts with original subject identifiers.}
\value{A coherent4cast_prediction data frame. Its unit column contains the original observed identifiers. The automatic_fibre attribute reports per-subject enumeration/sampling diagnostics and the path_weights attribute includes source-unit mappings.}
\details{Every window and follow-up endpoint must lie at or before the landmark; follow-up must end exactly there. A window straddling the landmark is rejected, not silently truncated using future knowledge. Terminal events in the supplied history are rejected because these forecasts condition on survival to the landmark. Candidate complete histories retain their path-specific posterior frailty distributions. Numeric propagation bounds exclude fibre and frailty quadrature error. An importance-sampled fibre remains an approximation. Event counts and marks must be complete under the declared observation model. Measured covariate history is restricted to times through the landmark and requires one actual measurement per unit/time. Covariate models do not turn latent quadrature event times into observed measurement times.}
\seealso{\code{\link{rc_coarsen_fibre}}, \code{\link{rc_predict_paths}}, \code{\link{predict.coherent4cast}}}
