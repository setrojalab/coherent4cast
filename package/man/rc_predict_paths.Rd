\name{rc_predict_paths}
\alias{rc_predict_paths}
\title{Forecast from an explicit mixture of pre-landmark histories}
\usage{rc_predict_paths(object, paths, owner, log_weights, landmark, horizon,
  type = c("risk", "burden", "survival", "rmst", "all"),
  quadrature = object$control$quadrature, ...)}
\arguments{
\item{object}{A fitted \code{recurcast} model.}
\item{paths}{Candidate complete histories, with the same columns and units as the fit. Every path must span exactly zero through the landmark and contain no terminal event. No future records are allowed.}
\item{owner}{Contiguous positive integers mapping each path to an observed subject, ordered by first appearance of path identifiers.}
\item{log_weights}{Finite log integration masses for the pre-landmark observation fibre. They must not encode future follow-up.}
\item{landmark}{Positive finite assessment time in the fitted time unit.}
\item{horizon}{Positive finite prediction horizons.}
\item{type}{Requested prognostic summaries.}
\item{quadrature}{Number of frailty quadrature nodes used for the applicable fitted family.}
\item{...}{Additional arguments passed to \code{predict.recurcast}, including declared future-covariate settings.}
}
\description{Updates the probability of each candidate past history and retains its own posterior frailty distribution when forecasting.}
\value{A \code{recurcast_prediction} data frame with one forecast per observed subject, horizon and target. The \code{path_weights} attribute reports posterior mixture weights. Numerical series and probability roundoff bounds exclude frailty quadrature and observation-fibre quadrature error.}
\details{Path weights are proportional to the path likelihood times its supplied integration mass. Each path is forecast separately using its posterior frailty distribution; forecasts are then averaged. Averaging distribution parameters first is generally incorrect. This helper validates temporal support, but scientific construction of the observation fibre remains the analyst's responsibility. It does not infer hidden paths from arbitrary missing-data patterns and does not use full-follow-up \code{rc_colon()} paths at an earlier landmark.}
\seealso{\code{\link{rc_filter_paths}}, \code{\link{rc_colon}}, \code{\link{predict.recurcast}}}
