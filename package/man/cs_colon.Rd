\name{cs_colon}
\alias{cs_colon}
\title{Day-coarsened illness-death paths for the public colon trial}
\usage{cs_colon(data = NULL, nodes = 4L, covariate = c("node4", "nodes"),
  node_median = NULL, breaks = c(0, 1, 3, Inf))}
\arguments{
\item{data}{Optional source data frame. NULL loads \code{survival::colon}.}
\item{nodes}{Gauss-Legendre nodes per continuous integration segment, an integer from 1 to 64.}
\item{covariate}{Retain supplied binary \code{node4}, or use continuous \code{nodes}/5 with a missingness indicator.}
\item{node_median}{Training-set median for missing continuous nodes. NULL computes the median in the supplied source data. Reuse the training median for held-out data.}
\item{breaks}{Increasing baseline breaks in years, beginning at zero and ending at infinity. The adapter splits quadrature at every discontinuity. Use the same breaks in fitting.}
}
\description{Returns legal quadrature paths for all patients, including same-day recurrence and death, without assigning invented within-day event times or dropping tied patients.}
\value{A list with \code{data}, \code{coarsening=list(owner,log_weights)}, \code{formula}, \code{event_types="1"}, \code{nonrecurrent="1"}, \code{breaks}, and detailed \code{metadata}. Each data identifier is a latent path; \code{source_id} retains the observed patient identifier. \code{owner} maps paths to observed patients in sorted source-ID order.}
\details{Events on reported day d are integrated over (d-1,d]. Censoring endpoints remain fixed. Same-day recurrence/death integrates over the triangle recurrence-time less than death-time. Recurrence after death is not an admissible alternative. Gaussian quadrature approximates the coarsened likelihood; its weights are integration masses, not probabilities normalized within subjects. With four nodes and default breaks, the canonical 929-patient dataset gives 7,447 paths. There are 468 recorded recurrences, 452 recorded deaths, and five same-day pairs. A patient whose recurrence surveillance ends before mortality follow-up is truncated at the earlier endpoint; this assumption is reported in metadata. Use an illness-death submodel with a nonrepeatable recurrence mark. These full-follow-up paths must not be treated as information available at an earlier landmark. Reconstruct a history-only coarsening fibre for prediction.}
\examples{
if (requireNamespace("survival", quietly=TRUE)) {
  x <- cs_colon(nodes=2L)
  x$metadata[c("patients", "same_day_ids")]
  # fit <- cs_fit(x$formula, x$data, time_unit="years",
  #   event_types="1", nonrecurrent="1", breaks=x$breaks,
  #   observation="coarsened", coarsening=x$coarsening)
}
}
\seealso{\code{\link{cs_predict_paths}}, \code{\link{cs_coarsen}}}
\references{The official survival documentation describes the paired endpoints and the unresolved nodes/node4 discrepancy: \url{https://stat.ethz.ch/R-manual/R-devel/library/survival/html/colon.html}.}
