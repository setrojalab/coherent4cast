\name{cs_colon_fibre}
\alias{cs_colon_fibre}
\title{Automatic coarsening fibres for the public colon trial}
\usage{cs_colon_fibre(data=NULL, nodes=4L, covariate=c("node4","nodes"),
  node_median=NULL, breaks=c(0,1,3,Inf), max_perms=120L,draws=2000L,
  seed=20261122L,max_paths=100000L)}
\arguments{
\item{data}{Optional canonical source data frame; NULL loads survival::colon.}
\item{nodes}{Gauss-Legendre nodes per integration segment.}
\item{covariate}{Preserve supplied node4 or use continuous nodes with a missingness indicator.}
\item{node_median}{Training-set median for missing continuous nodes, or NULL to compute it on supplied data.}
\item{breaks}{Baseline discontinuities in years.}
\item{max_perms,draws,seed,max_paths}{Automatic fibre controls; see \code{cs_coarsen_fibre}.}
}
\description{Transforms the source paired endpoints into event-time windows and calls the general automatic C++ fibre engine. With the default source data, all 929 patients and five same-day pairs are retained.}
\value{A \code{coherentSurv_fibre} object, suggested baseline formula, canonical dataset metadata, and \code{observed_data} containing the source event-window records. The latter can be supplied directly to \code{cs_fit(..., observation="coarsened")} without constructing paths manually.}
\details{A recorded day d defines the window (d-1,d] in days, converted using 365.25 days per year. Recurrence is nonrepeatable, terminal death is last, and a recorded recurrence is never discarded to create a death-only alternative. Early termination of recurrence surveillance truncates follow-up consistently with \code{cs_colon}. Binary node4 is retained as supplied; its known inconsistency with nodes is addressed by an explicit sensitivity choice. Full-follow-up fibres cannot be reused for earlier conditional forecasts.}
\examples{if (requireNamespace("survival", quietly=TRUE)) {
  # Keep the example small while including a recorded same-day pair.
  d <- survival::colon
  r <- d[d$etype==1 & d$status==1, ]
  terminal <- d[d$etype==2 & d$status==1, ]
  i <- match(r$id, terminal$id)
  tied <- r$id[!is.na(i) & r$time==terminal$time[i]]
  ids <- unique(c(head(unique(d$id), 2), head(tied, 1)))
  f <- cs_colon_fibre(d[!is.na(match(d$id, ids)), ], nodes=2)
  f$metadata[c("patients","recurrences","source_deaths","same_day_ids")]
}}
\seealso{\code{\link{cs_colon}}, \code{\link{cs_coarsen_fibre}}}
