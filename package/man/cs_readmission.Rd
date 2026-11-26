\name{cs_readmission}
\alias{cs_readmission}
\title{Readmission trajectories from the public source package}
\usage{cs_readmission(data = NULL, time_unit = c("years", "days"))}
\arguments{
\item{data}{Optional canonical source data frame. NULL loads the dataset from the suggested package \pkg{frailtypack}.}
\item{time_unit}{Years (365.25 days per year) or days.}
}
\description{Constructs exact counting-process trajectories without bundling clinical patient records in this package.}
\value{A list with \code{data}, \code{coarsening=NULL}, a suggested baseline \code{formula}, \code{event_types}, \code{nonrecurrent}, and provenance/count \code{metadata}. The output data include both \code{id} and \code{unit}; the default fit identifier is \code{id}.}
\details{The canonical data contain 403 patients, 861 interval rows, 458 recorded rehospitalizations, and 109 terminal deaths. Interval rows are not event counts. Sex, chemotherapy, Dukes stage, and Charlson indicators are fixed at the patient's first available interval. Later Charlson values remain in \code{recorded_charlson} but are not carried backward. This single recurrent mark represents rehospitalization, not distinct clinical diseases. Rehospitalization and death are prognostic endpoints; fitted contrasts are not treatment effects.}
\examples{
if (requireNamespace("frailtypack", quietly = TRUE)) {
  x <- cs_readmission()
  x$metadata[c("patients", "recurrences", "terminal_events")]
  # Fit explicitly when desired:
  # fit <- cs_fit(x$formula, x$data, time_unit="years", event_types="1")
}
}
\seealso{\code{\link{cs_colon}}, \code{\link{cs_fit}}}
\references{The source-package data dictionary is in the frailtypack reference manual, entry readmission: \url{https://CRAN.R-project.org/package=frailtypack}.}
