\name{rc_audit_report}
\alias{rc_audit_report}
\title{Export Fit Provenance and Numerical Diagnostics}
\description{Create a list, JSON or text audit from a fitted coherent4cast object, validated trajectory or prediction table. Subject records are not included. Fit context and export context remain distinct.}
\usage{
rc_audit_report(object, predictions = NULL, source_files = NULL,
  format = c("list", "json", "text"), file = NULL,
  include_session = TRUE)
}
\arguments{
\item{object}{A fitted \code{coherent4cast}, \code{coherent4cast_competing}, validated trajectory or \code{coherent4cast_prediction} object.}
\item{predictions}{Optional prediction data frame. Recognized numerical error columns are summarized without exporting unit identifiers or individual predictions.}
\item{source_files}{Explicit paths to existing source files to hash at export time. Directories are rejected and no directory or repository is scanned automatically.}
\item{format}{Return a nested R list, a JSON character string or readable text lines.}
\item{file}{Optional output path. List output is saved as RDS serialization version 3; JSON and text are written as UTF-8.}
\item{include_session}{Include captured fit session information, when available, and current export session information.}
}
\details{
New fitting workflows capture their context using a UTC timestamp, session information,
the supplied call with data-bearing expressions redacted, and a SHA-256 fingerprint of
\code{serialize(data, NULL, version = 3)}. SHA-256 uses \pkg{digest}; JSON encoding uses
\pkg{jsonlite}. The fingerprint covers serialized bytes and attributes. It does not
assert semantic equality after reordering columns or changing serialization/platform
details, and hashing is not an anonymization guarantee.

Legacy \code{audit$data_hash} values retain their recorded MD5 label. They are never
renamed SHA-256. For an old object without fit-time context, the historical timestamp
and session remain unavailable. A fresh fingerprint of the object's current data,
the current session and the export timestamp are recorded separately under
\code{export}; they cannot establish historical provenance.

Source-file hashes refer only to files explicitly supplied by the caller, at the
recorded export time. They do not establish that these files produced the fitted
object. Calls redact data, weights, subsets, coarsening and other data-bearing
arguments. Schema and model metadata may still contain analyst-supplied labels or
expressions, which should be reviewed before sharing.

Numerical bounds, quadrature diagnostics and covariate Monte Carlo errors remain
distinct from sampling uncertainty. The report preserves failed convergence flags.
It does not certify model assumptions, causal identification, coverage, safety,
regulatory acceptance or publication readiness. JSON null encodes missing or
nonfinite numerical values; inspect diagnostic flags and nonfinite counts.
}
\value{The selected list, JSON string or text lines, also written to \code{file} when supplied. List fields include \code{fit_context}, \code{model}, \code{convergence}, \code{numerical}, \code{source_files} and \code{export}.}
\examples{
d <- data.frame(id = 1:2, start = 0, stop = c(1, 2),
  event_type = "0", terminal = c(1, 0))
trajectory <- as_trajectory(d, time_unit = "years")
report <- rc_audit_report(trajectory, include_session = FALSE)
report$export$data_fingerprint$algorithm
report$fit_context$status
}
\seealso{\code{as_trajectory}, \code{coherent4cast}, \code{rc_surv}, \code{rc_competing}}
