\name{cs_validate_env}
\alias{cs_validate_env}
\alias{print.coherentSurv_environment}
\title{Offline Environment and Analytic Diagnostics}
\usage{
cs_validate_env(timings = FALSE, format = c("list", "json", "text"),
  file = NULL)
\method{print}{coherentSurv_environment}(x, ...)
}
\arguments{
\item{timings}{Include single-run fixture elapsed times and returned-object sizes.}
\item{format}{Return a structured list, JSON string or text lines.}
\item{file}{Optional explicit output path. List output is saved as RDS; JSON and
text are written as UTF-8. NULL writes no file.}
\item{x}{An environment diagnostic report.}
\item{...}{Additional print arguments; currently unused.}
}
\description{Reads installed runtime metadata and executes four bounded,
deterministic numerical fixtures using the maintained R/C++ engine. No patient
data, network access, package installations or external compiler commands are used.}
\details{The report records the loaded package version, installed DESCRIPTION
build metadata when available, required dependency versions, registered native
entry points, current native-file SHA-256 and session library versions. The
DESCRIPTION Built field describes R/platform/date; it does not identify the
compiler that produced the package. Actual package compiler identity and flags
are reported as unavailable when not recorded; a compiler on PATH is not used
to reconstruct historical build provenance.

The operational fixtures check a two-unit Poisson-boundary likelihood, scores
and information; dense and sparse two-state absorbing probabilities and
integrated survival; and a two-path mixture with Louis missing-information
correction. Each has a fixed analytic reference and absolute tolerance 1e-12.
A fixture error, warning, nonfinite result or schema mismatch is a failed check.
The default does not use random numbers or modify RNG state.

Statuses are passed, failed or unavailable. Overall status is failed if any
check fails, otherwise unavailable if provenance is missing, otherwise passed.
Operational status is reported separately so missing compiler provenance does
not obscure successful numerical checks. Optional elapsed times are descriptive
observations for these tiny fixtures, and object sizes exclude peak process
memory, allocator overhead and other runtime objects.

This is an unsigned, mutable diagnostic report, not IQ/OQ/PQ certification,
signed qualification, 21 CFR Part 11 compliance, a machine speed benchmark or
proof of clinical validity. Hashes describe current on-disk bytes and do not
authenticate historical build provenance or an in-memory binary. Environment
paths may identify local installations.}
\value{By default, a coherentSurv_environment list with installation metadata,
session versions, per-check statuses, analytic references, observed fixture
outputs and interpretation. Other formats serialize the same report.}
\seealso{\code{cs_audit_report}, \code{cs_loglik}, \code{cs_propagate},
\code{cs_propagate_sparse}, \code{cs_coarsen}}
\examples{
report <- cs_validate_env()
report$checks[, c("check", "status")]
report$operational_status
}
