\name{rc_coarsen_fibre}
\alias{rc_coarsen_fibre}
\title{Automatic legal histories for coarsened event times}
\usage{rc_coarsen_fibre(data, unit="id", lower="event_time_lower",
  upper="event_time_upper", start="start", stop="stop", event_type="event_type",
  terminal="terminal", followup="followup", event_types=NULL, time_unit,
  rounding_unit=NULL, nonrecurrent=character(), breaks=c(0,Inf), nodes=4L,
  max_perms=120L, draws=2000L, seed=20261122L, max_paths=100000L,
  format=c("auto","windows","counting"), covariates=NULL)}
\arguments{
\item{data}{Observed event records. Event counts and marks must be complete; times and ordering may be coarsened. No patient rows are silently removed.}
\item{unit,lower,upper,start,stop,event_type,terminal,followup}{Column names. Window input uses lower, upper and a common followup endpoint per subject. Counting input uses start and stop, with the maximum stop as followup. Mark zero denotes no recurrent event; terminal is binary. A row can record both a recurrent and terminal mark, which are expanded as two events.}
\item{event_types}{Declared recurrent labels, or NULL to infer the observed nonzero labels.}
\item{time_unit}{Required unit of all numeric times.}
\item{rounding_unit}{Optional positive recording resolution in the same time unit. Counting event dates t become (max(0,t-resolution),t]; point windows are expanded similarly. It is never inferred from the data. Without it, equal endpoints are exact event times and incompatible exact ties are rejected.}
\item{nonrecurrent}{Marks that can occur at most once.}
\item{breaks}{Baseline discontinuities, increasing from zero to infinity. Quadrature splits there and at all window endpoints.}
\item{nodes}{Gauss-Legendre nodes per continuous integration segment.}
\item{max_perms}{Enumeration budget, default 120 canonical orders per subject.}
\item{draws}{Attempted importance samples per subject when enumeration or quadrature exceeds its budget. Rejected attempts contribute zero; they remain in the Monte Carlo denominator.}
\item{seed}{Nonnegative deterministic seed. A local C++ generator leaves the R random-number stream unchanged.}
\item{max_paths}{Maximum generated paths per subject. A quadrature budget overflow triggers importance sampling, whose draws must fit this budget.}
\item{format}{Explicit windows, counting rows, or automatic detection by available columns.}
\item{covariates}{Subject-constant covariate columns; NULL selects all non-observation columns. Time-varying covariates require a separate known schedule and are not silently assigned to alternative event orders.}
}
\description{Enumerates legal event orders and integrates over their ordered time regions in a C++17 engine, or uses a disclosed importance-sampling approximation when the declared budget is exceeded.}
\value{A \code{recurcast_fibre} list containing expanded counting-process \code{data}, \code{coarsening=list(owner,log_weights)}, event definitions, baseline breaks, a per-subject \code{tie_audit}, attempted-draw identifiers and metadata. Column names for the expanded unit/start/stop/mark/terminal fields match the supplied names. Subject identifiers in expanded data are path identifiers; \code{.rc_source_unit} retains provenance.}
\details{The terminal event is last. Same-mark overlapping windows use a canonical earliest-deadline matching, so an identical ordered history is not counted more than once merely because several assignments to record windows are possible. Identical same-mark windows therefore produce an ordered simplex rather than duplicated factorial copies. Distinct marks can have different legal orders.

Exact point events carry unit integration mass. Continuous windows use nested Gaussian quadrature; enumeration is complete below the cap but the resulting integral remains a numerical approximation. The large-fibre fallback samples event choices and times sequentially, with inverse joint proposal density divided by attempted draws. Its likelihood, scores and Louis information are Monte Carlo approximations, not exact inference. The audit's volume MCSE concerns a constant integrand. Use \code{rc_fibre_mcse} with complete-path log likelihoods for conditional likelihood integration error.

This observation model does not infer how many unrecorded events occurred. Arbitrary gaps with unknown counts require an additional observation model. Forecasting must regenerate a fibre from information available through the landmark, not use future event records. At most 64 recorded events per subject are supported by this enumerator.}
\examples{
observed <- data.frame(id=c(1,1), event_time_lower=c(1,1),
  event_time_upper=c(2,2), event_type=c("1","0"),
  terminal=c(0,1), followup=2, treatment=0)
f <- rc_coarsen_fibre(observed, time_unit="years", event_types="1")
f$tie_audit
sum(exp(f$coarsening$log_weights)) # ordered triangle has area 1/2
}
\seealso{\code{\link{rc_coarsen}}, \code{\link{rc_fibre_mcse}}, \code{\link{rc_colon_fibre}}}
