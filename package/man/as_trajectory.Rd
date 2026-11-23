\name{as_trajectory}
\alias{as_trajectory}
\title{Validated event histories}
\description{Validate time ordering, eligibility and observation metadata without silently removing records.}
\usage{
as_trajectory(data, unit = "id", start = "start", stop = "stop", 
    event_type = "event_type", terminal = "terminal", 
    observed = NULL, event_types = NULL, time_unit, 
    observation = c("exact", "coarsened"), terminal_cause = NULL) 
}
\arguments{
\item{data}{A data frame of complete counting-process intervals; one row is one observed interval.}
\item{unit}{Name of the unit identifier column, or a selected unit identifier for plotting.}
\item{start}{Name of the interval start column.}
\item{stop}{Name of the interval end column.}
\item{event_type}{Name of the nonterminal event-label column. Zero or 'none' denotes no event.}
\item{terminal}{Name of the binary absorbing-event column.}
\item{observed}{Optional name of a binary monitoring indicator column. Unobserved intervals cannot be fitted as observed zeros.}
\item{event_types}{Declared nonterminal event labels, including types with no observed events.}
\item{time_unit}{Explicit nonempty unit of follow-up time, such as years or operating hours.}
\item{observation}{Exact complete histories or an explicitly enumerated coarsened-path experiment.}
\item{terminal_cause}{Optional name of a column containing the absorbing cause on terminal-event rows only. Nonterminal rows use missing, empty, zero or 'none'; observed cause labels are metadata and do not select a model.}
}
\details{Overlaps, exact simultaneous recurrence and terminal events, zero-length intervals and post-terminal records are rejected. A genuine compound jump requires a different event mark; rounded ties require an explicit coarsening model. The regime attribute summarizes observed recurrent and terminal marks; the observation_process attribute distinguishes exact, coarsened_complete and delayed_entry declarations. Additional fields retain gaps and unmonitored intervals. These tags do not infer scientific eligibility, a missingness mechanism or a fitted likelihood. Competing absorbing causes require an explicit terminal_cause column and are separate from recurrent event marks.}
\value{An rc_trajectory data frame, also inheriting surv_trajectory for compatibility, with schema and observation attributes.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
\examples{
d <- data.frame(id=c(1,1),start=c(0,1),stop=c(1,2),event_type=c(1,0),terminal=c(0,1))
z <- as_trajectory(d,event_types='1',time_unit='years')
build_history(z)
}
