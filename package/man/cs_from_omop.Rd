\name{cs_from_omop}
\alias{cs_from_omop}
\title{Map a reviewed OMOP event extract with explicit observation declarations}
\description{Map selected event records to contiguous trajectories under a
restricted complete-observation and exact-date convention. No clinical concept
set is supplied or discovered, and this is not a full or certified OMOP parser.}
\usage{
cs_from_omop(person, observation_period, events, concept_map, origin,
  time_unit, timing, ascertainment, provenance, death = NULL,
  concept_column = "condition_concept_id",
  date_column = "condition_start_date")
}
\arguments{
\item{person}{Data frame containing unique \code{person_id} values. Other
columns are not copied or interpreted.}
\item{observation_period}{Data frame with \code{person_id},
\code{observation_period_start_date} and \code{observation_period_end_date}.
Exactly one period for each and only each supplied person is required.
Dates must be nonmissing whole \code{Date} objects with start before end.}
\item{events}{Selected, reviewed episode data frame with \code{person_id} and
the declared concept/date columns. All rows must link to supplied persons and
map to declared concepts. Zero rows are allowed with correctly typed columns.
No unknown concept, duplicate or post-exit row is silently dropped.}
\item{concept_map}{Nonempty named character vector mapping unique source
concept values to recurrent-event labels. Labels cannot be \code{"0"},
\code{"none"} or empty. Several concepts may map to one label. Concept selection,
clinical review, eligibility and episode derivation remain the caller's work.}
\item{origin}{One whole \code{Date} per row of \code{person}, in that order,
explicitly declaring the process origins. Each must equal that person's
observation-period start. The adapter does not reset a previously running
process to a later observation start or support delayed entry.}
\item{time_unit}{Explicitly \code{"days"} or \code{"years"}; years use
365.25 days. Elapsed times are date differences, without adding an inclusive day.}
\item{timing}{Must be explicitly \code{"exact_day_boundaries"}: the analyst
asserts that treating the supplied dates as exact boundaries is appropriate.
This does not resolve rounded, imputed or uncertain event times. Such data need
a separately justified coarsening model.}
\item{ascertainment}{Nonempty analyst declaration describing complete
recurrent and terminal ascertainment, continuous monitoring from each process
origin, and the upstream reviewed episode rule. Table presence does not verify
these properties.}
\item{provenance}{Nonempty source/version/extraction and concept-review
description retained in mapping metadata. Avoid unnecessary identifying information.}
\item{death}{Optional data frame with \code{person_id} and whole
\code{Date} \code{death_date}, at most one record per person. All supplied deaths
must occur after origin and within the observation period. With \code{NULL},
the caller declares no terminal events during these intervals; absence of a
database death table alone is not evidence for that declaration.}
\item{concept_column, date_column}{Names of the event concept and event date
columns. Defaults target a selected \code{condition_occurrence}-style extract;
other domains require the caller's explicit selection and column declarations.}
}
\details{
Every event must occur strictly after origin and strictly before the earlier
of death and observation end. The final interval ends at that exit and is
terminal only for a supplied in-period death. Events on an administrative exit
date are also rejected by this deliberately narrow implementation.

Duplicate person/concept/date records are rejected; they are not assumed to
represent one clinical episode. Multiple events on the same person's date,
including different concepts or marks, and event/death ties are rejected without
jitter or invented ordering. Derive episodes using a reviewed upstream rule or
use a separate coarsening workflow. Multiple, overlapping or gapped observation
periods, unknown hidden-event counts, delayed entry and date imputation are not
supported. Even valid OMOP data may lie outside this mapping's scope.

Use character keys to preserve leading zeros and large identifiers. Ordinary
numeric integer keys up to \eqn{2^{53}} are accepted, but lost source precision
cannot be recovered. Event rows are ordered within each person; person order
follows \code{person}. Original event row numbers and source concepts are
retained. Inputs are not modified. No ADSL-like baseline covariate join or
vocabulary download is performed.

The \code{cs_mapping} attribute retains the reviewed map, origins, exits and
declarations. Existing \code{cs_schema} and observation metadata still apply.
Complete ascertainment, exact timing, censoring assumptions and clinical coding
remain analyst responsibilities. This mapping does not establish causal
identification, database conformance or regulatory qualification.
}
\value{A \code{cs_trajectory} data frame (also inheriting
\code{surv_trajectory}) with \code{cs_mapping} metadata.}
\references{OMOP Common Data Model v5.4 table documentation:
\url{https://ohdsi.github.io/CommonDataModel/cdm54.html}. The adapter supports
only the declared extract contract above, not the full standard.}
\seealso{\code{\link{cs_as_trajectory}}, \code{\link{cs_from_adam}}}
\examples{
# Synthetic local codes; these are not real clinical concepts.
p <- data.frame(person_id = c("001", "002"))
origin <- rep(as.Date("2020-01-01"), 2)
op <- data.frame(person_id = p$person_id,
  observation_period_start_date = origin,
  observation_period_end_date = origin + 10)
e <- data.frame(person_id = "001", condition_source_value = "EXAMPLE_A",
  condition_start_date = origin[1] + 2)
z <- cs_from_omop(p, op, e, concept_map = c(EXAMPLE_A = "A"),
  origin = origin, time_unit = "days", timing = "exact_day_boundaries",
  ascertainment = "Complete synthetic follow-up; one row per constructed episode",
  provenance = "Synthetic documentation fixture",
  concept_column = "condition_source_value")
z[, c("id", "start", "stop", "event_type", "terminal")]
}
