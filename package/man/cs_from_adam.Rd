\name{cs_from_adam}
\alias{cs_from_adam}
\title{Map a selected ADTTE parameter to single-endpoint trajectories}
\description{Validate a deliberately restricted ADaM-style ADTTE extract with
explicit duration, censoring-code and ascertainment declarations. This is a
selected-extract mapping, not a complete or certified CDISC parser.}
\usage{
cs_from_adam(adtte, paramcd, event_codes, censor_codes, time_unit,
  duration_definition, ascertainment, provenance,
  adsl = NULL, covariates = character())
}
\arguments{
\item{adtte}{Data frame containing \code{USUBJID}, \code{PARAMCD}, \code{AVAL}
and \code{CNSR}. There must be exactly one row per subject after parameter
selection. Optional \code{AVALU} values must match \code{time_unit}, ignoring case.}
\item{paramcd}{One explicitly selected \code{PARAMCD}; other parameters are not pooled.}
\item{event_codes, censor_codes}{Explicit nonempty, unique, disjoint vectors
of source \code{CNSR} codes. Their union must cover every selected code;
unobserved declared codes are permitted. Codes are matched as character strings.}
\item{time_unit}{One declared duration unit. No conversion is inferred.}
\item{duration_definition}{Nonempty text describing the source derivation of
\code{AVAL}, including origin and any inclusive-day convention. Positive finite
\code{AVAL} is retained exactly; dates are not used to rederive it.}
\item{ascertainment}{Nonempty analyst declaration describing endpoint and
censoring ascertainment and complete single-endpoint follow-up from the stated
process origin. Its scientific validity is not established by this function.}
\item{provenance}{Nonempty source/version/derivation description to retain in
mapping metadata. Avoid unnecessary identifying information.}
\item{adsl}{Optional ADSL-style data frame with unique \code{USUBJID}, covering
every selected subject. Extra ADSL subjects are permitted.}
\item{covariates}{Distinct ADSL columns to copy by subject identifier. These
must be appropriate baseline covariates; the adapter cannot verify their timing.
Missing values are preserved without imputation. Reserved trajectory and source
column names cannot be overwritten.}
}
\details{
The output has one interval from zero to \code{AVAL} for each selected subject.
The terminal indicator is membership in \code{event_codes}, never
\code{1 - CNSR}. A source parameter representing a different endpoint or a
reverse follow-up analysis must be selected separately with its own definition.
This adapter treats the selected event as absorbing for that single-endpoint
analysis; it does not assert that the endpoint is death or create recurrent marks.

Use character source identifiers to preserve leading zeros and large integer
identifiers. Ordinary numeric integer keys up to \eqn{2^{53}} are accepted, but
precision or zeros already lost during import cannot be recovered. Inputs are
not modified. Source parameter, censor code and original ADTTE row are retained.

The attached \code{cs_mapping} attribute records declarations and the ADSL
covariate selection; ordinary \code{cs_schema} observation metadata remain in
place. Mapping does not establish independent censoring, event completeness,
causal identification or regulatory qualification. Zero-duration observations,
multiple analysis rows per selected subject, recurrent histories and delayed
entry are outside this mapping. Fitting interfaces retain their own assumptions.
}
\value{A \code{cs_trajectory} data frame (also inheriting
\code{surv_trajectory}) with a \code{cs_mapping} metadata attribute.}
\seealso{\code{\link{cs_as_trajectory}}, \code{\link{cs_from_omop}}}
\examples{
# Wholly synthetic; no clinical code meaning is inferred.
d <- data.frame(USUBJID = c("001", "002", "003"), PARAMCD = "OS",
                AVAL = c(5, 10, 7), CNSR = c(0, 1, 2))
z <- cs_from_adam(d, paramcd = "OS", event_codes = 0,
  censor_codes = c(1, 2), time_unit = "days",
  duration_definition = "Synthetic elapsed days from assigned study origin",
  ascertainment = "Complete synthetic endpoint ascertainment until exit",
  provenance = "Synthetic documentation fixture")
z$terminal
attr(z, "cs_mapping")$paramcd
}
