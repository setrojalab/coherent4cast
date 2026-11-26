# These adapters validate declared, selected extracts. They do not discover
# clinical definitions or establish ascertainment from source-table presence.
.cs_standards_text <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x)))
    .cs_error("mapping_declaration", paste(name, "must be one nonempty character declaration."))
  x
}

.cs_standards_table <- function(x, columns, name, nonempty = TRUE) {
  if (!is.data.frame(x) || (nonempty && !nrow(x)) || anyDuplicated(names(x)) ||
      !all(columns %in% names(x)))
    .cs_error("mapping_columns", paste(name, "must be a data frame with unique column names and required columns:",
      paste(columns, collapse = ", ")))
}

.cs_standards_keys <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) && !is.object(x)) {
    if (any(!is.finite(x)) || any(x != floor(x)) || any(abs(x) > 2^53))
      .cs_error("mapping_key", paste(name, "numeric keys must be finite exactly representable integers; use character keys otherwise."))
  } else if (!is.character(x)) {
    .cs_error("mapping_key", paste(name, "must contain character, factor or ordinary integer-valued numeric keys."))
  }
  x <- as.character(x)
  if (anyNA(x) || any(!nzchar(x)) || any(trimws(x) != x))
    .cs_error("mapping_key", paste(name, "keys must be nonmissing, nonempty and free of outer whitespace."))
  x
}

.cs_standards_dates <- function(x, name) {
  if (!inherits(x, "Date") || anyNA(x) || any(!is.finite(as.numeric(x))) ||
      any(as.numeric(x) != floor(as.numeric(x))))
    .cs_error("mapping_date", paste(name, "must contain nonmissing whole Date values; parse and review source dates explicitly."))
  x
}

#' Map a selected ADTTE parameter with explicit censoring codes
#' @export
cs_from_adam <- function(adtte, paramcd, event_codes, censor_codes, time_unit,
                         duration_definition, ascertainment, provenance,
                         adsl = NULL, covariates = character()) {
  paramcd <- .cs_standards_text(paramcd, "paramcd")
  time_unit <- .cs_standards_text(time_unit, "time_unit")
  duration_definition <- .cs_standards_text(duration_definition, "duration_definition")
  ascertainment <- .cs_standards_text(ascertainment, "ascertainment")
  provenance <- .cs_standards_text(provenance, "provenance")
  .cs_standards_table(adtte, c("USUBJID", "PARAMCD", "AVAL", "CNSR"), "adtte")
  parameters <- .cs_standards_keys(adtte$PARAMCD, "PARAMCD")
  selected_rows <- which(parameters == paramcd)
  if (!length(selected_rows)) .cs_error("mapping_parameter", "The selected PARAMCD has no rows.")
  selected <- as.data.frame(adtte[selected_rows, , drop = FALSE])
  ids <- .cs_standards_keys(selected$USUBJID, "USUBJID")
  if (anyDuplicated(ids)) .cs_error("mapping_subject", "Require exactly one row per USUBJID for the selected PARAMCD.")
  event_codes <- .cs_standards_keys(event_codes, "event_codes")
  censor_codes <- .cs_standards_keys(censor_codes, "censor_codes")
  if (!length(event_codes) || !length(censor_codes) || anyDuplicated(event_codes) ||
      anyDuplicated(censor_codes) || length(intersect(event_codes, censor_codes)))
    .cs_error("mapping_censor", "Declare nonempty, distinct, disjoint event_codes and censor_codes.")
  codes <- .cs_standards_keys(selected$CNSR, "CNSR")
  if (any(!codes %in% c(event_codes, censor_codes)))
    .cs_error("mapping_censor", "Every selected CNSR value must belong to exactly one declared event or censor set.")
  if (!is.numeric(selected$AVAL) || is.object(selected$AVAL) ||
      any(!is.finite(selected$AVAL)) || any(selected$AVAL <= 0))
    .cs_error("mapping_duration", "Selected AVAL values must be positive finite numeric durations; zero-duration observations need another observation model.")
  if ("AVALU" %in% names(selected)) {
    units <- .cs_standards_keys(selected$AVALU, "AVALU")
    if (any(tolower(units) != tolower(time_unit)))
      .cs_error("mapping_duration", "Selected AVALU must match the declared time_unit (ignoring case); no unit conversion is inferred.")
  }
  if (!is.character(covariates) || anyNA(covariates) || any(!nzchar(covariates)) || anyDuplicated(covariates))
    .cs_error("mapping_covariates", "covariates must be distinct nonempty ADSL column names.")
  d <- data.frame(id = ids, start = 0, stop = selected$AVAL, event_type = "0",
    terminal = as.integer(codes %in% event_codes), source_parameter = paramcd,
    source_censor_code = codes, source_row = selected_rows, stringsAsFactors = FALSE)
  if (any(covariates %in% c(names(d), "USUBJID")) || any(startsWith(covariates, ".cs_")))
    .cs_error("mapping_covariates", "ADSL covariate names must not overwrite identifiers, trajectory or source columns.")
  if (is.null(adsl)) {
    if (length(covariates)) .cs_error("mapping_covariates", "Supply adsl when selecting ADSL covariates.")
  } else {
    .cs_standards_table(adsl, c("USUBJID", covariates), "adsl")
    adsl_ids <- .cs_standards_keys(adsl$USUBJID, "ADSL USUBJID")
    if (anyDuplicated(adsl_ids) || any(!ids %in% adsl_ids))
      .cs_error("mapping_subject", "ADSL must have unique USUBJID values and cover every selected ADTTE subject.")
    for (v in covariates) d[[v]] <- adsl[[v]][match(ids, adsl_ids)]
  }
  ans <- cs_as_trajectory(d, event_types = character(), time_unit = time_unit)
  attr(ans, "cs_mapping") <- list(adapter = "cs_from_adam", source = "selected ADTTE parameter",
    paramcd = paramcd, event_codes = event_codes, censor_codes = censor_codes,
    duration = "AVAL preserved without date rederivation or rescaling", duration_definition = duration_definition,
    time_unit = time_unit, ascertainment = ascertainment, provenance = provenance,
    source_rows = selected_rows, adsl_covariates = covariates,
    observation = "Analyst-declared complete single-endpoint follow-up from time zero; censoring codes do not establish independent censoring. Optional ADSL baseline covariates join by exact unique subject keys without imputation.",
    scope = "Validated mapping of a selected extract, not a complete or certified CDISC parser, recurrent-event adapter or causal identification procedure.")
  ans
}

#' Map reviewed OMOP event extracts under a declared exact-date convention
#' @export
cs_from_omop <- function(person, observation_period, events, concept_map, origin,
                         time_unit, timing, ascertainment, provenance, death = NULL,
                         concept_column = "condition_concept_id",
                         date_column = "condition_start_date") {
  time_unit <- .cs_standards_text(time_unit, "time_unit")
  timing <- .cs_standards_text(timing, "timing")
  ascertainment <- .cs_standards_text(ascertainment, "ascertainment")
  provenance <- .cs_standards_text(provenance, "provenance")
  concept_column <- .cs_standards_text(concept_column, "concept_column")
  date_column <- .cs_standards_text(date_column, "date_column")
  if (!time_unit %in% c("days", "years"))
    .cs_error("mapping_duration", "OMOP mapping supports explicitly declared 'days' or 'years' (365.25 days/year).")
  if (timing != "exact_day_boundaries")
    .cs_error("mapping_timing", "Explicitly declare timing='exact_day_boundaries' only when justified. Rounded, imputed or uncertain times require a separate coarsening model.")
  .cs_standards_table(person, "person_id", "person")
  .cs_standards_table(observation_period,
    c("person_id", "observation_period_start_date", "observation_period_end_date"), "observation_period")
  .cs_standards_table(events, c("person_id", concept_column, date_column), "events", nonempty = FALSE)
  ids <- .cs_standards_keys(person$person_id, "person_id")
  op_ids <- .cs_standards_keys(observation_period$person_id, "observation_period person_id")
  if (anyDuplicated(ids) || anyDuplicated(op_ids) || !setequal(ids, op_ids))
    .cs_error("mapping_period", "Require unique persons and exactly one observation period for each and only each supplied person; multiple or overlapping periods are unsupported.")
  origin <- .cs_standards_dates(origin, "origin")
  if (length(origin) != length(ids)) .cs_error("mapping_origin", "Supply one declared origin Date per person row, in that order.")
  op_start <- .cs_standards_dates(observation_period$observation_period_start_date, "observation_period_start_date")[match(ids, op_ids)]
  op_end <- .cs_standards_dates(observation_period$observation_period_end_date, "observation_period_end_date")[match(ids, op_ids)]
  if (any(op_end <= op_start) || any(origin != op_start))
    .cs_error("mapping_origin", "Each declared process origin must equal its observation-period start, with a later end. Delayed entry or discarded pre-origin history is unsupported.")
  if (!is.character(concept_map) || is.null(names(concept_map)) || !length(concept_map) ||
      anyNA(concept_map) || anyNA(names(concept_map)) || anyDuplicated(names(concept_map)) ||
      any(!nzchar(names(concept_map))) || any(trimws(names(concept_map)) != names(concept_map)) ||
      any(trimws(concept_map) != concept_map) || any(concept_map %in% c("", "0", "none")))
    .cs_error("mapping_concepts", "Supply a reviewed named character map from unique source concepts to nonzero recurrent-event labels.")
  event_ids <- .cs_standards_keys(events$person_id, "event person_id")
  concepts <- .cs_standards_keys(events[[concept_column]], "event concept")
  event_dates <- .cs_standards_dates(events[[date_column]], date_column)
  if (any(!event_ids %in% ids)) .cs_error("mapping_subject", "Every event must link to a supplied person.")
  if (any(!concepts %in% names(concept_map)))
    .cs_error("mapping_concepts", "An event concept is absent from the reviewed map; select the intended extract explicitly, without automatic dropping.")
  owner <- match(event_ids, ids)
  if (any(event_dates <= origin[owner] | event_dates > op_end[owner]))
    .cs_error("mapping_containment", "Events must occur after the declared origin and within the observation period.")
  if (anyDuplicated(data.frame(person = event_ids, concept = concepts, date = event_dates)))
    .cs_error("mapping_duplicate", "Duplicate person/concept/date records require a reviewed upstream episode rule; no automatic collapse is performed.")
  if (anyDuplicated(data.frame(person = event_ids, date = event_dates)))
    .cs_error("mapping_tie", "Multiple events on one person's date require a separate coarsening model; no ordering or jitter is invented.")
  exit <- op_end
  terminal <- rep(0L, length(ids))
  if (!is.null(death)) {
    .cs_standards_table(death, c("person_id", "death_date"), "death", nonempty = FALSE)
    death_ids <- .cs_standards_keys(death$person_id, "death person_id")
    death_dates <- .cs_standards_dates(death$death_date, "death_date")
    if (anyDuplicated(death_ids) || any(!death_ids %in% ids))
      .cs_error("mapping_subject", "Death records must link uniquely to supplied persons.")
    death_owner <- match(death_ids, ids)
    if (any(death_dates <= origin[death_owner] | death_dates > op_end[death_owner]))
      .cs_error("mapping_containment", "Every supplied death must occur after origin and within the declared observation period.")
    exit[death_owner] <- death_dates
    terminal[death_owner] <- 1L
  }
  if (any(event_dates >= exit[owner]))
    .cs_error("mapping_absorption", "Events must precede death or administrative exit strictly; terminal/exit ties and post-terminal events require another declared observation model.")
  divisor <- if (time_unit == "years") 365.25 else 1
  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    ix <- which(owner == i)
    ix <- ix[order(event_dates[ix])]
    times <- as.numeric(event_dates[ix] - origin[i]) / divisor
    end <- as.numeric(exit[i] - origin[i]) / divisor
    rows[[i]] <- data.frame(id = ids[i], start = c(0, times), stop = c(times, end),
      event_type = c(unname(concept_map[concepts[ix]]), "0"),
      terminal = c(rep(0L, length(ix)), terminal[i]), observed = TRUE,
      source_concept = c(concepts[ix], NA_character_), source_row = c(ix, NA_integer_),
      stringsAsFactors = FALSE)
  }
  ans <- cs_as_trajectory(do.call(rbind, rows), event_types = unique(unname(concept_map)),
    observed = "observed", time_unit = time_unit)
  attr(ans, "cs_mapping") <- list(adapter = "cs_from_omop", source = "reviewed selected OMOP-style event extract",
    concept_map = concept_map, concept_column = concept_column, date_column = date_column,
    origin = data.frame(person_id = ids, origin = origin, exit = exit), time_unit = time_unit,
    days_per_year = if (time_unit == "years") 365.25 else NA_real_, timing = timing,
    ascertainment = ascertainment, provenance = provenance,
    episode_rule = "One supplied row per reviewed episode; duplicate and same-day events rejected, never collapsed or jittered.",
    observation = "Analyst declares complete event and terminal ascertainment continuously from the process origin, which must equal observation-period start. Follow-up ends at supplied death or observation end; elapsed dates have no inclusive-day adjustment. Database observation periods do not establish complete ascertainment.",
    scope = "Selected-extract mapping, not a complete or certified OMOP CDM parser. No concept discovery, hidden-count reconstruction, intermittent gaps, delayed entry, date imputation or causal identification.")
  ans
}
