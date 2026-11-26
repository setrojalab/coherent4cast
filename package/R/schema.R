.cs_error <- function(code, message) {
  stop(structure(list(message = paste0(code, ": ", message), call = NULL, code = code),
                 class = c(paste0("cs_error_", gsub("[^A-Za-z0-9_]", "_", code)),
                           "coherentSurv_error", "error", "condition")))
}

.cs_warning <- function(code, message) {
  warning(structure(list(message = message, call = NULL, code = code),
                    class = c(paste0("cs_warning_", gsub("[^A-Za-z0-9_]", "_", code)),
                              "coherentSurv_warning", "warning", "condition")))
}

#' Construct validated counting-process trajectories
#' @export
cs_as_trajectory <- function(data, unit = "id", start = "start", stop = "stop",
                          event_type = "event_type", terminal = "terminal",
                          observed = NULL, event_types = NULL, time_unit,
                          observation = c("exact", "coarsened"), terminal_cause = NULL) {
  observation <- match.arg(observation)
  if (missing(time_unit) || length(time_unit) != 1L || is.na(time_unit) || !nzchar(time_unit))
    .cs_error("time_unit", "Declare one nonempty time unit.")
  if (!is.data.frame(data) || !nrow(data)) .cs_error("data", "data must be a nonempty data frame.")
  # Revalidation of an already declared trajectory retains its explicit cause
  # column. Raw nonterminal marks are never reinterpreted as absorbing causes.
  if (is.null(terminal_cause) && inherits(data, "surv_trajectory"))
    terminal_cause <- attr(data, "cs_schema")$terminal_cause
  if (!is.null(terminal_cause) && (!is.character(terminal_cause) || length(terminal_cause) != 1L ||
      is.na(terminal_cause) || !nzchar(terminal_cause)))
    .cs_error("terminal_cause", "terminal_cause must be NULL or one absorbing-cause column name.")
  required <- c(unit, start, stop, event_type, terminal, observed, terminal_cause)
  if (!all(required %in% names(data))) .cs_error("columns", "A declared data column is absent.")
  d <- data
  if (anyNA(d[[unit]]) || any(!nzchar(as.character(d[[unit]])))) .cs_error("unit", "Unit identifiers cannot be missing or empty.")
  if (!is.numeric(d[[start]]) || !is.numeric(d[[stop]]) || any(!is.finite(d[[start]])) || any(!is.finite(d[[stop]])))
    .cs_error("time", "Start and stop times must be finite numeric values.")
  if (any(d[[start]] < 0) || any(d[[stop]] <= d[[start]])) .cs_error("interval", "Require 0 <= start < stop; zero-length ties need a coarsening model.")
  if (anyNA(d[[terminal]]) || !all(d[[terminal]] %in% c(0, 1, FALSE, TRUE))) .cs_error("terminal", "terminal must contain only 0 and 1.")
  raw <- as.character(d[[event_type]])
  if (anyNA(raw)) .cs_error("event_type", "Event types cannot be missing.")
  none <- raw %in% c("0", "none", "")
  if (is.null(event_types)) event_types <- sort(unique(raw[!none]))
  event_types <- as.character(event_types)
  if (anyNA(event_types) || anyDuplicated(event_types) || any(event_types %in% c("0", "none", ""))) .cs_error("event_types", "Declare distinct nonzero nonterminal event labels.")
  if (any(!none & !raw %in% event_types)) .cs_error("event_type", "An event is outside the declared event types.")
  mark <- match(raw, event_types); mark[none] <- 0L
  if (any(mark > 0 & d[[terminal]] == 1)) .cs_error("tie", "One exact row cannot contain both recurrence and terminal event. Supply legal coarsened paths.")
  d$.cs_unit <- match(as.character(d[[unit]]), unique(as.character(d[[unit]])))
  d$.cs_start <- d[[start]]; d$.cs_stop <- d[[stop]]
  d$.cs_mark <- mark; d$.cs_terminal <- as.integer(d[[terminal]])
  causes <- character()
  if (!is.null(terminal_cause)) {
    cause <- as.character(d[[terminal_cause]])
    absent <- is.na(cause) | cause %in% c("", "0", "none")
    if (any(d$.cs_terminal == 1L & absent))
      .cs_error("terminal_cause", "Every absorbing-event row requires a nonmissing, nonzero cause label.")
    if (any(d$.cs_terminal == 0L & !absent))
      .cs_error("terminal_cause", "A cause label is permitted only on an absorbing-event row; recurrent marks are separate.")
    causes <- sort(unique(cause[d$.cs_terminal == 1L]))
  }
  if (!is.null(observed) && (anyNA(d[[observed]]) || !all(d[[observed]] %in% c(0, 1, FALSE, TRUE))))
    .cs_error("observation", "Observation indicators must contain only 0 and 1.")
  d$.cs_observed <- if (is.null(observed)) rep(TRUE, nrow(d)) else as.logical(d[[observed]])
  if (anyNA(d$.cs_observed)) .cs_error("observation", "Observation indicators cannot be missing.")
  ordering <- order(d$.cs_unit, d$.cs_start, d$.cs_stop)
  if (!identical(ordering, seq_len(nrow(d)))) .cs_error("ordering", "Rows must be grouped by unit and ordered in time. Sort explicitly before fitting.")
  gaps <- integer(); entry <- integer()
  for (rows in split(seq_len(nrow(d)), d$.cs_unit)) {
    if (d$.cs_start[rows[1L]] > 0) entry <- c(entry, rows[1L])
    if (length(rows) > 1L) {
      delta <- d$.cs_start[rows[-1L]] - d$.cs_stop[rows[-length(rows)]]
      if (any(delta < -1e-12)) .cs_error("overlap", "Intervals overlap within a unit.")
      gaps <- c(gaps, rows[-1L][delta > 1e-12])
    }
    terminal_rows <- which(d$.cs_terminal[rows] == 1L)
    if (length(terminal_rows) > 1L || (length(terminal_rows) && terminal_rows != length(rows))) .cs_error("post_terminal", "A unit has records after its terminal event.")
  }
  recurrence_counts <- as.integer(rowsum(as.integer(d$.cs_mark > 0L), d$.cs_unit))
  observed_pattern <- if (!is.null(terminal_cause) && length(causes) > 1L) "competing_terminal_causes_observed" else
    if (any(recurrence_counts > 1L)) "recurrent_nonterminal_events_observed" else
    if (any(recurrence_counts > 0L)) "at_most_one_nonterminal_event_per_unit_observed" else
    if (any(d$.cs_terminal > 0L)) "terminal_events_only_observed" else "no_events_observed"
  regime <- if (!is.null(terminal_cause) && length(causes) > 1L) "competing_risks" else
    if (any(recurrence_counts > 0L) && any(d$.cs_terminal > 0L)) "joint_recurrent_terminal" else
    if (any(recurrence_counts > 0L)) "recurrent_only" else "standard_survival"
  regime_details <- list(observed_pattern = observed_pattern,
    nonterminal_types_declared = event_types,
    nonterminal_types_observed = event_types[sort(unique(mark[mark > 0L]))],
    maximum_observed_nonterminal_events_per_unit = max(recurrence_counts),
    units_with_multiple_nonterminal_events = sum(recurrence_counts > 1L),
    terminal_events_observed = sum(d$.cs_terminal),
    terminal_cause_column = terminal_cause, terminal_causes_observed = causes,
    model_inferred = FALSE,
    scope = paste("Observed event counts do not determine whether recurrence is scientifically possible.",
      "Competing absorbing causes are described only when terminal_cause was explicitly declared and multiple causes occur.",
      "These labels do not select a likelihood, event-eligibility rule or causal model."))
  observation_process <- if (length(entry)) "delayed_entry" else
    if (observation == "coarsened") "coarsened_complete" else "exact"
  observation_process_details <- list(declaration = observation, monitoring_column = observed,
    unmonitored_rows = which(!d$.cs_observed), gap_rows = gaps, delayed_entry_rows = entry,
    continuous_observation_from_zero = !length(gaps) && !length(entry) && all(d$.cs_observed),
    missingness_mechanism = "not identified by schema validation",
    scope = paste("Exact/coarsened is an analyst declaration, not a verified MAR/CAR/MNAR assumption.",
      "An observed zero, an unmonitored interval, a between-row gap and delayed entry are distinct.",
      "Delayed entry is tagged here; fit and prediction interfaces retain their own supported-observation checks."))
  attr(d, "cs_schema") <- list(unit = unit, start = start, stop = stop, event_type = event_type,
       terminal = terminal, observed = observed, event_types = event_types, time_unit = time_unit,
       observation = observation, gaps = gaps, delayed_entry = entry,
       original_rows = nrow(d), dropped_rows = 0L, terminal_cause = terminal_cause,
       terminal_causes_observed = causes, regime = regime, regime_details = regime_details)
  attr(d, "regime") <- regime
  attr(d, "observation_process") <- observation_process
  attr(d, "observation_process_details") <- observation_process_details
  class(d) <- c("cs_trajectory", "surv_trajectory", "data.frame")
  d
}

#' Predictable counts and recency summaries
#' @export
cs_build_history <- function(data, cap = 2L, decay = NULL) {
  if (!inherits(data, "surv_trajectory")) .cs_error("schema", "Use cs_as_trajectory() first.")
  if (length(cap) != 1L || !is.finite(cap) || cap < 0 || cap != as.integer(cap)) .cs_error("cap", "cap must be a nonnegative integer.")
  if (!is.null(decay) && (length(decay) != 1L || !is.finite(decay) || decay < 0)) .cs_error("decay", "decay must be nonnegative and finite.")
  K <- length(attr(data, "cs_schema")$event_types)
  z<-tryCatch(cs_cpp_history(as.integer(data$.cs_unit),data$.cs_start,data$.cs_stop,as.integer(data$.cs_mark),K,cap,if(is.null(decay))-1 else decay),
    error=function(e).cs_error("history",conditionMessage(e)))
  counts<-z$counts;recency<-z$recency;gap<-z$gap
  if(K) {colnames(counts) <- paste0("count_", seq_len(K)); colnames(recency) <- paste0("recency_", seq_len(K)); colnames(gap) <- paste0("gap_", seq_len(K))}
  list(counts = counts, recency = if (is.null(decay)) NULL else recency, gap = gap,
       predictable = TRUE, cap = cap, decay = decay)
}

#' Controls for joint gamma estimation and forecasting
#' @export
cs_control <- function(maxit = 300L, tolerance = 1e-8, history_cap = 2L,
    phi_fixed = NULL, start_phi = 0.5, quadrature = 30L, propagation_tolerance = 1e-12,
    max_states = 50000L, coefficient_bounds = c(-12, 8), log_phi_bounds = c(-5, 2.5),
    quadrature_tolerance=1e-7,quadrature_max=1024L) {
  scalar <- function(x, lower, integer=FALSE) is.numeric(x) && length(x)==1L && is.finite(x) && x>=lower && (!integer || x==floor(x))
  interval <- function(x) is.numeric(x) && length(x)==2L && all(is.finite(x)) && x[1]<x[2]
  if (!scalar(maxit,1,TRUE) || !scalar(tolerance,.Machine$double.xmin) || !scalar(history_cap,0,TRUE) ||
      !scalar(quadrature,2,TRUE) || !scalar(max_states,2,TRUE) || !scalar(propagation_tolerance,.Machine$double.xmin) ||
      !scalar(start_phi,.Machine$double.xmin) || !interval(coefficient_bounds) || !interval(log_phi_bounds)||
      !scalar(quadrature_tolerance,.Machine$double.xmin)||!scalar(quadrature_max,quadrature,TRUE)||quadrature_max>1024)
    .cs_error("control", "Iteration, tolerance, dispersion, state and quadrature settings must be valid.")
  if (!is.null(phi_fixed) && (length(phi_fixed) != 1L || !is.finite(phi_fixed) || phi_fixed < 0)) .cs_error("dispersion", "phi_fixed must be nonnegative or NULL.")
  list(maxit=as.integer(maxit),tolerance=tolerance,history_cap=as.integer(history_cap),phi_fixed=phi_fixed,
       start_phi=start_phi,quadrature=as.integer(quadrature),propagation_tolerance=propagation_tolerance,
       max_states=as.integer(max_states),coefficient_bounds=coefficient_bounds,log_phi_bounds=log_phi_bounds,
       quadrature_tolerance=quadrature_tolerance,quadrature_max=as.integer(quadrature_max))
}
