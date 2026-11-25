.rc_source_data <- function(name, package, supplied) {
  if (!is.null(supplied)) {
    if (!is.data.frame(supplied) || !nrow(supplied))
      .rc_error("adapter_data", "Supply a nonempty data frame or NULL to load the canonical package data.")
    return(supplied)
  }
  if (!requireNamespace(package, quietly = TRUE))
    .rc_error("optional_package", paste0("Install the suggested package '", package, "' to load ", name, ", or supply its data frame explicitly."))
  # survival exports colon from its lazy-data database; recent versions do not
  # provide a separately indexed data("colon") entry.
  exported <- tryCatch(getExportedValue(package, name), error = function(e) NULL)
  if (is.data.frame(exported)) return(exported)
  envir <- new.env(parent = baseenv())
  utils::data(list = name, package = package, envir = envir)
  if (!exists(name, envir = envir, inherits = FALSE))
    .rc_error("adapter_data", "The named dataset was not found in its source package.")
  get(name, envir = envir, inherits = FALSE)
}

.rc_adapter_columns <- function(data, columns) {
  if (!all(columns %in% names(data)))
    .rc_error("adapter_columns", paste("Missing source columns:", paste(setdiff(columns, names(data)), collapse = ", ")))
}

#' Load the public colorectal readmission data as counting-process trajectories
#' @export
rc_readmission <- function(data = NULL, time_unit = c("years", "days")) {
  time_unit <- match.arg(time_unit)
  supplied <- !is.null(data)
  raw <- .rc_source_data("readmission", "frailtypack", data)
  .rc_adapter_columns(raw, c("id", "t.start", "t.stop", "event", "death", "sex", "chemo", "dukes", "charlson"))
  if (anyNA(raw$id)) .rc_error("adapter_id", "Source patient identifiers cannot be missing.")
  coding <- list(sex = c("Female", "Male"), chemo = c("Treated", "NonTreated"),
    dukes = c("A-B", "C", "D"), charlson = c("0", "1-2", "3"))
  for (v in names(coding)) if (anyNA(raw[[v]]) || any(!as.character(raw[[v]]) %in% coding[[v]]))
    .rc_error("adapter_covariates", paste("Source", v, "does not use the documented nonmissing coding."))
  raw <- raw[order(raw$id, raw$t.start, raw$t.stop), , drop = FALSE]
  divisor <- if (time_unit == "years") 365.25 else 1
  # Forecast covariates are fixed at each patient's first recorded interval.
  # Later Charlson measurements are retained separately, never moved backwards.
  first <- match(raw$id, raw$id)
  d <- data.frame(id = raw$id, unit = raw$id, start = raw$t.start / divisor,
    stop = raw$t.stop / divisor, event_type = as.character(raw$event),
    terminal = raw$death, female = as.integer(as.character(raw$sex[first]) == "Female"),
    chemo_treated = as.integer(as.character(raw$chemo[first]) == "Treated"),
    dukes_C = as.integer(as.character(raw$dukes[first]) == "C"),
    dukes_D = as.integer(as.character(raw$dukes[first]) == "D"),
    charlson_1_2 = as.integer(as.character(raw$charlson[first]) == "1-2"),
    charlson_3 = as.integer(as.character(raw$charlson[first]) == "3"),
    recorded_charlson = as.character(raw$charlson), stringsAsFactors = FALSE)
  checked <- as_trajectory(d, event_types = "1", time_unit = time_unit)
  if (length(attr(checked, "rc_schema")$gaps))
    .rc_error("adapter_gap", "The source contains an unrecorded gap; the exact adapter cannot fill it with zero events.")
  list(data = d, coarsening = NULL,
    formula = ~ female + chemo_treated + dukes_C + dukes_D + charlson_1_2 + charlson_3,
    event_types = "1", nonrecurrent = character(),
    metadata = list(dataset = "frailtypack::readmission", source = if (supplied) "supplied source data frame" else "installed canonical R package",
      source_version = if (requireNamespace("frailtypack", quietly = TRUE)) as.character(utils::packageVersion("frailtypack")) else NA_character_,
      patients = length(unique(d$id)), source_rows = nrow(raw), rows = nrow(d),
      recurrences = sum(d$event_type == "1"), terminal_events = sum(d$terminal),
      time_unit = time_unit, days_per_year = 365.25, observation = "exact",
      interpretation = "Rehospitalization episodes and terminal death; these are not distinct disease-type marks.",
      covariate_policy = "Baseline covariates use the first available interval. recorded_charlson preserves later source values without carrying them into the past.",
      data_status = "Public clinical data distributed by the source package; no patient rows are bundled in coherent4cast."))
}

.rc_legendre_unit <- function(nodes) {
  if (length(nodes) != 1L || !is.numeric(nodes) || !is.finite(nodes) ||
      nodes < 1 || nodes != as.integer(nodes) || nodes > 64)
    .rc_error("quadrature", "nodes must be an integer between 1 and 64.")
  if (nodes == 1L) return(list(x = 0.5, w = 1))
  J <- matrix(0, nodes, nodes)
  k <- seq_len(nodes - 1L); a <- k / sqrt(4 * k^2 - 1)
  J[cbind(k, k + 1L)] <- a; J[cbind(k + 1L, k)] <- a
  eig <- eigen(J, symmetric = TRUE); o <- order(eig$values)
  list(x = (eig$values[o] + 1) / 2, w = eig$vectors[1, o]^2)
}

#' Construct day-coarsened illness-death paths for the public colon trial
#' @export
rc_colon <- function(data = NULL, nodes = 4L, covariate = c("node4", "nodes"),
                     node_median = NULL, breaks = c(0, 1, 3, Inf)) {
  covariate <- match.arg(covariate); supplied <- !is.null(data)
  quadr <- .rc_legendre_unit(nodes)
  if (!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks) ||
      breaks[1L] != 0 || tail(breaks, 1L) != Inf || any(diff(breaks) <= 0))
    .rc_error("breaks", "Declare increasing baseline breaks in years from zero to Inf.")
  raw <- .rc_source_data("colon", "survival", data)
  .rc_adapter_columns(raw, c("id", "etype", "time", "status", "rx", "sex", "age", "node4", "nodes"))
  if (anyNA(raw$id) || anyNA(raw$etype) || any(!raw$etype %in% c(1, 2)))
    .rc_error("adapter_id", "Each patient must have a recurrence row (etype=1) and death row (etype=2).")
  if (!is.numeric(raw$time) || any(!is.finite(raw$time)) || any(raw$time <= 0) ||
      anyNA(raw$status) || any(!raw$status %in% c(0, 1)))
    .rc_error("adapter_time", "Follow-up days must be positive and finite and status must be binary.")
  if (any(raw$status == 1 & raw$time < 1))
    .rc_error("adapter_time", "A recorded event day below one cannot define the declared (day-1,day] interval.")
  if (anyNA(raw[, c("rx", "sex", "age", "node4")]) ||
      any(!is.finite(raw$age)) || any(!raw$sex %in% c(0, 1)) || any(!raw$node4 %in% c(0, 1)) ||
      any(!as.character(raw$rx) %in% c("Obs", "Lev", "Lev+5FU")))
    .rc_error("adapter_covariates", "Source treatment, sex, age and node4 must use the documented finite coding.")
  if (!is.numeric(raw$nodes) || any(!is.na(raw$nodes) & (!is.finite(raw$nodes) | raw$nodes < 0)))
    .rc_error("adapter_covariates", "nodes must contain nonnegative counts or documented missing values.")
  if (covariate == "nodes") {
    if (is.null(node_median)) node_median <- stats::median(raw$nodes, na.rm = TRUE)
    if (length(node_median) != 1L || !is.finite(node_median) || node_median < 0)
      .rc_error("node_median", "Supply a finite training-set median when node counts are all missing.")
  }
  cuts_days <- breaks[is.finite(breaks)] * 365.25
  integrate <- function(lo, hi) {
    if (lo >= hi) .rc_error("coarsening", "An event-time integration interval has zero or negative width.")
    cuts <- sort(unique(c(lo, hi, cuts_days[cuts_days > lo & cuts_days < hi])))
    x <- w <- numeric()
    for (j in seq_len(length(cuts) - 1L)) {
      width <- cuts[j + 1L] - cuts[j]
      x <- c(x, cuts[j] + width * quadr$x); w <- c(w, width * quadr$w)
    }
    list(x = x, w = w)
  }
  ids <- sort(unique(raw$id)); owners <- weights <- integer(); rows <- list()
  ties <- truncated <- numeric(); path_id <- 0L
  for (owner in seq_along(ids)) {
    g <- raw[raw$id == ids[owner], , drop = FALSE]
    if (nrow(g) != 2L || !setequal(g$etype, c(1, 2)))
      .rc_error("adapter_pair", "Exactly one etype=1 and one etype=2 row is required per patient.")
    r <- g[g$etype == 1, , drop = FALSE]; d <- g[g$etype == 2, , drop = FALSE]
    if (any(vapply(c("rx", "sex", "age", "node4", "nodes"), function(v) !identical(r[[v]], d[[v]]), logical(1))))
      .rc_error("adapter_covariates", "Paired endpoint rows disagree on baseline covariates.")
    end <- d$time; dead <- d$status; recurrence <- r$status
    if (!recurrence && r$time < end) { end <- r$time; dead <- 0; truncated <- c(truncated, ids[owner]) }
    if (recurrence && r$time > end)
      .rc_error("adapter_order", "Recorded recurrence lies after terminal follow-up.")
    x <- c(levamisole = as.integer(as.character(r$rx) == "Lev"),
      levamisole_5fu = as.integer(as.character(r$rx) == "Lev+5FU"),
      male = r$sex, age_per_10 = (r$age - 60) / 10)
    if (covariate == "node4") x <- c(x, node4 = r$node4) else
      x <- c(x, nodes_per_5 = if (is.na(r$nodes)) node_median / 5 else r$nodes / 5,
        nodes_missing = as.integer(is.na(r$nodes)))
    add_path <- function(rt, dt, weight) {
      if (!is.na(rt) && rt >= dt) .rc_error("adapter_order", "A quadrature path is not recurrence-before-death/follow-up.")
      path_id <<- path_id + 1L
      owners <<- c(owners, owner); weights <<- c(weights, weight)
      z <- if (!is.na(rt)) data.frame(id = path_id, unit = path_id, source_id = ids[owner],
        start = c(0, rt) / 365.25, stop = c(rt, dt) / 365.25,
        event_type = c("1", "0"), terminal = c(0, dead)) else
        data.frame(id = path_id, unit = path_id, source_id = ids[owner],
          start = 0, stop = dt / 365.25, event_type = "0", terminal = dead)
      for (v in names(x)) z[[v]] <- unname(x[v])
      rows[[path_id]] <<- z
    }
    if (recurrence && dead && r$time == d$time) {
      ties <- c(ties, ids[owner]); rs <- integrate(r$time - 1, r$time)
      for (j in seq_along(rs$x)) {
        ds <- integrate(rs$x[j], d$time)
        for (k in seq_along(ds$x)) add_path(rs$x[j], ds$x[k], rs$w[j] * ds$w[k] / 365.25^2)
      }
    } else {
      rs <- if (recurrence) integrate(r$time - 1, r$time) else list(x = NA_real_, w = 365.25)
      ds <- if (dead) integrate(end - 1, end) else list(x = end, w = 365.25)
      for (j in seq_along(rs$x)) for (k in seq_along(ds$x))
        add_path(rs$x[j], ds$x[k], rs$w[j] * ds$w[k] / 365.25^2)
    }
  }
  d <- do.call(rbind, rows); rownames(d) <- NULL
  checked <- as_trajectory(d, event_types = "1", time_unit = "years", observation = "coarsened")
  list(data = d, coarsening = list(owner = as.integer(owners), log_weights = log(weights)),
    formula = if (covariate == "node4") ~ levamisole + levamisole_5fu + male + age_per_10 + node4 else
      ~ levamisole + levamisole_5fu + male + age_per_10 + nodes_per_5 + nodes_missing,
    event_types = "1", nonrecurrent = "1", breaks = breaks,
    metadata = list(dataset = "survival::colon", source = if (supplied) "supplied source data frame" else "installed canonical R package",
      source_version = if (requireNamespace("survival", quietly = TRUE)) as.character(utils::packageVersion("survival")) else NA_character_,
      patients = length(ids), source_ids = ids, source_rows = nrow(raw), paths = path_id,
      recurrences = sum(raw$status[raw$etype == 1]), source_deaths = sum(raw$status[raw$etype == 2]),
      same_day_ids = ties, surveillance_truncated_ids = truncated,
      quadrature_nodes_per_segment = as.integer(nodes), breaks = breaks, time_unit = "years", days_per_year = 365.25,
      covariate = covariate, node_median = node_median,
      covariate_policy = "node4 is retained as supplied. Continuous nodes use the supplied training median, with a separate missingness indicator.",
      observation = "Event at reported day d is integrated over (d-1,d]; censoring endpoints are retained exactly.",
      legal_order = "Recorded recurrence precedes recorded death. Same-day pairs are integrated over the within-day triangle; post-death recurrence is excluded.",
      scope = "Illness-death submodel, not repeatable cancer recurrence. Gaussian quadrature approximates the coarsened likelihood; numerical weights are integration masses, not probabilities normalized within a patient.",
      forecast_scope = "These paths use full follow-up and must not be supplied as past information at an earlier landmark. Construct a landmark-specific observation fibre before filtering paths.",
      data_status = "Public clinical trial data distributed by survival; no patient rows are bundled in coherent4cast."))
}
