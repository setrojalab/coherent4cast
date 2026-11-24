.rc_audit_time <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")

.rc_audit_call <- function(call) {
  if (is.null(call)) return(NULL)
  if (!is.call(call)) return("Call was not captured as an R call; no expression is reconstructed.")
  values <- as.list(call)
  protected <- c("data", "newdata", "paths", "weights", "subset", "coarsening", "targeting_spec",
    "newdata_future", "covariate_history", "initial_laws", "recording", "labels",
    "spec", "propensity", "folds")
  named <- names(values)
  for (i in seq_along(values)[-1L]) {
    # match.call() names supplied arguments. Unnamed arbitrary expressions and
    # data-bearing arguments are not safe to reproduce as audit text.
    if (is.null(named) || !nzchar(named[i]) || named[i] %in% protected)
      values[[i]] <- quote(.audit_data_argument_redacted)
  }
  paste(deparse(as.call(values), width.cutoff = 500L), collapse = " ")
}

.rc_audit_session <- function() {
  list(R_version = R.version.string, platform = R.version$platform,
    recurcast_version = tryCatch(as.character(utils::packageVersion("recurcast")), error = function(e) NA_character_),
    session_info = capture.output(utils::sessionInfo()))
}

.rc_audit_data_fingerprint <- function(data) {
  raw <- serialize(data, NULL, version = 3L)
  list(algorithm = "SHA-256", value = digest::digest(raw, algo = "sha256", serialize = FALSE),
    representation = "R serialize(data, NULL, version = 3); SHA-256 of the resulting raw bytes",
    serialization_version = 3L, rows = if (is.data.frame(data)) nrow(data) else NULL,
    scope = "Integrity of this serialized object, including its attributes; not semantic equality across R versions, schema reorderings or platforms.")
}

# Called by the fitting interface when context is available; never reconstruct
# historical execution context at report-export time.
.rc_capture_audit_context <- function(data, call = NULL) {
  list(captured_at_utc = .rc_audit_time(), fit_session = .rc_audit_session(),
    data_fingerprint = .rc_audit_data_fingerprint(data), captured_call = .rc_audit_call(call))
}

.rc_audit_sources <- function(source_files, exported_at) {
  empty <- data.frame(path = character(), algorithm = character(), sha256 = character(),
    bytes = numeric(), modified_at_utc = character(), hashed_at_utc = character(), stringsAsFactors = FALSE)
  if (is.null(source_files) || !length(source_files)) return(empty)
  if (!is.character(source_files) || anyNA(source_files) || any(!nzchar(source_files)))
    .rc_error("audit_source", "source_files must contain explicitly supplied, nonmissing file paths.")
  if (any(!file.exists(source_files)) || any(dir.exists(source_files)))
    .rc_error("audit_source", "Every supplied source path must be an existing regular file; directories are not scanned.")
  paths <- unique(normalizePath(source_files, winslash = "/", mustWork = TRUE))
  info <- file.info(paths)
  hashes <- vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))
  data.frame(path = paths, algorithm = "SHA-256", sha256 = unname(hashes), bytes = info$size,
    modified_at_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    hashed_at_utc = exported_at, stringsAsFactors = FALSE, row.names = NULL)
}

.rc_audit_numerical_summary <- function(predictions) {
  if (is.null(predictions)) return(NULL)
  if (!is.data.frame(predictions)) .rc_error("audit_prediction", "predictions must be a prediction data frame.")
  fields <- intersect(c("series_error_bound", "probability_roundoff_bound", "frailty_compression_bound",
    "covariate_mcse", "interpolation_error"), names(predictions))
  if (!length(fields)) return(data.frame(metric = character(), finite = integer(), nonfinite = integer(),
    minimum = numeric(), maximum = numeric(), stringsAsFactors = FALSE))
  do.call(rbind, lapply(fields, function(name) {
    x <- predictions[[name]]
    if (!is.numeric(x)) .rc_error("audit_prediction", "Numerical error columns must be numeric.")
    finite <- is.finite(x)
    data.frame(metric = name, finite = sum(finite), nonfinite = sum(!finite),
      minimum = if (any(finite)) min(x[finite]) else NA_real_,
      maximum = if (any(finite)) max(x[finite]) else NA_real_, stringsAsFactors = FALSE)
  }))
}

.rc_audit_text <- function(report) {
  show <- function(x) if (is.null(x)) "Not captured / not supplied." else capture.output(print(x))
  c("RecurCast audit report", paste("Report schema:", report$report_schema_version),
    paste("Exported at UTC:", report$export$created_at_utc),
    paste("Fit context captured at UTC:", if (is.null(report$fit_context$captured_at_utc)) "not captured" else report$fit_context$captured_at_utc),
    "", "Captured call (data-bearing expressions redacted):", show(report$model$call),
    "", "Model and observation specification:", show(report$model[setdiff(names(report$model), "call")]),
    "", "Fit-time data fingerprint:", show(report$fit_context$data_fingerprint),
    "", "Legacy fingerprint (algorithm retained):", show(report$fit_context$legacy_fingerprint),
    "", "Current data fingerprint at export:", show(report$export$data_fingerprint),
    "", "Convergence diagnostics:", show(report$convergence),
    "", "Fitting numerical diagnostics:", show(report$numerical$fitting),
    "", "Prediction numerical error summaries:", show(report$numerical$prediction),
    "", "Explicitly supplied source files, hashed at export:", show(report$source_files),
    "", "Session captured with fit:", show(report$fit_context$session),
    "", "Session at export:", show(report$export$session),
    "", "Interpretation:", report$interpretation)
}

#' Export provenance and numerical diagnostics without exporting subject records
#' @export
rc_audit_report <- function(object, predictions = NULL, source_files = NULL,
                            format = c("list", "json", "text"), file = NULL,
                            include_session = TRUE) {
  format <- match.arg(format)
  if (!inherits(object, c("recurcast", "recurcast_competing", "surv_trajectory", "recurcast_prediction")))
    .rc_error("audit_object", "Use a fitted RecurCast object, validated trajectory or RecurCast prediction table.")
  if (!is.logical(include_session) || length(include_session) != 1L || is.na(include_session))
    .rc_error("audit_session", "include_session must be TRUE or FALSE.")
  if (!is.null(file) && (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)))
    .rc_error("audit_file", "file must be NULL or one output path.")
  is_data <- inherits(object, "surv_trajectory")
  is_prediction <- inherits(object, "recurcast_prediction")
  data <- if (is_data) object else if (!is_prediction && is.data.frame(object$data)) object$data else NULL
  audit <- if (is_data || is_prediction || is.null(object$audit)) list() else object$audit
  schema <- if (is_data) attr(object, "rc_schema") else if (!is_prediction) object$schema else NULL
  if (is_prediction && is.null(predictions)) predictions <- object
  exported_at <- .rc_audit_time()
  legacy <- if (is.null(audit$data_hash)) NULL else list(value = audit$data_hash,
    algorithm = if (is.null(audit$data_hash_algorithm)) "MD5 (legacy data_hash field; original encoding not captured)" else audit$data_hash_algorithm,
    scope = "Original fit audit value retained verbatim; it has not been relabelled as SHA-256.")
  schema_summary <- if (is.null(schema)) NULL else schema[intersect(names(schema),
    c("unit", "start", "stop", "event_type", "terminal", "observed", "event_types", "time_unit",
      "observation", "original_rows", "dropped_rows", "terminal_cause", "terminal_causes_observed", "regime"))]
  context <- list(captured_at_utc = audit$captured_at_utc,
    session = if (include_session) audit$fit_session else NULL,
    data_fingerprint = audit$data_fingerprint, legacy_fingerprint = legacy,
    status = if (is.null(audit$captured_at_utc)) "Historical fit timestamp/session were not captured; export context is not substituted." else "Context captured by the fitting workflow.")
  call <- if (!is.null(audit$captured_call)) audit$captured_call else .rc_audit_call(audit$call)
  formula <- if (!is_data && !is_prediction) object$formula else NULL
  if (is.null(formula) && !is_data && !is_prediction && !is.null(object$wrapper)) formula <- object$wrapper$formula
  model <- list(class = class(object), call = call,
    wrapper_call = .rc_audit_call(audit$wrapper_call),
    formula = if (inherits(formula, "formula")) paste(deparse(formula), collapse = " ") else NULL,
    schema = schema_summary,
    regime_details = if (is.null(schema)) NULL else schema$regime_details,
    observation_process = if (is.null(data)) NULL else attr(data, "observation_process"),
    observation_process_details = if (is.null(data)) NULL else attr(data, "observation_process_details"),
    independent_units = if (!is_data && !is_prediction) object$n else if (is_data) length(unique(data$.rc_unit)) else NULL,
    complete_paths = if (!is_data && !is_prediction) object$npaths else NULL,
    targeting = if (inherits(object, "recurcast_competing")) "none" else if (!is_data && !is_prediction) object$targeting else NULL,
    frailty = if (inherits(object, "recurcast_competing")) "none" else if (!is_data && !is_prediction) object$frailty else NULL,
    absorbing_causes = if (inherits(object, "recurcast_competing")) object$causes else NULL,
    seed = audit$seed, assumptions = audit$assumptions)
  numerical <- list(fitting = audit[intersect(names(audit), c("quadrature_refinement", "quadrature_stable", "numerical_error", "minimum_censor_survival"))],
    prediction = .rc_audit_numerical_summary(predictions),
    prediction_assumptions = if (is.null(predictions)) NULL else attr(predictions, "assumptions"),
    scope = "Numerical bounds and simulation errors are not sampling confidence intervals. Quadrature, time-grid, propagation, covariate Monte Carlo and parameter uncertainty remain distinct.")
  report <- list(report_schema_version = "1.0", report_type = "recurcast_audit_report",
    fit_context = context, model = model,
    convergence = if (!is_data && !is_prediction) object$convergence else NULL,
    numerical = numerical, source_files = .rc_audit_sources(source_files, exported_at),
    export = list(created_at_utc = exported_at, session = if (include_session) .rc_audit_session() else NULL,
      data_fingerprint = if (is.null(data)) NULL else .rc_audit_data_fingerprint(data),
      fingerprint_scope = "Current object data at export, not a reconstruction of historical fit data or source provenance."),
    interpretation = c("Subject-record tables and identifier columns are not exported automatically. Model expressions, category labels and explicitly supplied paths should still be reviewed before sharing. Hashes are integrity fingerprints, not anonymization guarantees.",
      "Only explicitly supplied source_files are read and hashed. These hashes describe files at export time, not proof that those files produced the fit.",
      "No success flag certifies model assumptions, causal identification, finite-sample coverage, publication acceptance or regulatory validity.",
      "A missing fit timestamp/session remains missing. Export metadata is stored separately.",
      "JSON null represents unavailable or nonfinite numeric values; consult convergence flags and nonfinite counts."))
  value <- switch(format, list = report,
    json = as.character(jsonlite::toJSON(report, auto_unbox = TRUE, null = "null", na = "null",
      dataframe = "rows", digits = NA, pretty = TRUE)), text = .rc_audit_text(report))
  if (!is.null(file)) {
    if (format == "list") saveRDS(value, file, version = 3L) else writeLines(enc2utf8(value), file, useBytes = TRUE)
  }
  value
}
