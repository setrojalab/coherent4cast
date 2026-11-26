.cs_env_fixture <- function(name, run, expected, timings = FALSE) {
  started <- if (timings) proc.time()[["elapsed"]] else NULL
  warnings <- character()
  result <- tryCatch(withCallingHandlers(run(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  elapsed <- if (timings) proc.time()[["elapsed"]] - started else NA_real_
  tolerance <- 1e-12
  if (inherits(result, "error")) {
    observed <- NULL; error <- NA_real_; status <- "failed"
    detail <- paste("Fixture raised an error:", conditionMessage(result))
  } else {
    observed <- result
    valid <- is.numeric(result) && identical(names(result), names(expected)) &&
      length(result) == length(expected) && all(is.finite(result))
    error <- if (valid) max(abs(result - expected)) else NA_real_
    status <- if (valid && error <= tolerance && !length(warnings)) "passed" else "failed"
    detail <- if (!valid) "Returned values are nonfinite or do not match the expected schema." else
      if (length(warnings)) paste("Fixture emitted warnings:", paste(warnings, collapse = " | ")) else
      if (status == "passed") "Agrees with the declared analytic reference." else "Analytic tolerance exceeded."
  }
  list(check = data.frame(check = name, category = "operational", status = status,
    detail = detail, maximum_absolute_error = error, tolerance = tolerance,
    elapsed_seconds = elapsed, stringsAsFactors = FALSE),
    reference = list(expected = expected, observed = observed,
      result_object_bytes = if (timings && !inherits(result, "error")) as.numeric(utils::object.size(result)) else NULL))
}

.cs_env_metadata <- function() {
  ns <- asNamespace("coherentSurv")
  path <- getNamespaceInfo(ns, "path")
  description <- tryCatch(read.dcf(file.path(path, "DESCRIPTION")), error = function(e) NULL)
  field <- function(name) if (!is.null(description) && name %in% colnames(description))
    unname(description[1L, name]) else NULL
  dll <- getLoadedDLLs()[["coherentSurv"]]
  routines <- if (is.null(dll)) NULL else getDLLRegisteredRoutines(dll)[[".Call"]]
  native_path <- if (is.null(dll)) NULL else dll[["path"]]
  native_hash <- if (is.null(native_path) || !file.exists(native_path)) NULL else
    tryCatch(digest::digest(file = native_path, algo = "sha256"), error = function(e) NULL)
  package_names <- c("coherentSurv", "Rcpp", "digest", "jsonlite")
  packages <- do.call(rbind, lapply(package_names, function(name) {
    # Namespace versions describe the code loaded in this process, not a different
    # installation that happens to appear first on a library search path.
    loaded <- name %in% loadedNamespaces()
    version <- if (loaded) as.character(getNamespaceVersion(name)) else
      tryCatch(as.character(utils::packageVersion(name)), error = function(e) NA_character_)
    data.frame(package = name, version = version, loaded = loaded,
      source = if (loaded) "loaded namespace" else "installed package metadata",
      status = if (is.na(version)) "unavailable" else "passed", stringsAsFactors = FALSE)
  }))
  # Built records R/platform/date; it does not identify the compiler used for the
  # installed binary. A compiler currently on PATH would not establish that fact.
  list(package = list(name = "coherentSurv", version = as.character(getNamespaceVersion(ns)),
      installed_path = path, built = field("Built"), needs_compilation = field("NeedsCompilation")),
    packages = packages,
    native = list(status = if (is.null(dll)) "unavailable" else "passed",
      path = native_path, file_sha256 = native_hash,
      fingerprint_scope = "Bytes currently on disk at the loaded DLL path; not a signed or historical build attestation.",
      dynamic_lookup = if (is.null(dll)) NULL else dll[["dynamicLookup"]],
      registered_calls = if (is.null(routines)) character() else names(routines)),
    build_compiler = list(status = "unavailable", compiler = NULL, flags = NULL,
      reason = "Actual compiler identity and flags are not recorded in the installed package metadata. No compiler is guessed or executed."))
}

.cs_env_text <- function(report) {
  c("coherentSurv environment diagnostic", paste("Recorded at UTC:", report$created_at_utc),
    paste("Overall diagnostic status:", report$status),
    paste("Operational fixtures:", report$operational_status),
    paste("Loaded package version:", report$installation$package$version),
    paste("R:", report$session$R_version), "", capture.output(print(report$checks, row.names = FALSE)),
    "", "Installation and native provenance:", capture.output(utils::str(report$installation)),
    "", "Session libraries:", capture.output(utils::str(report$session)),
    "", "Analytic references and observations:", capture.output(utils::str(report$fixtures)),
    "", "Interpretation:", report$interpretation)
}

#' Inspect the local runtime and execute bounded analytic diagnostics
#'
#' This offline diagnostic does not certify an installation or a study. It reads
#' installed metadata and executes four tiny deterministic native fixtures.
#' @param timings Include elapsed seconds and returned-object sizes for the tiny
#'   fixtures. These observations are not peak memory or comparative benchmarks.
#' @param format Return a structured list, JSON string or text lines.
#' @param file Optional output path. With format="list", write an RDS object;
#'   otherwise write JSON or text. No file is written when NULL.
#' @return A coherentSurv_environment report, or its requested serialization.
#' @export
cs_validate_env <- function(timings = FALSE, format = c("list", "json", "text"), file = NULL) {
  format <- match.arg(format)
  if (!is.logical(timings) || length(timings) != 1L || is.na(timings))
    .cs_error("environment_timings", "timings must be TRUE or FALSE.")
  if (!is.null(file) && (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)))
    .cs_error("environment_file", "file must be NULL or one nonempty output path.")
  # Rcpp's generated RNGScope may initialize a seed even for deterministic
  # kernels. Preserve both an existing state and the absence of a state.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_before <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", seed_before, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  metadata <- .cs_env_metadata()
  metadata_row <- function(name, status, detail) data.frame(check = name, category = "installation",
    status = status, detail = detail, maximum_absolute_error = NA_real_, tolerance = NA_real_,
    elapsed_seconds = NA_real_, stringsAsFactors = FALSE)
  required <- c("_coherentSurv_cpp_cs_loglik", "_coherentSurv_cpp_cs_coarsen",
                "_coherentSurv_cpp_cs_propagate", "_coherentSurv_cpp_cs_propagate_sparse")
  native_ok <- all(required %in% metadata$native$registered_calls) &&
    identical(metadata$native$dynamic_lookup, FALSE)
  installation <- rbind(
    metadata_row("installed_dependencies", if (all(metadata$packages$status == "passed")) "passed" else "unavailable",
      "Version metadata for coherentSurv, Rcpp, digest and jsonlite; no package is installed or downloaded."),
    metadata_row("registered_native_entry_points", if (native_ok) "passed" else "failed",
      "Required fixture entry points must be registered and dynamic symbol lookup disabled."),
    metadata_row("native_file_fingerprint", if (is.null(metadata$native$file_sha256)) "unavailable" else "passed",
      "SHA-256 of the file currently at the loaded DLL path; this does not authenticate the build."),
    metadata_row("package_build_compiler", metadata$build_compiler$status, metadata$build_compiler$reason))
  fixtures <- list(
    poisson_boundary = .cs_env_fixture("poisson_boundary", function() {
      z <- cs_loglik(log(.3), 1:2, c(1, 2), matrix(1, 2, 1), c(1, 0), phi_fixed = 0)
      c(subject_loglik = z$subject_loglik, scores = as.numeric(z$scores), information = as.numeric(z$information))
    }, c(subject_loglik = c(log(.3) - .3, -.6), scores = c(.7, -.6), information = .9), timings),
    dense_absorption = .cs_env_fixture("dense_absorption", function() {
      q <- matrix(c(-.3, .3, 0, 0), 2, 2, byrow = TRUE)
      z <- cs_propagate(c(1, 0), q, 2, reward = c(1, 0), tolerance = 1e-14)
      c(probability = z$probability, reward = z$reward)
    }, c(probability = c(exp(-.6), 1 - exp(-.6)), reward = -expm1(-.6) / .3), timings),
    sparse_absorption = .cs_env_fixture("sparse_absorption", function() {
      q <- cs_generator_csr(2, 1, 2, .3)
      z <- cs_propagate_sparse(c(1, 0), q, 2, reward = c(1, 0), tolerance = 1e-14,
                               max_terms = 100L, max_chunks = 10L)
      c(probability = z$probability, reward = z$reward)
    }, c(probability = c(exp(-.6), 1 - exp(-.6)), reward = -expm1(-.6) / .3), timings),
    louis_mixture = .cs_env_fixture("louis_mixture", function() {
      z <- cs_coarsen(log(c(.2, .8)), matrix(c(1, 3), 2, 1), array(c(4, 6), c(1, 1, 2)),
                      c(1, 1), log(c(.5, .5)))
      c(loglik = z$loglik, posterior = z$posterior, score = as.numeric(z$scores), information = as.numeric(z$information))
    }, c(loglik = log(.5), posterior = c(.2, .8), score = 2.6, information = 4.96), timings))
  operational <- do.call(rbind, lapply(fixtures, `[[`, "check"))
  checks <- rbind(installation, operational); rownames(checks) <- NULL
  status <- function(x) if (any(x == "failed")) "failed" else if (any(x == "unavailable")) "unavailable" else "passed"
  loaded <- sort(loadedNamespaces())
  software <- as.list(extSoftVersion())
  available <- function(value) list(status = if (is.null(value) || !length(value) ||
    anyNA(value) || !nzchar(value[1L])) "unavailable" else "passed", value = value)
  report <- structure(list(report_schema_version = "1.0", report_type = "coherentSurv_environment_diagnostic",
    created_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    status = status(checks$status), operational_status = status(operational$status),
    checks = checks, installation = metadata,
    session = list(R_version = R.version.string, platform = R.version$platform,
      locale = Sys.getlocale(), external_software = software,
      BLAS = available(software$BLAS),
      LAPACK = available(tryCatch(La_version(), error = function(e) NULL)),
      loaded_packages = data.frame(package = loaded,
        version = vapply(loaded, function(x) as.character(getNamespaceVersion(x)), character(1)), stringsAsFactors = FALSE)),
    fixtures = lapply(fixtures, `[[`, "reference"), timings_requested = timings,
    interpretation = c(
      "Read-only local diagnostics with fixed synthetic fixtures; no patient records, network access, installations or external compiler commands.",
      "A passed operational check establishes agreement only for these four fixtures at absolute tolerance 1e-12; it is not validation of every model or dataset.",
      "Unavailable provenance is retained explicitly and is not treated as either an observed failure or a successful qualification.",
      "Optional timings are single-run elapsed observations; zero can mean below timer resolution. Returned-object sizes are not peak resident memory. These are not machine speed benchmarks or performance qualification.",
      "This report is unsigned and mutable. It is not IQ/OQ/PQ certification, a signed qualification, 21 CFR Part 11 compliance, or a guarantee of clinical or regulatory fitness.",
      "File hashes describe current on-disk bytes, not the provenance of an in-memory library or proof that a source/compiler produced it. Environment paths may identify local installations.")),
    class = "coherentSurv_environment")
  value <- switch(format, list = report,
    json = as.character(jsonlite::toJSON(unclass(report), auto_unbox = TRUE, null = "null", na = "null",
                                       dataframe = "rows", digits = NA, pretty = TRUE)),
    text = .cs_env_text(report))
  if (!is.null(file)) {
    if (format == "list") saveRDS(value, file, version = 3L) else writeLines(enc2utf8(value), file, useBytes = TRUE)
  }
  value
}

#' @export
print.coherentSurv_environment <- function(x, ...) {
  cat(paste(.cs_env_text(x), collapse = "\n"), "\n")
  invisible(x)
}
