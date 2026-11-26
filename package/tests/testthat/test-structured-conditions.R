test_that("specific errors retain their original broad condition contract", {
  d <- data.frame(id = 1, start = 0, stop = 0, event_type = "0", terminal = 0)
  err <- tryCatch(cs_as_trajectory(d, time_unit = "day"), error = identity)
  expect_s3_class(err, "cs_error_interval")
  expect_s3_class(err, "coherentSurv_error")
  expect_s3_class(err, "error")
  expect_identical(err$code, "interval")
  expect_match(conditionMessage(err), "Require 0 <= start < stop", fixed = TRUE)
  handled <- tryCatch(cs_as_trajectory(d, time_unit = "day"),
    cs_error_interval = function(e) "specific handler")
  expect_identical(handled, "specific handler")
})

test_that("nonregular warnings can be caught without parsing prose", {
  warn <- get(".cs_warning", envir = environment(cs_fit))
  captured <- NULL
  withCallingHandlers(warn("nonregular_fit", "Diagnostic message"),
    cs_warning_nonregular_fit = function(w) {
      captured <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(captured, "coherentSurv_warning")
  expect_identical(captured$code, "nonregular_fit")
  expect_identical(conditionMessage(captured), "Diagnostic message")
})

test_that("audit exposes declared estimand and mapping without inventing them", {
  d <- cs_as_trajectory(data.frame(id = 1:2, start = 0, stop = 1,
    event_type = "0", terminal = c(1, 0)), time_unit = "day")
  declaration <- list(population = "Eligible adults", treatment = "Observed care",
    endpoint = "Death", intercurrent_events = "Death is absorbing",
    summary = "One-day terminal risk")
  fit <- structure(list(data = d, schema = attr(d, "cs_schema"),
    audit = list(estimand = declaration)), class = "coherentSurv")
  expect_identical(cs_audit_report(fit, include_session = FALSE)$model$estimand_declaration,
    declaration)
  expect_null(cs_audit_report(d, include_session = FALSE)$model$estimand_declaration)
  attr(d, "cs_mapping") <- list(source = "synthetic ADTTE", ascertainment = "Declared")
  expect_identical(cs_audit_report(d, include_session = FALSE)$model$mapping_provenance,
    attr(d, "cs_mapping"))
  attr(d, "cs_mapping")$origin <- data.frame(person_id = "PRIVATE_SUBJECT", origin = 0)
  safe <- cs_audit_report(d, include_session = FALSE, format = "json")
  expect_false(grepl("PRIVATE_SUBJECT", safe, fixed = TRUE))
  expect_false("origin" %in% names(cs_audit_report(d)$model$mapping_provenance))
})

test_that("causal estimation retains optional scientific declarations", {
  d <- cs_simulate(40, horizon = 1, rates = c(.9, .5, .4),
    dispersion = 0, censor_rate = 0, seed = 71)
  d$treatment <- as.integer(d$id %% 2); d$censored <- 0L
  declaration <- list(population = "Synthetic cohort", treatment = "Baseline 0 versus 1",
    endpoint = "Recurrent count and terminal risk", intercurrent_events = "Absorption at death",
    summary = "Population mean difference")
  fit <- cs_causal(~1, d, "treatment", 1, "censored", "years",
    propensity = .5, spec = list(complete_followup = TRUE, working = "pooled", grid = 16),
    estimand = declaration)
  expect_identical(fit$audit$estimand, declaration)
  expect_identical(summary(fit)$estimand_declaration, declaration)
  expect_identical(cs_audit_report(fit, include_session = FALSE)$model$estimand_declaration,
    declaration)
})
