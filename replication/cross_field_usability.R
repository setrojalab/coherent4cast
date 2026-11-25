# Run from the repository root with an installed coherent4cast 0.10.1.
# This is one synthetic cohort under six naming conventions, not six validation studies.
local_library <- Sys.getenv("COHERENT4CAST_LIBRARY", "")
if (nzchar(local_library)) {
  if (!dir.exists(local_library)) stop("COHERENT4CAST_LIBRARY is not an existing R library.")
  .libPaths(c(local_library, .libPaths()))
}
library(coherent4cast)
out <- Sys.getenv("COHERENT4CAST_CROSS_FIELD_OUTPUT", "replication/results/cross_field_usability")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(normalizePath(system.file(package = "coherent4cast"), winslash = "/"),
  file.path(out, "loaded_package_path.txt"))
criteria <- list()
check <- function(name, observed, requirement, passed, scope) {
  criteria[[length(criteria) + 1L]] <<- data.frame(check = name, observed = as.character(observed),
    requirement = requirement, passed = isTRUE(passed), scope = scope)
  write.csv(do.call(rbind, criteria), file.path(out, "criteria.csv"), row.names = FALSE)
}

# Independent subjects; independent simulated censoring; explicit time origin and unit.
d <- rc_simulate(n = 500, horizon = 4, rates = c(.5, .3, .12, .08),
  beta = c(.2, -.2, .1, .2), history = c(.05, .08, .12, .08),
  dispersion = .35, censor_rate = .03, history_cap = 1L,
  seed = 20261123L, time_unit = "years")
write.csv(d, file.path(out, "synthetic_histories.csv"), row.names = FALSE)
test_ids <- seq(5, 500, by = 5)
train <- d[!d$id %in% test_ids, ]
test <- d[d$id %in% test_ids, ]
control <- coherent4cast_control(maxit = 500L, tolerance = 1e-9, history_cap = 1L, quadrature = 30L)
warnings <- character()
fit <- withCallingHandlers(coherent4cast(~ x, train, event_types = c("1", "2", "3"),
  time_unit = "years", control = control, frailty_diagnostics = FALSE), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
saveRDS(fit, file.path(out, "training_fit.rds"))
write.csv(tidy(fit), file.path(out, "coefficients.csv"), row.names = FALSE)
write.csv(glance(fit), file.path(out, "fit_summary.csv"), row.names = FALSE)
check("Training optimizer regularity", fit$convergence$code, "converged and regular",
  fit$convergence$valid, "A failed check is retained; it is not repaired by changing the simulation seed.")

# The column role is unchanged. Event labels are examples, not causal assertions.
maps <- list(
  research = c("subject", "interval_start", "interval_end", "event_mark", "terminal_event", "baseline_predictor", "monitored"),
  hospital = c("patient", "care_start", "care_end", "morbidity_type", "death", "baseline_risk", "observed_care"),
  trial = c("participant", "followup_start", "followup_end", "adverse_event_type", "terminal_endpoint", "randomized_arm", "ascertained"),
  policy = c("service_user", "window_start", "window_end", "service_event", "permanent_exit", "baseline_need", "captured"),
  banking = c("account", "exposure_start", "exposure_end", "distress_event", "absorbing_default", "baseline_score", "reported"),
  engineering = c("asset", "operation_start", "operation_end", "fault_type", "decommissioning", "baseline_load", "instrumented"))
canonical <- c("id", "start", "stop", "event_type", "terminal", "x", "observed")
mapping <- do.call(rbind, lapply(names(maps), function(domain)
  data.frame(domain = domain, canonical_role = canonical, domain_column = maps[[domain]])))
write.csv(mapping, file.path(out, "schema_mapping.csv"), row.names = FALSE)
reference <- predict(fit, test, landmark = 1, horizon = c(.5, 1), type = "all")
write.csv(reference, file.path(out, "heldout_forecasts.csv"), row.names = FALSE)
parity <- list()
for (domain in names(maps)) {
  dd <- test[, canonical]; names(dd) <- maps[[domain]]
  named <- as_trajectory(dd, unit = maps[[domain]][1], start = maps[[domain]][2],
    stop = maps[[domain]][3], event_type = maps[[domain]][4], terminal = maps[[domain]][5],
    observed = maps[[domain]][7], event_types = c("1", "2", "3"), time_unit = "years")
  # Return to canonical columns before forecasting with the already fitted model.
  # A fit object retains the exact schema and formula used at training time.
  roundtrip <- as.data.frame(named)[, maps[[domain]]]; names(roundtrip) <- canonical
  p <- predict(fit, roundtrip, landmark = 1, horizon = c(.5, 1), type = "all")
  delta <- max(abs(p$estimate - reference$estimate))
  parity[[domain]] <- data.frame(domain = domain, rows = nrow(roundtrip),
    same_values = identical(unname(as.matrix(roundtrip)), unname(as.matrix(test[, canonical]))),
    max_forecast_difference = delta)
  check(paste(domain, "schema and prediction parity"), delta, "absolute difference <= 1e-12",
    delta <= 1e-12, "Naming/schema interoperability on the same synthetic cohort; not field validation.")
}
write.csv(do.call(rbind, parity), file.path(out, "domain_parity.csv"), row.names = FALSE)

# One direct fit through an entirely renamed public schema verifies dispatch as well.
renamed <- train[, canonical]; names(renamed) <- maps$engineering
fit_renamed <- withCallingHandlers(coherent4cast(~ baseline_load, renamed,
  unit = "asset", start = "operation_start", stop = "operation_end",
  event_type = "fault_type", terminal = "decommissioning", observed = "instrumented",
  event_types = c("1", "2", "3"), time_unit = "years", control = control,
  frailty_diagnostics = FALSE), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
check("Direct renamed-model coefficient parity", max(abs(coef(fit) - coef(fit_renamed))),
  "absolute difference <= 1e-8", max(abs(coef(fit) - coef(fit_renamed))) <= 1e-8,
  "Same numerical model and subjects; different input column names.")

r <- reference[reference$target == "risk", ]
s <- reference[reference$target == "survival", ]
check("Terminal probability and survival identity", max(abs(r$estimate + s$estimate - 1)),
  "absolute difference <= 1e-10", max(abs(r$estimate + s$estimate - 1)) <= 1e-10,
  "Common-law numerical identity; not empirical calibration.")
check("Forecast support", paste(range(reference$estimate), collapse = ", "),
  "risk/survival in [0,1], RMST in [0,horizon], burden >= 0",
  all(r$estimate >= 0 & r$estimate <= 1) && all(s$estimate >= 0 & s$estimate <= 1) &&
    all(reference$estimate[reference$target == "rmst"] >= 0) &&
    all(reference$estimate[reference$target == "rmst"] <= reference$horizon[reference$target == "rmst"] + 1e-10) &&
    all(reference$estimate[grepl("burden", reference$target)] >= 0),
  "No empirical superiority target is imposed.")
endpoint <- do.call(rbind, lapply(split(test, test$id), function(g)
  data.frame(id = g$id[1], time = tail(g$stop, 1), status = tail(g$terminal, 1))))
at <- reference[reference$target == "risk" & reference$horizon == 1, ]
ep <- endpoint[match(at$unit, as.character(endpoint$id)), ]
val <- validate(at$estimate, time = ep$time - 1, status = ep$status, horizon = 1, threshold = .2)
write.csv(val$metrics, file.path(out, "validation_metrics.csv"), row.names = FALSE)
write.csv(val$confusion, file.path(out, "weighted_confusion.csv"))
saveRDS(val, file.path(out, "validation_full.rds"))
write.csv(data.frame(predicted = at$estimate, time = ep$time - 1, status = ep$status),
  file.path(out, "validation_input.csv"), row.names = FALSE)
if (file.exists("replication/research_metrics.R")) {
  source("replication/research_metrics.R")
  report <- rc_research_report(at$estimate, time = ep$time - 1, status = ep$status,
    horizon = 1, thresholds = c(.2, .05, .1, .3, .5), bootstrap_reps = 200L, seed = 20261123L,
    false_positive_cost = 1, false_negative_cost = 5, loss_unit = "illustrative loss units")
  saveRDS(report, file.path(out, "companion_report.rds"))
  write.csv(report$metrics, file.path(out, "companion_metrics.csv"), row.names = FALSE)
  write.csv(report$thresholds, file.path(out, "companion_thresholds.csv"), row.names = FALSE)
  write.csv(report$precision_recall, file.path(out, "companion_precision_recall.csv"), row.names = FALSE)
  if (is.data.frame(report$bootstrap)) write.csv(report$bootstrap, file.path(out, "companion_bootstrap.csv"), row.names = FALSE)
  if (!is.null(report$bootstrap_draws)) write.csv(data.frame(replicate = seq_len(nrow(report$bootstrap_draws)),
    report$bootstrap_draws, check.names = FALSE), file.path(out, "companion_bootstrap_draws.csv"), row.names = FALSE)
  if (!is.null(report$bootstrap_errors)) write.csv(data.frame(replicate = seq_along(report$bootstrap_errors),
    error = report$bootstrap_errors), file.path(out, "companion_bootstrap_errors.csv"), row.names = FALSE)
  valuation <- rc_value_burden(reference, unit_cost = c(burden_1 = 1, burden_2 = 4, burden_3 = 10),
    cost_unit = "illustrative loss units")
  saveRDS(valuation, file.path(out, "companion_valuation.rds"))
  write.csv(valuation$total, file.path(out, "companion_valuation.csv"), row.names = FALSE)
  write.csv(valuation$components, file.path(out, "companion_valuation_components.csv"), row.names = FALSE)
}
check("Held-out alignment", length(at$estimate), "one aligned prediction per eligible held-out subject",
  !anyNA(ep$id) && !anyDuplicated(at$unit) && all(ep$time > 1), "Training and test subjects are disjoint.")
check("Metric reporting contract", nrow(val$metrics), "10 metric rows, with explicit defined flags",
  nrow(val$metrics) == 10L && identical(val$metrics$defined, is.finite(val$metrics$estimate)),
  "Undefined predictive values are retained as NA; no arbitrary success score.")
undefined <- validate(c(.1, .1, .1, .1), time = rep(2, 4), status = rep(0, 4), horizon = 1)
write.csv(undefined$metrics, file.path(out, "zero_case_metrics.csv"), row.names = FALSE)
check("Zero-case uncertainty is explicit", paste(undefined$metrics$metric[!undefined$metrics$defined], collapse = ", "),
  "AUROC, sensitivity and PPV undefined", all(!undefined$metrics$defined[match(c("AUROC", "sensitivity", "PPV"), undefined$metrics$metric)]),
  "A zero event count does not establish perfect discrimination or safety.")
unsupported <- tryCatch({validate(c(.1, .2), time = c(.1, .2), status = c(0, 0), horizon = 1); NULL}, error = identity)
check("No observed horizon outcome is rejected", if (is.null(unsupported)) "accepted" else unsupported$code,
  "structured validation_support error", inherits(unsupported, "coherent4cast_error") && identical(unsupported$code, "validation_support"),
  "All-censored-before-horizon data cannot provide a zero-error validation score.")

# Independent small finite-path MAR example; this is not an EIF for the frailty fit.
set.seed(20261123L)
features <- rbind(c(0, 0, 0, 0), c(1, 0, 0, 0), c(1, 1, 0, 0), c(1, 1, 1, 0), c(0, 0, 0, 1), c(1, 1, 1, 1))
colnames(features) <- c("burden_1", "burden_2", "burden_3", "risk")
z <- rbinom(400L, 1L, .5)
path <- vapply(z, function(x) sample.int(6, 1, prob = if (x == 0) c(.4,.3,.15,.08,.04,.03) else c(.25,.3,.2,.12,.06,.07)), integer(1))
recorded <- rbinom(length(z), 1, ifelse(z == 0, .8, .65))
target_data <- data.frame(id = seq_along(z), x = z, path = ifelse(recorded == 1, path, NA), recorded = recorded)
target <- coherent4cast(~ x, target_data, time_unit = "years", targeting = "finite_path",
  targeting_spec = list(path = "path", recorded = "recorded", strata = "x", features = features, horizon = 1))
write.csv(summary(target)$estimates, file.path(out, "finite_path_EIF.csv"), row.names = FALSE)
outer <- predict(target, interval = "projection", type = "all")
write.csv(outer, file.path(out, "finite_path_outer_intervals.csv"), row.names = FALSE)
check("Finite-path common-law agreement", max(abs(target$estimate - as.vector(crossprod(target$mixture, features)))),
  "absolute difference <= 1e-10", max(abs(target$estimate - as.vector(crossprod(target$mixture, features)))) <= 1e-10,
  "MAR finite-path estimands, separately specified from the continuous-time frailty likelihood.")
check("EIF score equation", max(abs(target$EIF_mean)), "convergence flag true",
  target$convergence$valid, "A solved score is a numerical check; it does not establish finite-sample efficiency.")
counts <- tabulate(path[recorded == 1], 6)
b <- rc_bounds(counts / sum(counts), numerator = features[, "risk"], kappa = 2,
  recording_fraction = .5 * .8 + .5 * .65, counts = counts)
write.csv(data.frame(plugin_lower = b$lower, plugin_upper = b$upper,
  outer_lower = b$outer$lower, outer_upper = b$outer$upper), file.path(out, "recording_bounds.csv"), row.names = FALSE)
saveRDS(b, file.path(out, "recording_bounds_full.rds"))
check("Bound nesting", paste(b$outer$lower, b$lower, b$upper, b$outer$upper, sep = ", "),
  "outer lower <= plugin lower <= plugin upper <= outer upper",
  b$outer$lower <= b$lower + 1e-10 && b$lower <= b$upper && b$upper <= b$outer$upper + 1e-10,
  "Declared finite support, known simulated marginal recording fraction and kappa; not unrestricted MNAR.")

bad <- test[test$id == test$id[1], ]; bad$start[1] <- .01
error <- tryCatch({predict(fit, bad, landmark = 1, horizon = 1); NULL}, error = identity)
check("Delayed entry is rejected", if (is.null(error)) "accepted" else error$code,
  "structured left_truncation error", inherits(error, "coherent4cast_error") && identical(error$code, "left_truncation"),
  "Unobserved pre-entry histories require an observation model.")
error2 <- tryCatch({coherent4cast(~ x, train, time_unit = "years", intervention = list(treatment = 1)); NULL}, error = identity)
check("Unsupported causal request is rejected", if (is.null(error2)) "accepted" else error2$code,
  "structured causal_scope error", inherits(error2, "coherent4cast_error") && identical(error2$code, "causal_scope"),
  "Prognostic scores do not identify intervention effects.")
writeLines(unique(warnings), file.path(out, "warnings.txt"))
writeLines(c("One synthetic cohort, six column-name mappings, one independent renamed fit.",
  "No six-domain external-validation claim; no threshold was tuned using this held-out cohort.",
  "Native validate() returns point estimates; companion uncertainty must be reported separately.",
  "Known synthetic recording fraction is supplied only in the finite-support demonstration.",
  "The finite-path target example is a different declared observation experiment from the main frailty fit."),
  file.path(out, "INTERPRETATION.txt"))
print(do.call(rbind, criteria), row.names = FALSE)
print(val$metrics, row.names = FALSE)
cat("Outputs:", normalizePath(out, winslash = "/"), "\n")
