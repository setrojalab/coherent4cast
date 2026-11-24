# Baseline treatment extension: Supplement S12, equations S12.1-S12.6.
# All subjects standardize both arms; only observed-arm records enter each
# arm's martingale correction. Treatment is never inferred from event marks.
rc_causal <- function(formula, data, treatment, horizon, censor, time_unit,
                      unit = "id", start = "start", stop = "stop",
                      event_type = "event_type", terminal = "terminal",
                      observed = NULL, event_types = NULL, folds = NULL,
                      propensity = NULL, spec = list(), seed = 1L) {
  call <- match.call()
  if (!inherits(formula, "formula") || length(formula) != 2L)
    .rc_error("causal_formula", "Use a one-sided formula for baseline confounder strata.")
  if (!is.list(spec) || length(setdiff(names(spec), c("strata", "nfolds", "grid",
      "history_cap", "working", "nonrecurrent", "complete_followup", "initial_rates",
      "censor_rates", "min_propensity", "min_censor_survival", "maxit", "tolerance"))))
    .rc_error("causal_spec", "Unknown causal specification field.")
  .rc_scalar(horizon, "horizon", 0, Inf, TRUE)
  value <- function(name, default) .rc_target_spec_value(spec, name, default)
  for (name in list(treatment, censor))
    if (!is.character(name) || length(name) != 1L || is.na(name) || !name %in% names(data))
      .rc_error("causal_schema", "Declare existing treatment and censor column names.")
  if (identical(treatment, censor)) .rc_error("causal_schema", "Treatment and censor must be different columns.")
  d <- as_trajectory(data, unit, start, stop, event_type, terminal, observed,
                     event_types, time_unit, "exact")
  schema <- attr(d, "rc_schema")
  if (length(schema$gaps) || length(schema$delayed_entry) || any(!d$.rc_observed))
    .rc_error("causal_observation", "Require contiguous exact histories from zero, with explicit right censoring.")
  if (anyNA(d[[treatment]]) || !all(d[[treatment]] %in% c(0, 1)) ||
      anyNA(d[[censor]]) || !all(d[[censor]] %in% c(0, 1)))
    .rc_error("causal_schema", "Treatment and right-censor indicators must be binary 0/1.")
  strata <- value("strata", all.vars(formula))
  if (!is.character(strata) || anyNA(strata) || anyDuplicated(strata) ||
      !all(strata %in% names(d)) || !all(all.vars(formula) %in% strata) ||
      any(c(treatment, censor, start, stop, event_type, terminal, unit) %in% strata))
    .rc_error("causal_strata", "Declare baseline confounder columns only; include every formula variable.")
  subjects <- split(seq_len(nrow(d)), d$.rc_unit)
  n <- length(subjects)
  if (n < 8L) .rc_error("causal_subjects", "At least eight independent subjects are required.")
  first <- vapply(subjects, `[`, integer(1), 1)
  for (ix in subjects) {
    if (length(unique(d[[treatment]][ix])) != 1L)
      .rc_error("causal_treatment", "Treatment must be assigned at baseline and remain fixed.")
    if (length(strata) && any(vapply(d[ix, strata, drop = FALSE],
        function(z) anyNA(z) || length(unique(z)) != 1L, logical(1))))
      .rc_error("causal_strata", "Confounder strata must be observed and fixed at baseline.")
    if (any(d[[censor]][head(ix, -1)] != 0) ||
        any(d[[censor]][ix] == 1 & (d$.rc_mark[ix] > 0 | d$.rc_terminal[ix] == 1)))
      .rc_error("causal_censor", "Censoring must end the history without an event at that time.")
    last <- tail(ix, 1)
    if (d$.rc_stop[last] < horizon && d$.rc_terminal[last] == 0 && d[[censor]][last] != 1)
      .rc_error("causal_censor", "Follow-up ending before the horizon needs an explicit censor or terminal indicator.")
  }
  w <- as.integer(d[[treatment]][first])
  if (!setequal(unique(w), 0:1)) .rc_error("causal_arms", "Both baseline treatment arms must be represented.")
  keys <- if (length(strata)) do.call(paste, c(d[first, strata, drop = FALSE], sep = "\r")) else rep("all", n)
  levels <- sort(unique(keys)); x <- match(keys, levels); nx <- length(levels)
  if (nx > 30L) .rc_error("causal_strata", "This implementation supports at most 30 prespecified discrete baseline strata.")
  nf <- value("nfolds", 5L); .rc_scalar(nf, "number of folds", 2, 20)
  if (nf != floor(nf)) .rc_error("causal_folds", "Number of folds must be integer.")
  if (is.character(folds) && length(folds) == 1L && folds %in% names(d)) {
    if (any(vapply(subjects, function(ix) length(unique(d[[folds]][ix])) != 1L, logical(1))))
      .rc_error("causal_folds", "One fold per subject is required.")
    folds <- d[[folds]][first]
  }
  if (is.null(folds)) folds <- .rc_target_rng(seed, {
    z <- integer(n)
    for (a in seq_len(nx)) for (arm in 0:1) {
      ix <- which(x == a & w == arm)
      z[ix] <- sample(rep(seq_len(nf), length.out = length(ix)))
    }
    z
  })
  if (length(folds) != n || anyNA(folds)) .rc_error("causal_folds", "Supply one fold label per subject in first-appearance order.")
  folds <- match(folds, sort(unique(folds))); nf <- max(folds)
  if (nf < 2L || nf > 20L) .rc_error("causal_folds", "Use between two and twenty folds.")
  group <- (folds - 1L) * nx + x; ng <- nf * nx
  L <- value("grid", 64L); cap <- value("history_cap", 2L)
  .rc_scalar(L, "grid", 4, 512); .rc_scalar(cap, "history cap", 0, 6)
  K <- length(schema$event_types); M <- K + 1L; J <- K + 2L
  if (L != floor(L) || cap != floor(cap) || K > 6L || (cap + 1)^K > 128)
    .rc_error("causal_state", "Use integer grid/cap, at most six recurrent marks and 128 sufficient living states.")
  states <- if (K) as.matrix(expand.grid(rep(list(0:cap), K)))[, rev(seq_len(K)), drop = FALSE] else matrix(numeric(), 1, 0)
  S <- nrow(states); powers <- if (K) (cap + 1)^rev(seq.int(0, K - 1L)) else numeric()
  nxt <- matrix(as.integer(S), S, M)
  if (K) for (q in seq_len(S)) for (k in seq_len(K)) {
    s <- states[q, ]; s[k] <- min(cap, s[k] + 1L); nxt[q, k] <- sum(s * powers)
  }
  storage.mode(nxt) <- "integer"
  history <- if (K) build_history(d, cap = cap)$counts else matrix(numeric(), nrow(d), 0)
  state <- if (K) as.integer(history %*% powers) + 1L else rep(1L, nrow(d))
  keep <- d$.rc_start < horizon
  dd <- d[keep, , drop = FALSE]; state <- state[keep]
  after <- dd$.rc_stop > horizon
  dd$.rc_stop <- pmin(dd$.rc_stop, horizon)
  dd$.rc_mark[after] <- dd$.rc_terminal[after] <- dd[[censor]][after] <- 0L
  u <- dd$.rc_unit; exposure <- dd$.rc_stop - dd$.rc_start
  event <- ifelse(dd$.rc_terminal == 1, K, ifelse(dd$.rc_mark > 0,
                  dd$.rc_mark - 1L, ifelse(dd[[censor]] == 1, M, -1L)))
  records <- cbind(u - 1L, group[u] - 1L, state - 1L, dd$.rc_start, dd$.rc_stop, event)
  if (!is.logical(value("complete_followup", FALSE)) || length(value("complete_followup", FALSE)) != 1L || is.na(value("complete_followup", FALSE)))
    .rc_error("causal_followup", "complete_followup must be TRUE or FALSE.")
  complete <- isTRUE(value("complete_followup", FALSE))
  if (complete && any(dd[[censor]] == 1 & dd$.rc_stop < horizon))
    .rc_error("causal_followup", "The complete-followup declaration conflicts with early right censoring.")
  nonrecurrent <- match(value("nonrecurrent", character()), schema$event_types)
  if (anyNA(nonrecurrent)) .rc_error("causal_eligibility", "A nonrecurrent mark is absent from event_types.")
  for (k in nonrecurrent)
    if (any(dd$.rc_mark == k & history[keep, k] > 0)) .rc_error("causal_eligibility", "A declared nonrecurrent mark recurs.")
  working <- value("working", "loglinear")
  if (length(working) != 1L || !working %in% c("loglinear", "pooled")) .rc_error("causal_working", "working must be loglinear or pooled.")
  min_e <- value("min_propensity", .01); min_G <- value("min_censor_survival", .001)
  .rc_scalar(min_e, "minimum propensity", 0, .5, TRUE)
  .rc_scalar(min_G, "minimum censor survival", 0, 1, TRUE)
  e <- numeric(ng); training <- vector("list", nf)
  for (f in seq_len(nf)) {
    training[[f]] <- which(folds != f)
    for (a in seq_len(nx)) {
      tr <- which(folds != f & x == a)
      if (!length(tr) || !setequal(unique(w[tr]), 0:1))
        .rc_error("causal_overlap", "Every training fold/stratum needs both treatment arms; no clipping or arm extrapolation is applied.")
      e[(f - 1L) * nx + a] <- (sum(w[tr]) + .5) / (length(tr) + 1)
    }
  }
  if (!is.null(propensity)) {
    if (!is.numeric(propensity) || any(!is.finite(propensity))) .rc_error("causal_propensity", "Supply finite numeric treatment-1 probabilities.")
    if (length(propensity) == 1L) e[] <- propensity else if (is.matrix(propensity) && identical(dim(propensity), c(nf, nx)))
      e <- as.vector(t(propensity)) else .rc_error("causal_propensity", "propensity must be one known probability or a fold-by-stratum matrix.")
  }
  if (any(e < min_e | e > 1 - min_e)) .rc_error("causal_positivity", "A treatment probability is outside the declared overlap range; no clipping is applied.")
  maxit <- value("maxit", 100L); tolerance <- value("tolerance", 1e-7 / sqrt(n))
  .rc_scalar(maxit, "maximum iterations", 0, 10000); .rc_scalar(tolerance, "tolerance", 0, Inf, TRUE)
  if (maxit != floor(maxit)) .rc_error("causal_iterations", "maxit must be integer.")
  for (field in c("initial_rates", "censor_rates")) if (!is.null(spec[[field]]) &&
      (!is.list(spec[[field]]) || !setequal(names(spec[[field]]), c("0", "1"))))
    .rc_error("causal_nuisance", "Supplied initial_rates/censor_rates must be named lists with entries '0' and '1'.")
  targets <- c(if (K) paste0("burden_", schema$event_types) else character(), "risk", "rmst")
  arms <- vector("list", 2L); names(arms) <- c("0", "1")
  for (arm in 0:1) {
    rates <- array(0, c(M, S, L, ng)); censor_rates <- numeric(ng)
    for (f in seq_len(nf)) for (a in seq_len(nx)) {
      g <- (f - 1L) * nx + a
      ii <- which(folds[u] != f & x[u] == a & w[u] == arm)
      risk_time <- sum(exposure[ii])
      if (!length(ii) || risk_time <= 0) .rc_error("causal_training", "Each arm needs training exposure in every fold/stratum.")
      censor_rates[g] <- if (complete) 0 else (sum(event[ii] == M) + .5) / (risk_time + 1)
      for (k in seq_len(M)) {
        if (working == "pooled" || K == 0L) lambda <- rep((sum(event[ii] == k - 1L) + .5) / (risk_time + 1), S) else {
          tt <- vapply(seq_len(S), function(s) sum(exposure[ii[state[ii] == s]]), numeric(1)) + 1 / S
          cc <- tabulate(state[ii[event[ii] == k - 1L]], S) + .5 / S
          X <- cbind(1, states); init <- c(log(sum(cc) / sum(tt)), rep(0, K))
          fn <- function(b) {eta <- as.vector(X %*% b); sum(tt * exp(eta) - cc * eta) + 1e-6 * sum(b[-1]^2)}
          gr <- function(b) as.vector(crossprod(X, tt * exp(as.vector(X %*% b)) - cc)) + c(0, 2e-6 * b[-1])
          opt <- stats::optim(init, fn, gr, method = "L-BFGS-B", lower = rep(-10, K + 1), upper = rep(5, K + 1), control = list(maxit = 200))
          if (opt$convergence != 0L) .rc_error("causal_nuisance_convergence", "A training intensity fit did not converge.")
          lambda <- exp(as.vector(X %*% opt$par))
        }
        if (k %in% nonrecurrent) lambda[states[, k] > 0] <- 0
        for (b in seq_len(L)) rates[k, , b, g] <- lambda
      }
    }
    if (!is.null(spec$initial_rates)) {
      supplied <- spec$initial_rates[[as.character(arm)]]
      if (!is.numeric(supplied) || !identical(dim(supplied), dim(rates))) .rc_error("causal_nuisance", "Arm rates must have dimensions marks by living states by grid by fold-stratum.")
      rates <- supplied
    }
    if (!is.null(spec$censor_rates)) {
      supplied <- spec$censor_rates[[as.character(arm)]]
      if (!is.numeric(supplied) || length(supplied) != ng) .rc_error("causal_nuisance", "Arm censor_rates need one value per fold-stratum.")
      censor_rates <- supplied
    }
    if (any(!is.finite(rates)) || any(rates < 0) || any(!is.finite(censor_rates)) || any(censor_rates < 0))
      .rc_error("causal_nuisance", "Intensities and censor hazards must be finite and nonnegative.")
    for (k in nonrecurrent) if (any(rates[k, states[, k] > 0, , , drop = FALSE] != 0))
      .rc_error("causal_eligibility", "Supplied rates violate nonrecurrent eligibility.")
    if (complete && any(censor_rates != 0)) .rc_error("causal_followup", "Complete follow-up requires zero censor hazards.")
    if (exp(-max(censor_rates) * horizon) < min_G) .rc_error("causal_positivity", "Censor survival falls below the declared minimum; no clipping is applied.")
    arm_e <- if (arm == 1L) e else 1 - e
    selected <- w[u] == arm
    z <- .rc_core_call(rc_cpp_causal_target(as.numeric(rates), censor_rates, arm_e, nxt,
          as.integer(group - 1L), records[selected, , drop = FALSE], as.integer(L), horizon, as.integer(maxit), tolerance))
    z$initial_rates <- rates; z$rates <- array(z$rates, dim(rates))
    z$values <- array(z$values, c(J, S + 1L, L + 1L, ng))
    z$censor_rates <- censor_rates; z$propensity <- arm_e
    names(z$estimate) <- names(z$se) <- targets; colnames(z$EIF) <- targets
    counts <- c(vapply(seq_len(K), function(k) sum(dd$.rc_mark[selected] == k), numeric(1)), sum(dd$.rc_terminal[selected]))
    z$event_counts <- stats::setNames(counts, c(if (K) paste0("burden_", schema$event_types) else character(), "risk"))
    z$regular <- c(counts > 0, tail(counts, 1) > 0) & z$se > 0 & isTRUE(z$converged)
    z$EIF_mean <- colMeans(z$EIF); z$EIF_variance <- apply(z$EIF, 2, stats::var)
    arms[[as.character(arm)]] <- z
  }
  influence <- cbind(arms[["0"]]$EIF, arms[["1"]]$EIF)
  colnames(influence) <- c(paste0("arm0_", targets), paste0("arm1_", targets))
  covariance <- stats::cov(influence) / n
  effect_eif <- arms[["1"]]$EIF - arms[["0"]]$EIF
  effect <- arms[["1"]]$estimate - arms[["0"]]$estimate
  valid <- all(arms[[1]]$regular & arms[[2]]$regular)
  assumptions <- c("static baseline binary treatment; consistency and no interference",
    "conditional treatment exchangeability given declared baseline strata (not empirically testable)",
    "positive treatment probabilities and conditionally independent right censoring with positive survival",
    "independent subjects; sufficient retained count state; no simultaneous distinct jumps",
    "correct event intensities OR correct treatment-probability times censor-survival product for consistency",
    "EIF asymptotics require nuisance convergence, product-rate and targeting regularity conditions",
    "finite-state/time-grid approximation must be assessed; no finite-sample coverage guarantee")
  ans <- list(arms = arms, estimate = effect, EIF = effect_eif,
      covariance = stats::cov(effect_eif) / n, se = sqrt(diag(stats::cov(effect_eif) / n)),
      joint_EIF = influence, joint_covariance = covariance, n = n, data = d, schema = schema,
      treatment = treatment, censor = censor, target_names = targets, horizon = horizon,
      time_unit = time_unit, formula = formula, folds = folds, group = group,
      subject_stratum = x, strata = levels, strata_names = strata, subject_treatment = w,
      training_subjects = training, states = states, next_state = nxt, grid = L,
      targeting = "baseline_causal", spec = spec, assumptions = assumptions,
      convergence = list(valid = valid, score_solved = all(vapply(arms, function(z) z$converged, logical(1))),
          code = if (valid) "regular_asymptotic" else "inspect_arm_score_or_boundary",
          score = max(vapply(arms, function(z) z$score, numeric(1)))),
      diagnostics = list(minimum_arm_propensity = min(e, 1 - e),
          maximum_inverse_treatment_weight = max(1 / e, 1 / (1 - e)),
          minimum_censor_survival = min(vapply(arms, function(z) exp(-max(z$censor_rates) * horizon), numeric(1))),
          interpolation_error = max(vapply(arms, function(z) z$interpolation_error, numeric(1))),
          propensity_source = if (is.null(propensity)) "cross-fitted stratum proportions with 0.5 pseudocount" else "analyst-supplied known or out-of-fold probabilities"),
      audit = list(call = call, seed = seed, assumptions = assumptions))
  ans$audit$numerical_error <- ans$diagnostics$interpolation_error
  ans$audit$minimum_censor_survival <- ans$diagnostics$minimum_censor_survival
  ans$audit$captured_at_utc <- .rc_audit_time()
  ans$audit$fit_session <- .rc_audit_session()
  ans$audit$data_fingerprint <- .rc_audit_data_fingerprint(d)
  class(ans) <- c("recurcast_causal", "recurcast")
  if (!valid) warning("Some causal coordinates did not satisfy score, event-count or variance criteria; inspect diagnostics. Regular intervals are withheld for those coordinates.", call. = FALSE)
  ans
}

rc_causal_contrast <- function(object, contrast = c("difference", "ratio", "log_ratio"),
                               level = .95, interval = c("wald", "multiplier", "none"),
                               draws = 2000L, seed = 1L) {
  if (!inherits(object, "recurcast_causal")) .rc_error("causal_object", "Use an rc_causal fit.")
  contrast <- match.arg(contrast); interval <- match.arg(interval)
  .rc_scalar(level, "confidence level", 0, 1, TRUE)
  a <- object$arms[["1"]]$estimate; b <- object$arms[["0"]]$estimate
  D1 <- object$arms[["1"]]$EIF; D0 <- object$arms[["0"]]$EIF
  regular <- object$arms[["0"]]$regular & object$arms[["1"]]$regular
  if (contrast == "difference") {estimate <- a - b; influence <- D1 - D0} else {
    if (any(a <= 0 | b <= 0)) .rc_error("causal_ratio", "Ratio and log-ratio contrasts require strictly positive arm means.")
    influence <- sweep(D1, 2, a, "/") - sweep(D0, 2, b, "/")
    estimate <- if (contrast == "ratio") a / b else log(a / b)
  }
  cov <- stats::cov(influence) / object$n; se <- sqrt(diag(cov))
  valid <- regular & is.finite(se) & se > 0
  lower <- upper <- rep(NA_real_, length(estimate)); critical <- stats::qnorm((1 + level) / 2)
  if (interval == "multiplier" && any(valid)) {
    .rc_scalar(draws, "multiplier draws", 100, 100000)
    if (draws != floor(draws)) .rc_error("causal_draws", "draws must be integer.")
    # This small covariance factor gives the same conditional Gaussian law as
    # subject-level multipliers without allocating draws by cohort-size weights.
    spectrum <- eigen(cov[valid, valid, drop = FALSE], symmetric = TRUE)
    root <- sweep(spectrum$vectors, 2, sqrt(pmax(0, spectrum$values)), "*")
    perturb <- .rc_target_rng(seed, matrix(stats::rnorm(draws * sum(valid)), draws, sum(valid)) %*% t(root))
    maxima <- apply(abs(sweep(perturb, 2, se[valid], "/")), 1, max)
    critical <- unname(stats::quantile(maxima, level, type = 8))
  }
  if (interval != "none") {
    location <- if (contrast == "ratio") log(estimate) else estimate
    lower[valid] <- location[valid] - critical * se[valid]
    upper[valid] <- location[valid] + critical * se[valid]
    if (contrast == "ratio") {lower <- exp(lower); upper <- exp(upper)}
  }
  result <- data.frame(target = object$target_names, contrast = contrast, arm0 = unname(b), arm1 = unname(a),
      estimate = unname(estimate), std_error = unname(if (contrast == "ratio") estimate * se else se),
      lower = lower, upper = upper, inference_status = if (interval == "none") "not_requested" else
        ifelse(valid, if (interval == "multiplier") "simultaneous_asymptotic_regular_coordinates" else "pointwise_asymptotic", "withheld_nonregular"))
  attr(result, "EIF") <- if (contrast == "ratio") sweep(influence, 2, estimate, "*") else influence
  attr(result, "covariance") <- stats::cov(attr(result, "EIF")) / object$n
  attr(result, "critical_value") <- critical
  attr(result, "assumptions") <- object$assumptions
  result
}

predict.recurcast_causal <- function(object, contrast = "difference", ...) rc_causal_contrast(object, contrast, ...)
tidy.recurcast_causal <- function(x, ...) rc_causal_contrast(x, ...)
coef.recurcast_causal <- function(object, ...) object$estimate
vcov.recurcast_causal <- function(object, ...) object$covariance
logLik.recurcast_causal <- function(object, ...) .rc_error("causal_loglik", "The cross-fitted causal fit does not define one joint parametric maximized log likelihood.")
glance.recurcast_causal <- function(x, ...) data.frame(n = x$n, horizon = x$horizon,
    valid = x$convergence$valid, score = x$convergence$score,
    minimum_propensity = x$diagnostics$minimum_arm_propensity,
    minimum_censor_survival = x$diagnostics$minimum_censor_survival,
    interpolation_error = x$diagnostics$interpolation_error)
summary.recurcast_causal <- function(object, ...) list(estimates = rc_causal_contrast(object, ...),
    diagnostics = object$diagnostics, assumptions = object$assumptions, convergence = object$convergence)
print.recurcast_causal <- function(x, ...) {
  cat("Baseline-treatment causal extension (conditional on stated identification assumptions)\n")
  print(rc_causal_contrast(x, ...), row.names = FALSE)
  invisible(x)
}
