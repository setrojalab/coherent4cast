.cs_inference_stop <- function(message, code = "invalid_inference_input") {
  stop(structure(list(message = message, call = NULL, code = code),
                 class = c("coherentSurv_error", "error", "condition")))
}
.cs_scalar <- function(x, name, lower = -Inf, upper = Inf, open = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      if (open) x <= lower || x >= upper else x < lower || x > upper)
    .cs_inference_stop(paste("Invalid", name))
  invisible(x)
}
.cs_core_call <- function(expr) {
  tryCatch(expr, error = function(e) {
    if (inherits(e, "coherentSurv_error")) stop(e)
    .cs_inference_stop(conditionMessage(e), "numerical_inference_failure")
  })
}

#' Cross-fitted targeting of a finite path law
#' @export
cs_target <- function(path, recorded, strata, features, folds = NULL,
                      initial_laws = NULL, recording = NULL, control = list()) {
  call <- match.call(); n <- length(path); features <- as.matrix(features)
  if (!is.numeric(path) || n < 2L || length(recorded) != n || length(strata) != n ||
      anyNA(recorded) || anyNA(strata) || !all(recorded %in% c(0, 1)) ||
      !is.numeric(features) || any(!is.finite(features)) || !nrow(features) || !ncol(features))
    .cs_inference_stop("One path, recording indicator and observed stratum are required per subject")
  if (length(setdiff(names(control), c("folds", "smoothing", "recording_floor", "maxit", "tolerance"))))
    .cs_inference_stop("Unknown targeting control")
  get_control <- function(name, default) if (is.null(control[[name]])) default else control[[name]]
  nfold <- get_control("folds", 5L)
  .cs_scalar(nfold, "number of folds", 2, Inf)
  if (nfold != as.integer(nfold)) .cs_inference_stop("number of folds must be integer")
  if (is.null(folds)) folds <- (seq_len(n) - 1L) %% min(nfold, n) + 1L
  if (length(folds) != n || anyNA(folds)) .cs_inference_stop("Invalid subject folds")
  folds <- match(folds, unique(folds)); nf <- max(folds)
  if (nf < 2L) .cs_inference_stop("At least two cross-fitting folds are required")
  levels <- sort(unique(strata)); x <- match(strata, levels); ns <- length(levels)
  j <- nrow(features); d <- ncol(features); obs <- which(recorded == 1)
  if (!length(obs) || anyNA(path[obs]) || any(!is.finite(path[obs])) ||
      any(path[obs] != as.integer(path[obs])) || any(path[obs] < 1 | path[obs] > j))
    .cs_inference_stop("Recorded paths must be integer indices in the complete support")
  smoothing <- get_control("smoothing", .5)
  floor <- get_control("recording_floor", .05)
  .cs_scalar(smoothing, "smoothing", 0, Inf, TRUE)
  .cs_scalar(floor, "recording floor", 0, .5, TRUE)
  g <- nf * ns; count <- matrix(0, g, j); q0 <- count; pi <- numeric(g); size <- integer(g)
  group <- (folds - 1L) * ns + x
  for (f in seq_len(nf)) for (a in seq_len(ns)) {
    id <- (f - 1L) * ns + a; tr <- folds != f & x == a; va <- which(group == id)
    if (!any(tr)) .cs_inference_stop("A training fold has no subjects in an observed stratum")
    observed_train <- which(tr & recorded == 1); size[id] <- length(va)
    count[id, ] <- tabulate(path[va[recorded[va] == 1]], j)
    q0[id, ] <- (tabulate(path[observed_train], j) + smoothing / j) /
      (length(observed_train) + smoothing)
    pi[id] <- max(floor, min(.99, (sum(recorded[tr]) + .5) / (sum(tr) + 1)))
  }
  if (!is.null(initial_laws)) {
    initial_laws <- as.matrix(initial_laws)
    if (!all(dim(initial_laws) == c(g, j))) .cs_inference_stop("initial_laws must have fold-by-stratum rows and complete-path columns")
    q0 <- initial_laws
  }
  if (!is.null(recording)) {
    if (length(recording) != g) .cs_inference_stop("recording needs one probability per fold-by-stratum group")
    pi <- recording
  }
  tolerance <- get_control("tolerance", 1e-6 / sqrt(n)); maxit <- get_control("maxit", 60L)
  .cs_scalar(tolerance, "targeting tolerance", 0, Inf, TRUE)
  .cs_scalar(maxit, "maximum iterations", 1, 10000)
  z <- .cs_core_call(cs_cpp_target(count, q0, features, pi, n, as.integer(maxit), tolerance))
  estimate <- as.vector(crossprod(size, z$mu) / n)
  influence <- sweep(z$mu[group, , drop = FALSE], 2L, estimate)
  influence[obs, ] <- influence[obs, , drop = FALSE] +
    (features[path[obs], , drop = FALSE] - z$mu[group[obs], , drop = FALSE]) / pi[group[obs]]
  covariance <- crossprod(influence) / (n * (n - 1))
  probs <- matrix(0, ns, j)
  for (a in seq_len(ns)) {
    ix <- seq.int(a, g, by = ns)
    probs[a, ] <- as.vector(crossprod(size[ix], z$q[ix, , drop = FALSE])) / sum(size[ix])
  }
  names(estimate) <- colnames(features)
  mixture <- as.vector(crossprod(tabulate(x, ns) / n, probs))
  direct <- as.vector(crossprod(mixture, features))
  if (max(abs(estimate - direct)) > 1e-10) .cs_inference_stop("Common-law target check failed", "law_incoherence")
  structure(list(estimate = estimate, se = sqrt(diag(covariance)), covariance = covariance,
                 EIF = influence, EIF_mean = colMeans(influence), EIF_variance = apply(influence, 2L, stats::var),
                 probabilities = probs, mixture = mixture, group_probabilities = z$q, initial_laws = q0,
                 folds = folds, strata = levels, recording = pi, subject_recording = pi[group],
                 epsilon = z$epsilon, iterations = z$iterations,
                 convergence = list(valid = z$converged && max(abs(colMeans(influence))) < tolerance,
                                    score = max(abs(colMeans(influence))), code = if (z$converged) "converged" else "target_score_not_solved"),
                 assumptions = c("finite complete-path support", "MAR conditional on observed strata", "positive recording probabilities", "subject-separated nuisance fitting"),
                 call = call), class = "cs_target")
}

#' Sharp finite-support recording-sensitivity bounds
#' @export
cs_bounds <- function(p, numerator, denominator = rep(1, length(numerator)),
                      kappa, recording_fraction, counts = NULL, alpha = .05,
                      observation_map = seq_along(numerator)) {
  if (length(numerator) > 64L) return(cs_bounds_sparse(p = p,
      numerator = numerator, denominator = denominator, kappa = kappa,
      recording_fraction = recording_fraction, counts = counts, alpha = alpha,
      observation_map = observation_map))
  .cs_scalar(alpha, "alpha", 0, 1, TRUE)
  if (!is.numeric(p) || !is.numeric(numerator) || !is.numeric(denominator) ||
      length(numerator) != length(denominator) || length(observation_map) != length(numerator) ||
      anyNA(observation_map) || any(observation_map != as.integer(observation_map)) ||
      any(observation_map < 1 | observation_map > length(p)))
    .cs_inference_stop("Invalid support or coarsening-fibre map")
  .cs_scalar(kappa, "kappa", 1, Inf)
  .cs_scalar(recording_fraction, "recording_fraction", 0, 1, TRUE)
  rho <- recording_fraction; lower <- upper <- rho * p
  fallback_messages<-character()
  run <- function(lo, hi) tryCatch(.cs_core_call(cs_cpp_bounds(p, numerator, denominator, kappa, rho, lo, hi, as.integer(observation_map - 1L))),error=function(e) {
    z<-.cs_core_call(cpp_cs_bounds_sparse(p,numerator,denominator,kappa,rho,lo,hi,as.integer(observation_map-1L),1e-10,1000L))
    fallback_messages<<-c(fallback_messages,conditionMessage(e));z$solver_fallback<-"Structural sparse solver after dense numerical failure";z
  })
  plugin <- run(lower, upper); outer <- NULL; region <- NULL
  if (!is.null(counts)) {
    if (length(counts) != length(p) || any(!is.finite(counts)) || any(counts < 0 | counts != floor(counts)) || sum(counts) < 1 || max(abs(counts / sum(counts) - p)) > 1e-9)
      .cs_inference_stop("counts must be nonnegative integers whose proportions equal p")
    n <- sum(counts); m <- length(counts); level <- alpha / (2 * m)
    lo <- ifelse(counts == 0, 0, stats::qbeta(level, pmax(counts, 1), n - counts + 1))
    hi <- ifelse(counts == n, 1, stats::qbeta(1 - level, counts + 1, pmax(n - counts, 1)))
    region <- cbind(lower = rho * lo, upper = rho * hi)
    min_den <- vapply(seq_len(m), function(h) min(denominator[observation_map == h]), numeric(1))
    if (sum(region[, 1L] * min_den) > 0) outer <- run(region[, 1L], region[, 2L]) else {
      if (any(numerator[denominator == 0] != 0)) .cs_inference_stop("Cannot form a bounded conditional support envelope")
      values <- numerator[denominator > 0] / denominator[denominator > 0]
      outer <- list(lower = min(values), upper = max(values), method = "support envelope; denominator lower bound is zero")
    }
  }
  structure(list(lower = plugin$lower, upper = plugin$upper, plugin = plugin,
                 outer = outer, recorded_confidence_region = region, kappa = kappa,
                 recording_fraction = rho, observation_map = observation_map, alpha = alpha,
                 assumptions = c("known recording fraction", "finite complete support and declared coarsening fibres", "marginal recording odds bound", "positive conditioning denominator"),
                 solver_fallback_messages=fallback_messages,
                 method = if(length(fallback_messages))"Charnes-Cooper bounds; recorded structural-sparse fallback after dense numerical failure"else "Charnes-Cooper transformation with checked two-phase simplex"), class = "cs_bounds")
}

#' Weighted conformal candidate-set inversion
#' @export
cs_conformal <- function(scores, candidate_scores, weights = rep(1, length(scores)),
                         candidate_weights = rep(1, length(candidate_scores)),
                         weight_lower = NULL, weight_upper = NULL,
                         candidate_weight_upper = NULL, alpha = .05, band_error = 0,
                         weight_status = c("known", "estimated", "certified")) {
  weight_status <- match.arg(weight_status); .cs_scalar(alpha, "alpha", 0, 1, TRUE)
  .cs_scalar(band_error, "band error", 0, alpha)
  if (band_error >= alpha) .cs_inference_stop("band_error must be smaller than alpha")
  if (weight_status == "certified") {
    if (is.null(weight_lower) || is.null(weight_upper) || is.null(candidate_weight_upper))
      .cs_inference_stop("Certified inversion requires all calibration and candidate weight bounds")
    lo <- weight_lower; hi <- weight_upper; cw <- candidate_weight_upper
  } else {
    if (!is.null(weight_lower) || !is.null(weight_upper) || !is.null(candidate_weight_upper) || band_error != 0)
      .cs_inference_stop("Weight intervals and a band-error allocation require weight_status='certified'")
    lo <- hi <- weights; cw <- candidate_weights
  }
  pvalue <- .cs_core_call(cs_cpp_conformal(scores, candidate_scores, lo, hi, cw))
  structure(list(candidates = data.frame(score = candidate_scores, p_value = pvalue,
                                          included = pvalue > alpha - band_error),
                 alpha = alpha, rank_alpha = alpha - band_error, band_error = band_error,
                 weight_status = weight_status, coverage_guarantee = if (weight_status == "estimated") "No exact guarantee from estimated weights alone" else "Marginal guarantee conditional on stated sampling and weight assumptions",
                 assumptions = c("independent training and calibration", "appropriate weighted exchangeability", "candidate mass included", if (weight_status == "certified") "simultaneous weight band has stated error probability")), class = "cs_conformal")
}

#' Exact known-dispersion gamma-frailty pivot
#' @export
cs_pivot <- function(times, dispersion, horizons, alpha = .05) {
  .cs_scalar(dispersion, "dispersion", 0, Inf); .cs_scalar(alpha, "alpha", 0, 1, TRUE)
  if (!length(horizons) || any(!is.finite(horizons)) || any(horizons <= 0)) .cs_inference_stop("horizons must be finite and positive")
  quantiles <- stats::qgamma(c(alpha / 2, 1 - alpha / 2), shape = length(times))
  rate <- .cs_core_call(cs_cpp_pivot(times, dispersion, quantiles))
  risk <- function(lambda) if (dispersion == 0) -expm1(-lambda * horizons) else -expm1(-log1p(dispersion * lambda * horizons) / dispersion)
  list(rate_interval = stats::setNames(rate, c("lower", "upper")),
       risk = data.frame(horizon = horizons, lower = risk(rate[1L]), upper = risk(rate[2L])),
       alpha = alpha, dispersion = dispersion,
       scope = "Known dispersion, constant baseline and uncensored independent terminal times; simultaneous over all positive horizons")
}

#' Cornish-Fisher standardized-mean quantiles
#' @export
cs_cornish_fisher <- function(probability, n, skewness = 2, excess_kurtosis = 6, order = 2L) {
  if (!length(probability) || any(!is.finite(probability)) || any(probability <= 0 | probability >= 1)) .cs_inference_stop("probability must be inside (0,1)")
  .cs_scalar(n, "sample size", 1, Inf); .cs_scalar(order, "order", 0, 2)
  if (order != as.integer(order)) .cs_inference_stop("order must be 0, 1 or 2")
  .cs_core_call(cs_cpp_cf(stats::qnorm(probability), n, skewness, excess_kurtosis, as.integer(order)))
}

#' Observable-iteration identification diagnostic
#' @export
cs_identify <- function(iterate, survival, contrast, iterations = 20L) {
  if (!is.function(iterate)) .cs_inference_stop("iterate must be the supplied observable contraction")
  .cs_scalar(survival, "survival", 0, 1, TRUE); .cs_scalar(contrast, "contrast", 1, Inf, TRUE)
  .cs_scalar(iterations, "iterations", 1, 10000)
  if (iterations != as.integer(iterations)) .cs_inference_stop("iterations must be integer")
  out <- data.frame(iteration = seq_len(iterations), survival = NA_real_, baseline = NA_real_)
  value <- survival
  for (i in seq_len(iterations)) {
    value <- iterate(value)
    .cs_scalar(value, "iterated survival", 0, 1)
    recovered <- exp(i * log(contrast)) * (1 - value)
    if (!is.finite(recovered)) .cs_inference_stop("Identification iteration exceeds numerical range")
    out$survival[i] <- value; out$baseline[i] <- recovered
  }
  attr(out, "scope") <- "Population contraction diagnostic; common mean-one frailty independent of covariates and nonunit multiplicative contrast required"
  attr(out, "rounding_limit") <- any(out$survival == 1)
  out
}

#' Frequentist Gaussian simulation for simultaneous nonlinear forecast bands
#' @export
cs_band <- function(theta, covariance, forecast, draws = 1000L, alpha = .05, seed = 1L,
                    scale = c("logit", "identity")) {
  scale <- match.arg(scale); .cs_scalar(alpha, "alpha", 0, 1, TRUE)
  .cs_scalar(draws, "draws", 2, 1000000)
  if (draws != as.integer(draws) || !is.function(forecast) || !is.numeric(theta) || any(!is.finite(theta))) .cs_inference_stop("Invalid forecast, theta or draw count")
  covariance <- as.matrix(covariance)
  if (!identical(dim(covariance), c(length(theta), length(theta))) || any(!is.finite(covariance)) || max(abs(covariance - t(covariance))) > 1e-9) .cs_inference_stop("covariance must be a finite symmetric parameter covariance matrix")
  eig <- eigen(covariance, symmetric = TRUE)
  if (min(eig$values) < -1e-10 * max(1, max(eig$values))) .cs_inference_stop("covariance is not positive semidefinite")
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (present) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (present) assign(".Random.seed", old, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed); p <- length(theta); normals <- matrix(stats::rnorm(draws * p), draws, p)
  parameter_draws <- sweep(normals %*% (sqrt(pmax(eig$values, 0)) * t(eig$vectors)), 2L, theta, "+")
  point <- forecast(theta)
  if (!is.numeric(point) || !length(point) || any(!is.finite(point))) .cs_inference_stop("forecast returned an invalid point estimate")
  values <- t(matrix(vapply(seq_len(draws), function(i) forecast(parameter_draws[i, ]), numeric(length(point))), nrow = length(point), ncol = draws))
  if (any(!is.finite(values))) .cs_inference_stop("A forecast draw was nonfinite; no draws were silently removed", "forecast_draw_failure")
  if (scale == "logit") {
    if (any(point <= 0 | point >= 1) || any(values <= 0 | values >= 1)) .cs_inference_stop("logit bands require all point and draw probabilities strictly inside (0,1)", "forecast_boundary")
    center <- stats::qlogis(point); transformed <- stats::qlogis(values); inverse <- stats::plogis
  } else {center <- point; transformed <- values; inverse <- identity}
  se <- apply(transformed, 2L, stats::sd)
  if (any(se <= 0)) .cs_inference_stop("A forecast coordinate has zero Gaussian draw variance")
  max_t <- apply(abs(sweep(sweep(transformed, 2L, center), 2L, se, "/")), 1L, max)
  critical <- as.numeric(stats::quantile(max_t, 1 - alpha, type = 1))
  list(estimate = point, lower = inverse(center - critical * se), upper = inverse(center + critical * se),
       transformed_se = se, critical = critical, alpha = alpha, draws = draws, seed = seed,
       scale = scale, method = "Frequentist sandwich-Gaussian max-t approximation",
       scope = "Uses the supplied estimator and covariance; neither an exact confidence projection nor a Bayesian posterior")
}
