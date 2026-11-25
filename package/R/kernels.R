.rc_kernel_guard <- function(code) {
  tryCatch(force(code), error = function(e) {
    if (inherits(e, "coherent4cast_error")) stop(e)
    stop(structure(list(message = conditionMessage(e), call = NULL,
                        code = "KERNEL_INPUT_OR_NUMERICAL_ERROR"),
                   class = c("coherent4cast_error", "error", "condition")))
  })
}

#' Exact shared-gamma likelihood and analytic derivatives
#'
#' Evaluate the finite-basis complete-history likelihood used in the theory
#' manuscript. Expanded rows represent eligible marked-event exposures.
#' @param theta Numeric coefficient vector followed by log frailty variance when
#'   `phi_fixed` is `NULL`.
#' @param unit Contiguous positive integer subject indices, one per expanded row.
#' @param exposure Nonnegative exposure duration for each row.
#' @param design Numeric design matrix with one row per expanded row.
#' @param event Binary event indicator for each expanded row.
#' @param phi_fixed Nonnegative fixed frailty variance, or `NULL` to estimate it.
#' @param information Compute observed information, including each subject's matrix.
#' @return A list containing the log likelihood, negative-mean objective and
#'   gradient, subject log likelihoods and scores, cumulative hazards and counts.
#'   Information matrices are included when requested.
#' @export
rc_loglik <- function(theta, unit, exposure, design, event, phi_fixed = NULL,
                      information = TRUE) {
 .rc_kernel_guard({
  if (!is.matrix(design) || !is.numeric(design) || !ncol(design) ||
      !is.numeric(theta) || any(!is.finite(theta))) {
    stop("design must be a numeric matrix and theta must be finite", call. = FALSE)
  }
  if (!is.numeric(unit) || any(!is.finite(unit)) || any(unit < 1 | unit != floor(unit)) ||
      length(unit) != nrow(design) || !length(unit) ||
      !identical(sort(unique(as.integer(unit))), seq_len(max(unit)))) {
    stop("unit must contain contiguous positive integer subject indices", call. = FALSE)
  }
  if (!is.numeric(exposure) || length(exposure) != length(unit) ||
      !is.numeric(event) || length(event) != length(unit)) {
    stop("exposure and event must be numeric vectors matching unit", call. = FALSE)
  }
  if (!is.logical(information) || length(information) != 1L || is.na(information)) {
    stop("information must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(phi_fixed) && (!is.numeric(phi_fixed) || length(phi_fixed) != 1L ||
      !is.finite(phi_fixed) || phi_fixed < 0)) {
    stop("phi_fixed must be NULL or one finite nonnegative variance", call. = FALSE)
  }
  cpp_rc_loglik(theta, as.integer(unit), exposure, design, event,
               if (is.null(phi_fixed)) -1 else phi_fixed, information)
 })
}

# Aggregate already enumerated legal paths. Quadrature weights are integration
# weights and must not be normalized separately within observed subjects.
rc_coarsen <- function(loglik, scores = NULL, information = NULL, owner,
                       log_weights) {
 .rc_kernel_guard({
  if (is.list(loglik)) {
    scores <- loglik$scores
    information <- loglik$information_individual
    loglik <- loglik$subject_loglik
  }
  if (!is.numeric(loglik) || !is.matrix(scores) || !is.numeric(scores) ||
      nrow(scores) != length(loglik) || !is.numeric(owner) ||
      length(owner) != length(loglik) || any(!is.finite(owner)) ||
      any(owner < 1 | owner != floor(owner)) || !length(owner) ||
      !identical(sort(unique(as.integer(owner))), seq_len(max(owner)))) {
    stop("legal paths require compatible log likelihoods, scores and owner indices", call. = FALSE)
  }
  if (!is.numeric(log_weights) || length(log_weights) != length(owner)) {
    stop("log_weights must contain one log integration weight per legal path", call. = FALSE)
  }
  if (!is.null(information) && (!is.numeric(information) ||
      !identical(dim(information), c(ncol(scores), ncol(scores), nrow(scores))))) {
    stop("path information must have dimensions parameter by parameter by path", call. = FALSE)
  }
  cpp_rc_coarsen(loglik, scores, if (is.null(information)) numeric() else as.numeric(information),
                as.integer(owner), log_weights)
 })
}

#' Propagate a finite-state law by ordered uniformization
#'
#' A single generator or ordered list of generators defines a piecewise constant
#' process. A latent frailty value must remain fixed across its entire ordered
#' propagation; mixture averaging belongs outside this function.
#' @param probability Initial row probability vector.
#' @param generator Square generator matrix or ordered list of such matrices.
#' @param duration Segment duration or vector of durations.
#' @param reward Optional state by reward matrix (or vector for one reward), or
#'   list of matrices matching the segments.
#' @param tolerance Positive Poisson-series truncation tolerance, less than one.
#' @param max_terms Maximum terms per scaled uniformization calculation.
#' @return A list with the final probability, integrated reward, series truncation
#'   bounds, a separate floating-point roundoff bound, and computational counts.
#' @export
rc_propagate <- function(probability, generator, duration, reward = NULL,
                         tolerance = 1e-12, max_terms = 10000L) {
 .rc_kernel_guard({
  if (!is.numeric(probability) || !length(probability) || any(!is.finite(probability))) {
    stop("probability must be a finite numeric vector", call. = FALSE)
  }
  generators <- if (is.matrix(generator)) list(generator) else generator
  if (!is.list(generators) || !length(generators) ||
      !is.numeric(duration) || length(duration) != length(generators) ||
      any(!is.finite(duration) | duration < 0)) {
    stop("supply one nonnegative duration per ordered generator", call. = FALSE)
  }
  n <- length(probability)
  if (any(!vapply(generators, function(q) is.matrix(q) && is.numeric(q) &&
                  identical(dim(q), c(n, n)), logical(1)))) {
    stop("every generator must be a numeric square matrix matching probability", call. = FALSE)
  }
  if (is.null(reward)) reward <- matrix(numeric(), nrow = n, ncol = 0L)
  rewards <- if (is.list(reward)) reward else rep(list(reward), length(generators))
  if (length(rewards) != length(generators)) stop("one reward matrix is required per segment", call. = FALSE)
  rewards <- lapply(rewards, function(x) if (is.numeric(x) && !is.matrix(x)) matrix(x, ncol = 1L) else x)
  if (any(!vapply(rewards, function(r) is.matrix(r) && is.numeric(r) && nrow(r) == n, logical(1))) ||
      length(unique(vapply(rewards, ncol, integer(1)))) != 1L) {
    stop("reward matrices must share dimensions state by reward", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) ||
      tolerance <= 0 || tolerance >= 1 || !is.numeric(max_terms) || length(max_terms) != 1L ||
      !is.finite(max_terms) || max_terms < 2 || max_terms != floor(max_terms)) {
    stop("invalid truncation tolerance or maximum term count", call. = FALSE)
  }
  accumulated <- numeric(ncol(rewards[[1L]])); error <- 0; roundoff <- 0
  reward_error <- accumulated; terms <- squarings <- 0
  for (i in seq_along(generators)) {
    z <- cpp_rc_propagate(probability, generators[[i]], duration[i], rewards[[i]],
                          tolerance / length(generators), as.integer(max_terms))
    if (length(reward_error)) {
      reward_error <- reward_error + z$reward_error_bound +
        error * duration[i] * apply(abs(rewards[[i]]), 2L, max)
    }
    accumulated <- accumulated + z$reward
    probability <- z$probability
    error <- error + z$error_bound
    roundoff <- roundoff + z$roundoff_bound
    terms <- terms + z$terms; squarings <- squarings + z$squarings
  }
  list(probability = probability, reward = accumulated, error_bound = error,
       reward_error_bound = reward_error, roundoff_bound = min(1, roundoff),
       terms = terms, squarings = squarings)
 })
}
