#' Construct a compressed sparse row generator from directed transitions
#' @export
cs_generator_csr <- function(n, from = integer(), to = integer(), rate = numeric()) {
 .cs_kernel_guard({
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n) || n > .Machine$integer.max - 1L)
    stop("n must be a positive integer", call. = FALSE)
  if (!is.numeric(from) || !is.numeric(to) || !is.numeric(rate) ||
      length(from) != length(to) || length(from) != length(rate) ||
      any(!is.finite(from)) || any(!is.finite(to)) || any(!is.finite(rate)) ||
      any(from < 1 | from > n | from != floor(from)) ||
      any(to < 1 | to > n | to != floor(to)) || any(from == to) || any(rate < 0))
    stop("supply equally sized 1-based, off-diagonal from/to indices and nonnegative finite rates", call. = FALSE)
  keep <- rate > 0; from <- as.integer(from[keep]); to <- as.integer(to[keep]); rate <- rate[keep]
  o <- order(from, to); from <- from[o]; to <- to[o]; rate <- rate[o]
  counts <- tabulate(from, nbins = as.integer(n))
  structure(list(n = as.integer(n), offsets = as.integer(c(0, cumsum(counts))),
                 columns = to - 1L, rates = rate,
                 storage = "CSR off-diagonal rates; diagonal equals negative row sum"),
            class = "cs_generator_csr")
 })
}

#' Sparse ordered uniformization with integrated rewards
#' @export
cs_propagate_sparse <- function(probability, generator, duration, reward = NULL,
                               tolerance = 1e-12, max_terms = 10000L,
                               max_chunks = 100000L) {
 .cs_kernel_guard({
  if (!is.numeric(probability) || !length(probability) || any(!is.finite(probability)))
    stop("probability must be a finite numeric vector", call. = FALSE)
  generators <- if (inherits(generator, "cs_generator_csr")) list(generator) else generator
  if (!is.list(generators) || !length(generators) || !is.numeric(duration) ||
      length(duration) != length(generators) || any(!is.finite(duration) | duration < 0))
    stop("supply one nonnegative duration per ordered CSR generator", call. = FALSE)
  n <- length(probability)
  if (any(!vapply(generators, function(q) inherits(q, "cs_generator_csr") && identical(q$n, as.integer(n)), logical(1))))
    stop("each CSR generator must match the initial probability vector", call. = FALSE)
  if (is.null(reward)) reward <- matrix(numeric(), nrow = n, ncol = 0L)
  rewards <- if (is.list(reward)) reward else rep(list(reward), length(generators))
  rewards <- lapply(rewards, function(x) if (is.numeric(x) && !is.matrix(x)) matrix(x, ncol = 1L) else x)
  if (length(rewards) != length(generators) ||
      any(!vapply(rewards, function(x) is.matrix(x) && is.numeric(x) && nrow(x) == n, logical(1))) ||
      length(unique(vapply(rewards, ncol, integer(1)))) != 1L)
    stop("reward matrices must share dimensions state by reward", call. = FALSE)
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1 ||
      length(max_terms) != 1L || !is.finite(max_terms) || max_terms < 2 || max_terms != floor(max_terms) ||
      length(max_chunks) != 1L || !is.finite(max_chunks) || max_chunks < 1 || max_chunks != floor(max_chunks))
    stop("invalid sparse propagation tolerance or budget", call. = FALSE)
  total <- numeric(ncol(rewards[[1L]])); reward_error <- total
  error <- roundoff <- terms <- chunks <- workspace <- 0
  for (i in seq_along(generators)) {
    q <- generators[[i]]
    z <- cpp_cs_propagate_sparse(probability, q$offsets, q$columns, q$rates,
       duration[i], rewards[[i]], tolerance / length(generators), as.integer(max_terms), as.integer(max_chunks))
    if (length(total)) reward_error <- reward_error + z$reward_error_bound +
      error * duration[i] * apply(abs(rewards[[i]]), 2L, max)
    total <- total + z$reward; probability <- z$probability
    error <- error + z$error_bound; roundoff <- roundoff + z$roundoff_bound
    terms <- terms + z$terms; chunks <- chunks + z$chunks; workspace <- max(workspace, z$workspace_bytes)
  }
  list(probability = probability, reward = total, error_bound = error,
       reward_error_bound = reward_error, roundoff_bound = min(1, roundoff),
       terms = terms, chunks = chunks, workspace_bytes = workspace,
       method = "CSR row-vector uniformization; ordered segments; certified series truncation")
 })
}

#' Scalable sharp finite-support sensitivity bounds
#' @export
cs_bounds_sparse <- function(p, numerator, denominator = rep(1, length(numerator)),
                             kappa, recording_fraction, counts = NULL, alpha = .05,
                             observation_map = seq_along(numerator),
                             tolerance = 1e-10, max_iterations = 1000L) {
 .cs_kernel_guard({
  .cs_scalar(alpha, "alpha", 0, 1, TRUE)
  if (!is.numeric(p) || !is.numeric(numerator) || !is.numeric(denominator) ||
      length(numerator) != length(denominator) || length(observation_map) != length(numerator) ||
      anyNA(observation_map) || any(observation_map != as.integer(observation_map)) ||
      any(observation_map < 1 | observation_map > length(p)))
    stop("invalid complete support or observation-fibre map", call. = FALSE)
  .cs_scalar(kappa, "kappa", 1, Inf)
  .cs_scalar(recording_fraction, "recording_fraction", 0, 1, TRUE)
  .cs_scalar(tolerance, "tolerance", 0, Inf, TRUE)
  .cs_scalar(max_iterations, "max_iterations", 1, .Machine$integer.max)
  if (max_iterations != floor(max_iterations)) stop("max_iterations must be integer", call. = FALSE)
  rho <- recording_fraction
  run <- function(lo, hi) cpp_cs_bounds_sparse(p, numerator, denominator, kappa, rho,
      lo, hi, as.integer(observation_map - 1L), tolerance, as.integer(max_iterations))
  plugin <- run(rho * p, rho * p); outer <- region <- NULL
  if (!is.null(counts)) {
    if (length(counts) != length(p) || any(!is.finite(counts)) ||
        any(counts < 0 | counts != floor(counts)) || sum(counts) < 1 ||
        max(abs(counts / sum(counts) - p)) > 1e-9)
      stop("counts must be nonnegative integers whose proportions equal p", call. = FALSE)
    n <- sum(counts); m <- length(counts); level <- alpha / (2 * m)
    lo <- ifelse(counts == 0, 0, stats::qbeta(level, pmax(counts, 1), n - counts + 1))
    hi <- ifelse(counts == n, 1, stats::qbeta(1 - level, counts + 1, pmax(n - counts, 1)))
    region <- cbind(lower = rho * lo, upper = rho * hi)
    # One pass avoids an M by J scan with tens of thousands of support atoms.
    min_den <- rep(Inf, m)
    for (k in seq_along(denominator)) min_den[observation_map[k]] <- min(min_den[observation_map[k]], denominator[k])
    if (sum(region[, 1L] * min_den) > 0) outer <- run(region[, 1L], region[, 2L]) else {
      if (any(numerator[denominator == 0] != 0)) stop("cannot form a bounded conditional support envelope", call. = FALSE)
      values <- numerator[denominator > 0] / denominator[denominator > 0]
      outer <- list(lower = min(values), upper = max(values),
                    method = "support envelope; denominator lower bound is zero")
    }
  }
  structure(list(lower = plugin$lower, upper = plugin$upper, plugin = plugin,
      outer = outer, recorded_confidence_region = region, kappa = kappa,
      recording_fraction = rho, observation_map = observation_map, alpha = alpha,
      assumptions = c("known recording fraction", "finite complete support and declared coarsening fibres",
                      "marginal recording odds bound", "positive conditioning denominator"),
      method = "Exact structural reduction of Charnes-Cooper bounds; sorted fibre allocation and fractional dual certificates"),
      class = "cs_bounds")
 })
}
