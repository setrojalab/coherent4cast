# Scientific boundary cases supplement, rather than repeat, the oracle fixtures.

test_that("external nuisance laws retain subject folds and numerical failures remain visible", {
  a <- read.csv(test_path("fixtures", "target_subjects.csv"))
  features <- as.matrix(read.csv(test_path("fixtures", "target_features.csv")))
  base <- rc_target(a$path, a$recorded, a$strata, features, a$fold)
  supplied <- rc_target(a$path, a$recorded, a$strata, features, a$fold,
                        initial_laws = base$initial_laws, recording = base$recording)
  expect_equal(supplied$estimate, base$estimate, tolerance = 1e-13)
  expect_equal(supplied$EIF, base$EIF, tolerance = 1e-12)
  expect_identical(supplied$folds, base$folds)
  expect_equal(supplied$recording, base$recording)
  limited <- rc_target(a$path, a$recorded, a$strata, features, a$fold,
                        control = list(maxit = 1, tolerance = 1e-15))
  expect_false(limited$convergence$valid)
  expect_identical(limited$convergence$code, "target_score_not_solved")
  expect_gt(limited$convergence$score, 1e-15)
  expect_equal(rowSums(limited$probabilities), rep(1, nrow(limited$probabilities)), tolerance = 1e-13)
  expect_error(rc_target(a$path, a$recorded, a$strata, features, a$fold,
                         initial_laws = matrix(1 / 13, 9, 13)), class = "recurcast_error")
  wrong_mass <- base$initial_laws * .9
  expect_error(rc_target(a$path, a$recorded, a$strata, features, a$fold,
                         initial_laws = wrong_mass), "sum to one", class = "recurcast_error")
  zero_mass <- base$initial_laws; zero_mass[1, 1] <- 0
  expect_error(rc_target(a$path, a$recorded, a$strata, features, a$fold,
                         initial_laws = zero_mass), "strictly positive", class = "recurcast_error")
  expect_error(rc_target(a$path, a$recorded, a$strata, features, a$fold,
                         recording = .7), class = "recurcast_error")
  expect_error(rc_target(a$path, a$recorded, a$strata, features, a$fold,
                         control = list(silent_drop_failures = TRUE)), class = "recurcast_error")
  expect_error(rc_target(a$path, a$recorded, a$strata, features,
                         control = list(folds = 2.5)), class = "recurcast_error")
  expect_error(rc_target(a$path, a$recorded, a$strata, features, rep(1, nrow(a))), class = "recurcast_error")
  expect_error(rc_target(a$path, a$recorded, a$strata, features, rep(NA, nrow(a))), class = "recurcast_error")
  expect_error(rc_target(a$path, a$recorded, a$strata, features, 1:4), class = "recurcast_error")
  bad_path <- a$path; bad_path[which(a$recorded == 1)[1]] <- nrow(features) + 1
  expect_error(rc_target(bad_path, a$recorded, a$strata, features, a$fold), class = "recurcast_error")
  # A stratum occurring in only one fold has no training data for that validation stratum.
  expect_error(rc_target(c(1, 2, 1, 2), rep(1, 4), c(0, 0, 1, 1),
                         matrix(c(0, 1), 2), c(1, 1, 2, 2)), "training fold", class = "recurcast_error")
  # Two identical target columns cannot define two identifiable fluctuation directions.
  dependent <- cbind(features[, 1], features[, 1])
  expect_error(rc_target(a$path, a$recorded, a$strata, dependent, a$fold), "singular", class = "recurcast_error")
})

test_that("coarsening and conditioning restrictions are checked before reporting bounds", {
  expect_error(rc_bounds(c(.5, .5), c(0, 1, 0), kappa = 2, recording_fraction = .7,
                         observation_map = c(1, 2, 3)), class = "recurcast_error")
  expect_error(rc_bounds(c(.5, .5), c(0, 1), kappa = 2, recording_fraction = .7,
                         observation_map = c(1, 1)), class = "recurcast_error")
  expect_error(rc_bounds(c(.5, .5), c(0, 1), kappa = 2, recording_fraction = .7,
                         counts = c(2, 3)), "proportions equal p", class = "recurcast_error")
  expect_error(rc_bounds(c(.5, .5), c(0, 1), kappa = 2, recording_fraction = .7,
                         counts = c(.5, .5)), class = "recurcast_error")
  expect_error(rc_bounds(c(.5, .5), c(0, 1), kappa = 2, recording_fraction = 1), class = "recurcast_error")
  # A coarse category can hide every conditioning event; its lower denominator is zero.
  expect_error(rc_bounds(c(.5, .5), c(0, 1, 0), c(0, 1, 0), 2, .7,
                         observation_map = c(1, 2, 2)), "denominator", class = "recurcast_error")
  # At a subnormal error level, CP lower endpoints round to zero. The declared
  # support envelope must remain conservative instead of claiming a sharp LP endpoint.
  tiny_alpha <- .Machine$double.xmin * .Machine$double.eps
  wide <- rc_bounds(c(.5, .5), c(0, 1), kappa = 2, recording_fraction = .7,
                    counts = c(1, 1), alpha = tiny_alpha)
  expect_identical(wide$outer$method, "support envelope; denominator lower bound is zero")
  expect_equal(c(wide$outer$lower, wide$outer$upper), c(0, 1))
  # A numerator positive outside its conditioning event has no bounded ratio envelope.
  expect_error(rc_bounds(c(.5, .5), c(1, 1), c(0, 1), 2, .7,
                         counts = c(1, 1), alpha = tiny_alpha),
               "bounded conditional support envelope", class = "recurcast_error")
})

test_that("conformal certificates require coherent weight bands and an error budget", {
  sc <- c(.1, .5, .9)
  expect_error(rc_conformal(sc, .7, band_error = .05), "smaller than alpha", class = "recurcast_error")
  expect_error(rc_conformal(sc, .7, weight_status = "certified"), "requires all", class = "recurcast_error")
  expect_error(rc_conformal(sc, .7, weights = c(1, 2)), "dimension", class = "recurcast_error")
  expect_error(rc_conformal(sc, .7, candidate_weights = 0), "positive", class = "recurcast_error")
  expect_error(rc_conformal(sc, NaN), "finite", class = "recurcast_error")
  expect_error(rc_conformal(sc, .7, weight_lower = c(1, 2, 1), weight_upper = c(1, 1, 1),
                            candidate_weight_upper = 1, weight_status = "certified"),
               "interval", class = "recurcast_error")
  zero_lower <- rc_conformal(sc, .7, weight_lower = rep(0, 3), weight_upper = rep(1, 3),
                             candidate_weight_upper = 1, weight_status = "certified")
  expect_equal(zero_lower$candidates$p_value, 1)
  # Large but representable normalizing masses preserve the same probability rank.
  huge <- rc_conformal(sc, .7, weights = rep(1e150, 3), candidate_weights = 1e150)
  expect_equal(huge$candidates$p_value, .5)
  extreme <- rc_conformal(sc, .7, weights = c(1e308, 5e307, 1e308), candidate_weights = 1e308)
  expect_equal(extreme$candidates$p_value, 4 / 7, tolerance = 1e-15)
  rescaled <- rc_conformal(sc, .7, weights = c(2, 1, 2), candidate_weights = 2)
  expect_equal(extreme$candidates$p_value, rescaled$candidates$p_value, tolerance = 1e-15)
  irrelevant_upper <- rc_conformal(sc, .7, weight_lower = rep(1e-100, 3),
    weight_upper = c(1e308, 1e308, 1e-100), candidate_weight_upper = 1e-100,
    weight_status = "certified")
  expect_equal(irrelevant_upper$candidates$p_value, .5, tolerance = 1e-15)
})

test_that("pivot and expansion inputs respect their inferential parameter domains", {
  expect_error(rc_pivot(c(1, 2, 3), .7, c(0, 1)), "positive", class = "recurcast_error")
  expect_error(rc_pivot(c(1, 2, 3), -1, 1), class = "recurcast_error")
  expect_error(rc_pivot(c(1, 2, 3), .7, 1, alpha = 1), class = "recurcast_error")
  expect_error(rc_cornish_fisher(c(0, 1), 100), "inside", class = "recurcast_error")
  expect_error(rc_cornish_fisher(.5, 100, order = 1.5), "0, 1 or 2", class = "recurcast_error")
  expect_error(rc_cornish_fisher(.5, 100, skewness = NA_real_), class = "recurcast_error")
  # The symmetric first-order correction vanishes when the parent skewness is zero.
  expect_equal(rc_cornish_fisher(c(.025, .975), 100, skewness = 0, order = 1),
               qnorm(c(.025, .975)), tolerance = 1e-14)
})

test_that("identification iteration reports invalid contractions and finite precision", {
  expect_error(rc_identify(.8, .5, 2), "observable contraction", class = "recurcast_error")
  expect_error(rc_identify(identity, .5, 2, iterations = 2.5), "integer", class = "recurcast_error")
  expect_error(rc_identify(function(s) -s, .5, 2), "iterated survival", class = "recurcast_error")
  # A supplied operator that does not contract cannot silently produce an infinite baseline.
  expect_error(rc_identify(identity, .5, 1e100, iterations = 4), "numerical range", class = "recurcast_error")
  rounded <- rc_identify(function(s) 1, .5, 2, iterations = 2)
  expect_true(attr(rounded, "rounding_limit"))
  expect_equal(rounded$baseline, c(0, 0))
})

test_that("identity bands propagate correlated forecasts and RNG creation is reversible", {
  forecast <- function(b) c(b[1] + b[2], b[1] - b[2])
  covariance <- matrix(c(.04, .01, .01, .09), 2)
  band <- rc_band(c(1, 2), covariance, forecast, draws = 200, seed = 15, scale = "identity")
  expect_equal(band$estimate, c(3, -1))
  expect_true(all(band$lower < band$estimate & band$upper > band$estimate))
  expect_equal(band$upper - band$estimate, band$estimate - band$lower, tolerance = 1e-14)
  without_seed <- function() {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", saved, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
    if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
    z <- rc_band(0, matrix(.04), plogis, draws = 25, seed = 19)
    expect_length(z$lower, 1)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  }
  without_seed()
})

test_that("nonlinear bands reject invalid covariance, degenerate or failed forecasts", {
  expect_error(rc_band(0, matrix(.1), plogis, draws = 2.5), class = "recurcast_error")
  expect_error(rc_band(c(0, 0), matrix(c(1, .2, .3, 1), 2), function(x) .5),
               "symmetric", class = "recurcast_error")
  expect_error(rc_band(0, matrix(.1), function(x) NA_real_), "point estimate", class = "recurcast_error")
  fails_draws <- function(x) if (identical(as.numeric(x), 0)) .5 else NA_real_
  err <- tryCatch(rc_band(0, matrix(.1), fails_draws, draws = 25), error = identity)
  expect_s3_class(err, "recurcast_error")
  expect_identical(err$code, "forecast_draw_failure")
  expect_error(rc_band(0, matrix(0), plogis, draws = 25), "zero Gaussian draw variance", class = "recurcast_error")
  expect_error(rc_band(0, matrix(.1), function(x) 1, draws = 25), "strictly inside", class = "recurcast_error")
})


test_that("likelihood calculations remain stable at the vanishing-frailty boundary", {
  design <- matrix(1, 2, 1)
  poisson <- rc_loglik(0, 1:2, c(1, 2), design, c(1, 0), phi_fixed = 0)
  near <- rc_loglik(c(0, log(1e-8)), 1:2, c(1, 2), design, c(1, 0))
  expect_equal(near$loglik, poisson$loglik, tolerance = 1e-7)
  expect_equal(near$scores[, 1], poisson$scores[, 1], tolerance = 1e-7)
  expect_true(all(is.finite(near$information)))
  expect_true(all(is.finite(near$scores[, 2])))
  expect_error(rc_loglik(c(0, -1e6), 1:2, c(1, 2), design, c(1, 0)),
               "dispersion outside numerical range", class = "recurcast_error")
  expect_error(rc_loglik(c(0, 0), 1:2, c(1, 2), design, c(1, 0), phi_fixed = .5),
               "dimensions", class = "recurcast_error")
  expect_error(rc_loglik(0, 1:2, c(1, 2), c(1, 1), c(1, 0), phi_fixed = 0), class = "recurcast_error")
  expect_error(rc_loglik(0, 1:2, 1, design, c(1, 0), phi_fixed = 0), class = "recurcast_error")
  expect_error(rc_loglik(0, 1:2, c(1, 2), design, c(1, 0), phi_fixed = 0, information = NA), class = "recurcast_error")
})

test_that("coarsening can aggregate scores without requesting information matrices", {
  coarsen <- getFromNamespace("rc_coarsen", "recurcast")
  weights <- c(.25, .75)
  z <- coarsen(log(c(.2, .4)), matrix(c(1, 3), 2, 1), owner = c(1, 1),
                log_weights = log(weights))
  expect_equal(z$loglik, log(.25 * .2 + .75 * .4), tolerance = 1e-14)
  expect_equal(as.vector(z$scores), (.25 * .2 * 1 + .75 * .4 * 3) / (.25 * .2 + .75 * .4), tolerance = 1e-14)
  expect_null(z$information)
  expect_error(coarsen(log(c(.2, .4)), matrix(c(1, 3), 2, 1), owner = c(1, 3),
                       log_weights = log(weights)), "owner indices", class = "recurcast_error")
  expect_error(coarsen(log(c(.2, .4)), matrix(c(1, 3), 2, 1), owner = c(1, 1),
                       log_weights = 0), "one log integration weight", class = "recurcast_error")
  expect_error(coarsen(log(c(.2, .4)), matrix(c(1, 3), 2, 1), information = array(1, c(2, 2, 2)),
                       owner = c(1, 1), log_weights = log(weights)), "dimensions", class = "recurcast_error")
})

test_that("ordered propagation rejects incompatible probability and reward dimensions", {
  q <- matrix(c(-.2, 0, .2, 0), 2, 2)
  expect_error(rc_propagate(c(NA_real_, 0), q, 1), "finite numeric", class = "recurcast_error")
  expect_error(rc_propagate(c(1, 0), diag(3), 1), "square matrix", class = "recurcast_error")
  expect_error(rc_propagate(c(1, 0), list(q, q), c(1, 1), reward = list(c(1, 0))),
               "one reward matrix", class = "recurcast_error")
  expect_error(rc_propagate(c(1, 0), q, 1, reward = matrix(1, 3, 1)),
               "share dimensions", class = "recurcast_error")
  # Each segment must expose the same reward coordinates, even if values differ.
  expect_error(rc_propagate(c(1, 0), list(q, q), c(1, 1), reward = list(matrix(1, 2, 1), matrix(1, 2, 2))),
               "share dimensions", class = "recurcast_error")
})
