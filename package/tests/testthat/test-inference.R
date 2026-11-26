fixture <- function(name) read.csv(test_path("fixtures", name), check.names = FALSE)

test_that("cross-fitted targeting agrees with the preserved Python oracle", {
  a <- fixture("target_subjects.csv"); f <- as.matrix(fixture("target_features.csv"))
  expected <- fixture("target_expected.csv")
  z <- cs_target(a$path, a$recorded, a$strata, f, a$fold)
  expect_true(z$convergence$valid)
  expect_equal(unname(z$estimate), expected$estimate, tolerance = 1e-10)
  expect_equal(unname(z$se), expected$se, tolerance = 1e-10)
  expect_equal(unname(z$EIF), unname(as.matrix(fixture("target_eif.csv"))), tolerance = 1e-9)
  expect_equal(unname(z$probabilities), unname(as.matrix(fixture("target_laws.csv"))), tolerance = 1e-10)
  expect_lt(max(abs(z$EIF_mean)), 1e-6 / sqrt(nrow(a)))
  expect_equal(rowSums(z$group_probabilities), rep(1, nrow(z$group_probabilities)), tolerance = 1e-13)
  expect_equal(as.vector(crossprod(z$mixture, f)), unname(z$estimate), tolerance = 1e-13)
  a$path[a$recorded == 0] <- 1e9
  hidden <- cs_target(a$path, a$recorded, a$strata, f, a$fold)
  expect_identical(z$estimate, hidden$estimate)
  expect_error(cs_target(a$path, a$recorded, a$strata, f, a$fold,
                         recording = rep(0, 10)), class = "coherentSurv_error")
  expect_error(cs_target(a$path, rep(0, nrow(a)), a$strata, f, a$fold), class = "coherentSurv_error")
  expect_error(cs_target(a$path, a$recorded, a$strata, f, control = list(folds = 0)), class = "coherentSurv_error")
})

test_that("Charnes-Cooper bounds agree with HiGHS and independent vertices", {
  tab <- fixture("bounds_expected.csv")
  for (i in seq_len(nrow(tab))) {
    j <- tab$J[i]; p <- j:1 / sum(j:1)
    z <- cs_bounds(p, as.numeric(seq_len(j) == j), as.numeric(seq_len(j) > 1),
                   tab$kappa[i], .7)
    expect_equal(z$lower, tab$lower[i], tolerance = 1e-12)
    expect_equal(z$upper, tab$upper[i], tolerance = 1e-12)
    expect_lt(z$plugin$primal_dual_gap, 1e-12)
    expect_lt(z$plugin$primal_residual, 1e-12)
    expect_equal(sum(z$plugin$lower_law), 1, tolerance = 1e-12)
    expect_equal(sum(z$plugin$upper_law), 1, tolerance = 1e-12)
  }
  z <- cs_bounds(c(.5, .5), c(0, 1, 0), kappa = 1, recording_fraction = .7,
                 observation_map = c(1, 2, 2))
  expect_equal(c(z$lower, z$upper), c(0, .5), tolerance = 1e-12)
  expect_error(cs_bounds(c(.5, .5), c(0, 1), kappa = .5, recording_fraction = .7), class = "coherentSurv_error")
  expect_error(cs_bounds(c(.5, .5), c(0, 1), denominator = c(0, 0), kappa = 2, recording_fraction = .7), class = "coherentSurv_error")
})

test_that("confidence projection agrees with preserved simultaneous CP calculations", {
  tab <- fixture("projection_expected.csv")
  for (i in seq_len(nrow(tab))) {
    counts <- as.numeric(strsplit(tab$counts[i], ";", fixed = TRUE)[[1]])
    j <- length(counts)
    z <- cs_bounds(counts / sum(counts), as.numeric(seq_len(j) == j), as.numeric(seq_len(j) > 1),
                   tab$kappa[i], .7, counts = counts)
    expect_equal(z$outer$lower, tab$lower[i], tolerance = 2e-11)
    expect_equal(z$outer$upper, tab$upper[i], tolerance = 2e-11)
    expect_lte(z$outer$lower, z$lower + 1e-12)
    expect_gte(z$outer$upper, z$upper - 1e-12)
  }
  zero <- cs_bounds(c(1, 0), c(0, 1), kappa = 2, recording_fraction = .7, counts = c(10, 0))
  expect_equal(zero$lower, 0)
  expect_gt(zero$outer$upper, 0)
})

test_that("weighted ranks include candidate mass, ties and conservative weight bands", {
  known <- cs_conformal(c(.1, .5, .9), c(.5, .7), alpha = .1)
  expect_equal(known$candidates$p_value, c(.75, .5))
  scaled <- cs_conformal(c(.1, .5, .9), .7, rep(3, 3), 3)
  expect_equal(scaled$candidates$p_value, .5)
  certified <- cs_conformal(c(.1, .5, .9), .7, weight_lower = rep(1, 3),
                           weight_upper = rep(2, 3), candidate_weight_upper = 2,
                           band_error = .01, weight_status = "certified")
  expect_equal(certified$candidates$p_value, 2 / 3)
  expect_equal(certified$rank_alpha, .04)
  inf <- cs_conformal(c(.1, .9), .7, weight_lower = c(1, 1),
                      weight_upper = c(2, Inf), candidate_weight_upper = 2,
                      weight_status = "certified")
  expect_equal(inf$candidates$p_value, 1)
  estimated <- cs_conformal(c(.1, .9), .7, weight_status = "estimated")
  expect_match(estimated$coverage_guarantee, "No exact")
  expect_error(cs_conformal(c(.1, .9), .7, band_error = .01), class = "coherentSurv_error")
  expect_error(cs_conformal(c(.1, .9), .7, weights = c(-1, 2)), class = "coherentSurv_error")
})

test_that("the known-dispersion pivot includes its analytic exponential boundary", {
  z <- cs_pivot(c(1, 2, 3, 4), 0, c(.5, 1, 2))
  expect_equal(unname(z$rate_interval), qgamma(c(.025, .975), 4) / 10, tolerance = 1e-12)
  expect_true(all(diff(z$risk$lower) > 0))
  z2 <- cs_pivot(c(1, 2, 3, 4), .7, c(.5, 1, 2))
  expect_equal(vapply(z2$rate_interval, function(lam) sum(log1p(.7 * lam * c(1, 2, 3, 4))) / .7, numeric(1)),
               setNames(qgamma(c(.025, .975), 4), c("lower", "upper")), tolerance = 1e-11)
  expect_error(cs_pivot(c(0, 2), .7, 1), class = "coherentSurv_error")
  expect_equal(cs_cornish_fisher(.025, 20), -1.7395760556022086, tolerance = 1e-12)
  expect_equal(cs_cornish_fisher(c(.025, .975), 80, order = 0), qnorm(c(.025, .975)), tolerance = 1e-14)
})

test_that("population identification iteration recovers a known gamma baseline", {
  phi <- .7; baseline <- .8
  initial <- (1 + phi * baseline)^(-1 / phi)
  contraction <- function(s) (1 + (s^(-phi) - 1) / 2)^(-1 / phi)
  z <- cs_identify(contraction, initial, 2, iterations = 15)
  expect_lt(abs(tail(z$baseline, 1) - baseline), 2e-5)
  expect_true(all(diff(abs(z$baseline - baseline)) < 0))
})

test_that("Gaussian nonlinear bands preserve RNG state and fail visibly at boundaries", {
  set.seed(51); before <- .Random.seed
  z <- cs_band(0, matrix(.04), function(b) plogis(b + c(-1, 0, 1)), draws = 200, seed = 19)
  expect_identical(.Random.seed, before)
  expect_true(all(z$lower > 0 & z$upper < 1 & z$lower < z$upper))
  one <- cs_band(0, matrix(.04), function(b) plogis(b), draws = 200, seed = 19)
  expect_length(one$lower, 1)
  expect_true(one$lower < .5 && one$upper > .5)
  expect_match(z$scope, "neither an exact")
  expect_error(cs_band(0, matrix(-1), plogis), class = "coherentSurv_error")
  expect_error(cs_band(0, matrix(.04), function(b) 0), class = "coherentSurv_error")
})
