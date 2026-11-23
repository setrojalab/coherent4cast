test_that("shared-gamma scores and information reproduce the V5 Python oracle", {
  fixture <- read.csv(test_path("fixtures", "gamma_design.csv"))
  design <- as.matrix(fixture[c("x1", "x2")])
  for (label in c("estimated", "fixed", "poisson")) {
    phi <- switch(label, estimated = NULL, fixed = 0.6, poisson = 0)
    theta <- if (is.null(phi)) c(-0.7, 0.25, log(0.6)) else c(-0.7, 0.25)
    fit <- rc_loglik(theta, fixture$unit, fixture$exposure, design, fixture$event,
                     phi_fixed = phi)
    expected <- function(kind) as.matrix(read.csv(test_path("fixtures", paste0("gamma_", label, "_", kind, ".csv"))))
    expect_equal(unname(fit$scores), unname(expected("scores")), tolerance = 1e-12)
    expect_equal(unname(fit$information), unname(expected("information")), tolerance = 1e-12)
    expect_equal(unname(fit$subject_loglik), as.numeric(expected("loglik")), tolerance = 1e-12)
    expect_equal(fit$information, apply(fit$information_individual, c(1, 2), sum), tolerance = 1e-12)
    simple <- rc_loglik(theta, fixture$unit, fixture$exposure, design, fixture$event,
                        phi_fixed = phi, information = FALSE)
    expect_equal(simple$scores, fit$scores)
    expect_null(simple$information)
  }
})

test_that("legal-path mixture uses quadrature weights and Louis information", {
  d <- read.csv(test_path("fixtures", "gamma_design.csv"))
  base <- function(theta) rc_loglik(theta, d$unit, d$exposure,
                                   as.matrix(d[c("x1", "x2")]), d$event)
  theta <- c(-0.7, 0.25, log(0.6))
  mix <- function(theta) recurcast:::rc_coarsen(base(theta), owner = c(1, 1, 2, 2),
                                               log_weights = log(c(.2, .3, .5, .7)))
  fit <- mix(theta)
  expected <- log(sum(exp(base(theta)$subject_loglik[1:2]) * c(.2, .3))) +
    log(sum(exp(base(theta)$subject_loglik[3:4]) * c(.5, .7)))
  expect_equal(fit$loglik, expected, tolerance = 1e-12)
  expect_equal(rowsum(fit$posterior, c(1, 1, 2, 2)), matrix(1, 2, 1, dimnames = list(c("1", "2"), NULL)))
  h <- 2e-4
  numerical <- sapply(seq_along(theta), function(j) {
    score <- function(k) {v <- theta; v[j] <- v[j] + k * h; colSums(mix(v)$scores)}
    -(score(-2) - 8 * score(-1) + 8 * score(1) - score(2)) / (12 * h)
  })
  expect_equal(fit$information, numerical, tolerance = 1e-9)
})

test_that("ordered propagation and integrated rewards match dense matrix exponentials", {
  q1 <- matrix(c(-.5, .3, .2, 0, -.4, .4, 0, 0, 0), 3, byrow = TRUE)
  q2 <- matrix(c(-.8, .1, .7, 0, -.2, .2, 0, 0, 0), 3, byrow = TRUE)
  r1 <- cbind(c(.3, 0, 0), c(1, 1, 0))
  r2 <- cbind(c(.1, 0, 0), c(1, 1, 0))
  one <- rc_propagate(c(1, 0, 0), q1, 2, r1)
  two <- rc_propagate(c(1, 0, 0), list(q1, q2), c(.7, 1.3), list(r1, r2))
  expected <- function(s) as.numeric(read.csv(test_path("fixtures", paste0("propagation_", s, ".csv"))))
  expect_equal(c(one$probability, one$reward), expected("single"), tolerance = 1e-12)
  expect_equal(c(two$probability, two$reward), expected("ordered"), tolerance = 1e-12)
  expect_lt(one$error_bound, 1e-12)
  expect_equal(sum(two$probability), 1, tolerance = 1e-14)
  extreme <- rc_propagate(c(1, 0), matrix(c(-exp(30), exp(30), 0, 0), 2, byrow = TRUE), 1, c(1, 0))
  expect_true(all(is.finite(unlist(extreme))))
  expect_equal(extreme$probability, c(0, 1), tolerance = 1e-12)
  zero <- rc_propagate(1, matrix(0, 1, 1), 3, 2)
  expect_equal(zero$reward, 6)
})

test_that("kernel guards fail at the R boundary instead of producing NaN", {
  expect_error(rc_loglik(0, 1, 1, matrix(1), 1, phi_fixed = -1), "phi_fixed")
  expect_error(rc_loglik(0, 2, 1, matrix(1), 1, phi_fixed = 0), "unit")
  expect_error(rc_loglik(0, 1, -1, matrix(1), 1, phi_fixed = 0), "exposure")
  expect_error(rc_loglik(0, 1, 1, matrix(NA_real_), 1, phi_fixed = 0), "design")
  expect_error(rc_propagate(c(1, 0), matrix(c(-.5, .4, 0, 0), 2, byrow = TRUE), 1), "sum to zero")
  expect_error(rc_propagate(c(1, 0), diag(2), 1), "generator")
  expect_error(rc_propagate(1, matrix(0), -1), "duration")
  expect_error(rc_propagate(1, matrix(0), 1, max_terms = 0), "term")
})
