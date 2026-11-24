test_that("CSR propagation preserves ordered probability and reward calculations", {
  q1 <- rc_generator_csr(3, c(1, 1, 2), c(2, 3, 3), c(.3, .2, .4))
  q2 <- rc_generator_csr(3, c(1, 1, 2), c(2, 3, 3), c(.1, .7, .2))
  d1 <- matrix(c(-.5, .3, .2, 0, -.4, .4, 0, 0, 0), 3, byrow = TRUE)
  d2 <- matrix(c(-.8, .1, .7, 0, -.2, .2, 0, 0, 0), 3, byrow = TRUE)
  reward <- list(cbind(c(1, 1, 0), c(.3, -.1, 0)),
                 cbind(c(1, 1, 0), c(.1, -.2, 0)))
  s <- rc_propagate_sparse(c(1, 0, 0), list(q1, q2), c(.7, 35.3), reward)
  d <- rc_propagate(c(1, 0, 0), list(d1, d2), c(.7, 35.3), reward)
  expect_equal(s$probability, d$probability, tolerance = 1e-11)
  expect_equal(s$reward, d$reward, tolerance = 1e-10)
  expect_equal(sum(s$probability), 1, tolerance = 1e-14)
  expect_lte(s$error_bound, 1e-12)
  expect_gt(s$chunks, 2)
  qdup <- rc_generator_csr(2, c(1, 1), c(2, 2), c(.2, .3))
  expect_equal(rc_propagate_sparse(c(1, 0), qdup, 2)$probability[1], exp(-1), tolerance = 1e-12)
  zero <- rc_generator_csr(1)
  expect_equal(rc_propagate_sparse(1, zero, 3, 2)$reward, 6)
})

test_that("CSR has explicit malformed-input and computational-budget errors", {
  expect_error(rc_generator_csr(3, 1, 1, .1), class = "recurcast_error")
  expect_error(rc_generator_csr(3, 1, 4, .1), class = "recurcast_error")
  q <- rc_generator_csr(3, c(1, 2), c(2, 3), c(.1, .2))
  q$offsets <- c(0L, 100L, 2L, 2L)
  expect_error(rc_propagate_sparse(c(1, 0, 0), q, 1), class = "recurcast_error")
  q <- rc_generator_csr(2, 1, 2, 100)
  expect_error(rc_propagate_sparse(c(1, 0), q, 1, max_chunks = 1), "max_chunks")
})

test_that("two-atom fractional bounds match exact hand enumeration", {
  z <- rc_bounds_sparse(c(.5, .5), c(0, 1), c(1, 2), 2, .5)
  # Recorded a=(1/4,1/4); b_2 ranges from 1/8 to 3/8.
  # Thus q_2 ranges from 3/8 to 5/8, and q_2/(1+q_2) is monotone.
  expect_equal(z$lower, 3/11, tolerance = 1e-12)
  expect_equal(z$upper, 5/13, tolerance = 1e-12)
  expect_equal(z$plugin$lower_law, c(5/8, 3/8), tolerance = 1e-12)
  expect_equal(z$plugin$upper_law, c(3/8, 5/8), tolerance = 1e-12)
  expect_lte(z$plugin$gap, 1e-10)
  expect_lte(z$plugin$residual, 1e-12)
})

test_that("sparse bounds reproduce dense coarsened projection and MAR boundary", {
  p <- c(.3, .4, .3); map <- c(1, 1, 2, 2, 3, 3)
  a <- c(.1, .7, .4, .2, .8, .6); d <- c(1, .9, .7, 1, .9, .8)
  for (kappa in c(1, 2, 4)) {
    sparse <- rc_bounds_sparse(p, a, d, kappa, .6, counts = c(30, 40, 30), observation_map = map)
    dense <- rc_bounds(p, a, d, kappa, .6, counts = c(30, 40, 30), observation_map = map)
    expect_equal(c(sparse$lower, sparse$upper), c(dense$lower, dense$upper), tolerance = 1e-9)
    expect_equal(c(sparse$outer$lower, sparse$outer$upper), c(dense$outer$lower, dense$outer$upper), tolerance = 1e-9)
  }
  mar <- rc_bounds_sparse(c(.3, .7), c(0, 1), kappa = 1, recording_fraction = .6)
  expect_equal(c(mar$lower, mar$upper), c(.7, .7), tolerance = 1e-12)
  ties <- rc_bounds_sparse(c(.5, .5), rep(.5, 4), rep(1, 4), 3, .6, observation_map = c(1, 1, 2, 2))
  expect_equal(c(ties$lower, ties$upper), c(.5, .5), tolerance = 1e-12)
})

test_that("fractional primal and dual certificates can be independently reconstructed", {
  p <- c(.3, .4, .3); map <- c(1, 1, 2, 2, 3, 3)
  num <- c(.1, .7, .4, .2, .8, .6); den <- c(1, .9, .7, 1, .9, .8)
  z <- rc_bounds_sparse(p, num, den, 3, .6, counts = c(30, 40, 30), observation_map = map)
  l <- (.4/.6)/3; u <- (.4/.6)*3; Y <- (.4-l*.6)/(u-l)
  for (which in c("lower_certificate", "upper_certificate")) {
    e <- z$outer[[which]]; c <- e$objective_sign*num-e$linearization_point*den
    best <- vapply(seq_along(p), function(h) max(c[map == h]), numeric(1))
    dh <- (1+l)*best + pmax(0, (u-l)*best-e$dual_high_mass) - e$dual_mass
    ci <- z$recorded_confidence_region
    dual <- .6*e$dual_mass + Y*e$dual_high_mass + sum(pmax(dh, 0)*ci[, 2] + pmin(dh, 0)*ci[, 1])
    expect_equal(dual, e$linear_dual, tolerance = 1e-12)
    expect_equal(sum(c*e$law), e$linear_primal, tolerance = 1e-12)
    expect_gte(dual + 1e-12, e$linear_primal)
    expect_equal(sum(e$recorded), .6, tolerance = 1e-12)
    expect_equal(sum(e$unrecorded), .4, tolerance = 1e-12)
    expect_true(all(e$unrecorded >= l*e$recorded-1e-12))
    expect_true(all(e$unrecorded <= u*e$recorded+1e-12))
    expect_lte(e$certificate_lower, e$value+1e-12)
    expect_gte(e$certificate_upper, e$value-1e-12)
  }
})

test_that("50,000 atoms and 50,000 states use the sparse route", {
  J <- 50000L; M <- 1000L; k <- seq_len(J)-1L
  map <- k %% M+1L; p <- rep(1/M, M)
  den <- .2+.8*((k*104729) %% 50023)/50023
  num <- den*(.1+.8*((k*7919) %% 50021)/50021)
  z <- rc_bounds_sparse(p, num, den, 3, .6, observation_map = map)
  expect_lte(z$plugin$gap, 1e-8)
  expect_lte(z$plugin$residual, 1e-10)
  expect_gt(z$lower, 0); expect_lt(z$upper, 1)
  q <- rc_generator_csr(J, seq_len(J-1L), seq_len(J-1L)+1L, rep(.1, J-1L))
  prob <- numeric(J); prob[1] <- 1
  out <- rc_propagate_sparse(prob, q, 10, rep(1, J))
  expect_equal(out$probability[1], exp(-1), tolerance = 1e-12)
  expect_equal(out$reward, 10, tolerance = 1e-10)
  expect_lt(out$workspace_bytes, 10e6)
  expect_lt(as.numeric(object.size(q)), 2e6)
})
