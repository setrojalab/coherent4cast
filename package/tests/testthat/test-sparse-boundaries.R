test_that("sparse semigroup identity and absorbing rewards hold at zero time", {
  q <- rc_generator_csr(2, 1, 2, .7)
  initial <- c(.4, .6)
  zero <- rc_propagate_sparse(initial, q, 0, cbind(c(1, 0), c(-2, 3)))
  expect_equal(zero$probability, initial, tolerance = 1e-15)
  expect_equal(zero$reward, c(0, 0), tolerance = 1e-15)
  absorbed <- rc_propagate_sparse(c(0, 1), q, 2, cbind(c(1, 0), c(-2, 3)))
  expect_equal(absorbed$probability, c(0, 1), tolerance = 1e-15)
  expect_equal(absorbed$reward, c(0, 6), tolerance = 1e-12)
})

test_that("CSR input contracts prevent malformed generators and silent truncation", {
  expect_error(rc_generator_csr(0), class = "coherent4cast_error")
  expect_error(rc_generator_csr(2, 1, 2, NA_real_), class = "coherent4cast_error")
  expect_error(rc_generator_csr(2, 1.5, 2, .3), class = "coherent4cast_error")
  q <- rc_generator_csr(2, 1, 2, .7)
  expect_error(rc_propagate_sparse(c(NA_real_, 0), q, 1), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), list(q, q), 1), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), q, -1), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0, 0), q, 1), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), q, 1, matrix(1, 3, 1)), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), list(q, q), c(1, 1),
                                  list(matrix(1, 2, 1), matrix(1, 2, 2))),
               class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), q, 1, tolerance = 0), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), q, 1, max_chunks = 1.5), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(.6, .6), q, 1), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(-.1, 1.1), q, 1), class = "coherent4cast_error")
  expect_error(rc_propagate_sparse(c(1, 0), q, 5, max_terms = 2), class = "coherent4cast_error")
  malformed <- q; malformed$columns[1] <- 2L
  expect_error(rc_propagate_sparse(c(1, 0), malformed, 1), class = "coherent4cast_error")
  malformed <- q; malformed$rates[1] <- -.7
  expect_error(rc_propagate_sparse(c(1, 0), malformed, 1), class = "coherent4cast_error")
})

test_that("sparse confidence projection preserves a conservative support envelope", {
  # Subnormal alpha makes binomial lower endpoints round to zero. The output
  # must use the support envelope, not label an undefined ratio as a sharp LP.
  alpha <- .Machine$double.xmin * .Machine$double.eps
  result <- rc_bounds_sparse(c(.5, .5), c(0, 1), kappa = 2,
                             recording_fraction = .7, counts = c(1, 1), alpha = alpha)
  expect_equal(c(result$outer$lower, result$outer$upper), c(0, 1))
  expect_identical(result$outer$method, "support envelope; denominator lower bound is zero")
  expect_error(rc_bounds_sparse(c(.5, .5), c(1, 1), c(0, 1), 2, .7,
                               counts = c(1, 1), alpha = alpha),
               "bounded conditional support envelope", class = "coherent4cast_error")
  expect_error(rc_bounds_sparse(c(.5, .5), c(0, 1, 0), c(0, 1, 0), 2, .7,
                               observation_map = c(1, 2, 2)),
               "denominator", class = "coherent4cast_error")
})

test_that("sparse projection rejects mismatched sampling counts and fibre maps", {
  expect_error(rc_bounds_sparse(c(.5, .5), c(0, 1), kappa = 2,
                               recording_fraction = .7, observation_map = c(1, 3)),
               class = "coherent4cast_error")
  expect_error(rc_bounds_sparse(c(.5, .5), c(0, 1), kappa = 2,
                               recording_fraction = .7, counts = c(2, 3)),
               "proportions equal p", class = "coherent4cast_error")
  expect_error(rc_bounds_sparse(c(.5, .5), c(0, 1), kappa = 2,
                               recording_fraction = .7, counts = c(.5, .5)),
               class = "coherent4cast_error")
  expect_error(rc_bounds_sparse(c(.5, .5), c(0, 1), kappa = 2,
                               recording_fraction = .7, max_iterations = 1.5),
               "integer", class = "coherent4cast_error")
})
