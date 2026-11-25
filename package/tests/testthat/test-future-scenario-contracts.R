test_that("longitudinal submodels enforce measured-data and RNG contracts", {
  bad <- data.frame(id = rep(1:3, each = 2), time = rep(0:1, 3), x = 1:6)
  expect_error(rc_covariate_model(bad, variables = character()), "covariate_model", class = "coherent4cast_error")
  bad$x[1] <- Inf
  expect_error(rc_covariate_model(bad, variables = "x"), "finite numeric", class = "coherent4cast_error")
  local_seed <- getFromNamespace(".rc_with_seed", "coherent4cast")
  expect_error(local_seed(-1, 0), "seed", class = "coherent4cast_error")
  set.seed(9191); before <- .Random.seed
  a <- local_seed(123, stats::runif(4))
  expect_identical(.Random.seed, before)
  expect_equal(a, local_seed(123, stats::runif(4)))
  # Test the documented isolation contract when the caller has no RNG state.
  rm(".Random.seed", envir = .GlobalEnv)
  b <- local_seed(123, stats::runif(4))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  assign(".Random.seed", before, .GlobalEnv)
  expect_equal(a, b)
  draws <- getFromNamespace(".rc_covariate_draws", "coherent4cast")
  past <- data.frame(time = 0, x = 1)
  expect_error(draws(NULL, past, 1, "time", 2, 1), "covariate_model", class = "coherent4cast_error")
  model <- structure(list(variables = "x", center = 0, scale = 1,
      random_slope = FALSE, models = list(x = list(beta = c(0, 0), D = matrix(1), sigma2 = 1))),
      class = "rc_covariate_model")
  expect_error(draws(model, past["time"], 1, "time", 2, 1), "Past history", class = "coherent4cast_error")
  # Conditional N(1/2,1/2) random intercept, with no slope: both future times agree.
  sampled <- draws(model, past, c(1, 2), "time", 2000, 11)
  values <- vapply(sampled, function(z) z[1, 1], numeric(1))
  expect_equal(mean(values), .5, tolerance = .05)
  expect_equal(stats::var(values), .5, tolerance = .05)
  expect_true(all(vapply(sampled, function(z) z[1, 1] == z[2, 1], logical(1))))
  model$models$x$D <- matrix(-.1)
  expect_error(draws(model, past, 1, "time", 2, 1), "positive semidefinite", class = "coherent4cast_error")
})

test_that("future schedules and stochastic scenarios are distinct and validated", {
  data <- rc_simulate(30, horizon = 2, rates = c(.4, .2), beta = c(.1, .1),
                      history = c(.1, .1), seed = 726)
  fit <- suppressWarnings(coherent4cast(~x, data, time_unit = "years",
      control = coherent4cast_control(phi_fixed = .4), frailty_diagnostics = FALSE))
  past <- data.frame(id = "A", start = 0, stop = 1, event_type = 0, terminal = 0, x = 0)
  p <- function(...) predict(fit, past, landmark = 1, horizon = .2, quadrature = 2, ...)
  schedule <- data.frame(id = "A", time = 1.1, x = 1)
  expect_error(p(newdata_future = schedule), "requires", class = "coherent4cast_error")
  expect_error(p(covariate_evolution = "locf", newdata_future = schedule["time"]), "Future schedules", class = "coherent4cast_error")
  expect_error(p(covariate_evolution = "locf", newdata_future = rbind(schedule, schedule)), "duplicate future", class = "coherent4cast_error")
  frozen <- p()
  expect_equal(p(covariate_evolution = "locf")$estimate, frozen$estimate)
  expect_error(p(covariate_evolution = "joint_submodel", newdata_future = schedule), "separate forecast", class = "coherent4cast_error")
  expect_error(p(covariate_evolution = "joint_submodel"), "fitted rc_covariate_model", class = "coherent4cast_error")
  model <- structure(list(variables = "z"), class = "rc_covariate_model")
  expect_error(p(covariate_evolution = "joint_submodel", covariate_model = model), "event-model formula", class = "coherent4cast_error")
  model$variables <- "x"
  expect_error(p(covariate_evolution = "joint_submodel", covariate_model = model, covariate_draws = 1), "at least two", class = "coherent4cast_error")
  expect_error(p(covariate_evolution = "joint_submodel", covariate_model = model), "covariate_history", class = "coherent4cast_error")
  measure <- data.frame(id = c("A", "A"), time = c(0, 0), x = c(0, 0))
  expect_error(p(covariate_evolution = "joint_submodel", covariate_model = model, covariate_history = measure), "unique measurement", class = "coherent4cast_error")
})
