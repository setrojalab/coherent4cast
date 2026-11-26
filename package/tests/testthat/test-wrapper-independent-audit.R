p3_competing_data <- function() data.frame(id = 1:12,
  time = c(.2, .3, .4, .5, .6, .8, 1, 1.2, 1.3, 1.5, 1.7, 2),
  status = c(1, 2, 0, 1, 0, 2, 1, 0, 2, 0, 1, 0), x = rep(0:1, 6))

test_that("weighted competing likelihood and joint sandwich match direct subject formulas", {
  d <- p3_competing_data(); w <- rep(c(.5, 1.5, 2), 4)
  f <- cs_competing(~1, d, unit = "id", time_unit = "years", weights = w, frailty_diagnostics = FALSE)
  events <- cbind(as.integer(d$status == 1), as.integer(d$status == 2))
  counts <- colSums(events * w); rates <- counts / sum(w * d$time)
  loglik <- sum(sweep(events, 2, log(rates), "*") * w) - sum(w * d$time) * sum(rates)
  score <- events - outer(d$time, rates)
  H <- diag(counts); inv <- solve(H)
  oracle <- inv %*% crossprod(score * w) %*% inv
  expect_equal(unname(exp(coef(f))), rates, tolerance = 1e-7)
  expect_equal(as.numeric(logLik(f)), loglik, tolerance = 1e-8)
  expect_equal(unname(vcov(f)), oracle, tolerance = 1e-7)
  expect_true(abs(oracle[1, 2]) > 1e-5)
})

test_that("valid changes strictly after a landmark cannot change competing forecasts", {
  f <- cs_competing(~x, p3_competing_data(), unit = "id", time_unit = "years", frailty_diagnostics = FALSE)
  first <- data.frame(id = 100, left = 0, right = 2, status = 0, x = 0)
  train <- transform(p3_competing_data(), left = 0, right = time)
  fs <- cs_competing(~x, train, unit = "id", start = "left", time = "right", time_unit = "years", frailty_diagnostics = FALSE)
  past <- first; past$right <- .5
  future <- first; future$left <- .5; future$x <- 8; future$status <- 2
  a <- predict(fs, rbind(past, future), landmark = .25, horizon = .5, type = "all")
  future$x <- -8; future$status <- 1; future$right <- 3
  b <- predict(fs, rbind(past, future), landmark = .25, horizon = .5, type = "all")
  expect_equal(a$estimate, b$estimate, tolerance = 1e-13)
  profile <- predict(f, data.frame(id = 100, x = 0), landmark = .25, horizon = .5, type = "all")
  expect_equal(a$estimate, profile$estimate, tolerance = 1e-8)
  expect_match(attr(profile, "profile_assumption"), "event-free")
})

test_that("factor-coded competing validation uses event-time weights for competing controls", {
  d <- p3_competing_data(); d$status <- factor(ifelse(d$status == 0, "censored", paste0("c", d$status)))
  f <- cs_competing(~1, d, censor_code = "censored", time_unit = "years", frailty_diagnostics = FALSE)
  status <- factor(c("c1", "c2", "censored", "censored", "c1"), levels = levels(d$status))
  time <- c(.2, .4, .6, 2, 2); p <- c(.8, .7, .4, .1, .2)
  # All failures use G(T-); the early competing cause remains a known noncase.
  g <- list(event = c(.8, .5, 1, 1, 1), horizon = rep(.25, 5))
  v <- validate(f, time, status, horizon = 1, cause = "c1", predicted = p,
    censor_survival = g, threshold = .5)
  expected_weights <- c(1.25, 2, 0, 4, 4)
  y <- c(1, 0, 0, 0, 0)
  expect_equal(v$weights, expected_weights)
  expect_equal(v$metrics$estimate[v$metrics$metric == "Brier"], sum(expected_weights * (y-p)^2)/5, tolerance = 1e-14)
  expect_equal(v$usable, 4)
  expect_equal(v$controls, 3)
  expect_equal(v$competing_events_by_horizon, 1)
  sub <- validate(f, time, status, horizon = 1, cause = "c1", predicted = p,
    censor_survival = g, strata = c("A", "A", "B", "B", "B"))
  expect_equal(sub$A$weights, c(1.25, 2))
  expect_match(sub$A$endpoint, "observed noncases")
})

test_that("terminal eligibility and malformed survival formulas fail explicitly", {
  d <- p3_competing_data(); f <- cs_competing(~1, d, unit = "id", time_unit = "years", frailty_diagnostics = FALSE)
  ended <- data.frame(id = c(100, 101), time = c(.4, .3), status = c(1, 0), x = 0)
  expect_error(predict(f, ended, landmark = .5, horizon = 1), "No unit", class = "coherentSurv_error")
  expect_error(cs_surv(survival::Surv(time, status == 1) ~ x, d, time_unit = "years"), "column names", class = "coherentSurv_error")
  binary <- transform(d, status = as.integer(status > 0))
  binary$entry <- 0
  fit <- cs_surv(survival::Surv(entry, time, status) ~ x, binary, unit = "id", time_unit = "years", frailty_diagnostics = FALSE)
  expect_identical(fit$wrapper$start, "entry")
  expect_identical(fit$frailty_requested, "none")
  expect_equal(fit$phi, 0)
})
