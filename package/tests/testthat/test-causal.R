causal_demo_data <- function(n = 400L) {
  # Deterministic quantiles give an independent analytic exponential comparison.
  x <- rep(0:1, each = n / 2)
  w <- rep(rep(0:1, each = n / 4), 2)
  u <- rep((seq_len(n / 4) - .5) / (n / 4), 4)
  rate <- (.55 + .25 * x) * exp(-.4 * w)
  tt <- -log(u) / rate
  data.frame(id = seq_len(n), start = 0, stop = pmin(tt, 1),
    event_type = "none", terminal = as.integer(tt <= 1), censored = 0,
    x = x, treatment = w)
}

test_that("S12 martingale terms agree with independent exponential algebra", {
  n <- 8L; w <- rep(0:1, each = 4); tt <- rep(c(.1, .4, .7, 1), 2)
  death <- rep(c(1, 0, 1, 0), 2); censor <- .3; lambda <- .9; e <- .4
  rows <- cbind(0:7, 0, 0, 0, tt, ifelse(death == 1, 0, 1))
  tau <- 1; grid <- 256L
  z <- coherent4cast:::rc_cpp_causal_target(rep(lambda, grid), censor, e,
     matrix(1L, 1, 1), rep(0L, n), rows[w == 1, , drop = FALSE], grid, tau, 0L, 1e-10)
  A <- function(rate, t) expm1(rate * t) / rate
  risk <- -expm1(-lambda * tau); rmst <- risk / lambda
  B <- exp(-lambda * (tau - tt))
  risk_res <- (death * exp(censor * tt) * B - lambda * exp(-lambda * tau) * A(lambda + censor, tt)) / e
  rmst_res <- (-death * exp(censor * tt) * (1 - B) / lambda + A(censor, tt) - exp(-lambda * tau) * A(lambda + censor, tt)) / e
  truth <- cbind(ifelse(w == 1, risk_res, 0), ifelse(w == 1, rmst_res, 0))
  expect_equal(z$estimate, c(risk, rmst), tolerance = 1e-11)
  expect_equal(unname(z$EIF), truth, tolerance = 2e-9)
  expect_equal(z$covariance, cov(truth) / n, tolerance = 2e-9)
  expect_equal(z$se, sqrt(diag(cov(truth) / n)), tolerance = 2e-9)
  expect_lt(z$interpolation_error, 1e-8)
})

test_that("one-arm probability-one reduction reproduces the old native kernel", {
  rows <- cbind(0:5, 0, 0, 0, c(.1, .3, .5, .7, .9, 1), c(0, 0, 0, 0, 0, -1))
  args <- list(rates = rep(.8, 32), censor = .2, next_state = matrix(1L, 1, 1),
      subject_group = rep(0L, 6), records = rows, bins = 32L, horizon = 1, maxit = 40L, tolerance = 1e-8)
  old <- do.call(coherent4cast:::rc_cpp_continuous_target, args)
  new <- do.call(coherent4cast:::rc_cpp_causal_target, c(args, list(propensity = 1)))
  for (field in c("estimate", "initial", "EIF", "covariance", "rates", "score_trace", "values"))
    expect_equal(new[[field]], old[[field]], tolerance = 1e-13)
})

test_that("causal public API standardizes both arms and retains joint covariance", {
  d <- causal_demo_data(); set.seed(778); before <- .Random.seed
  fit <- rc_causal(~x, d, "treatment", 1, "censored", "year", propensity = .5,
       spec = list(complete_followup = TRUE, grid = 32L))
  expect_identical(before, .Random.seed)
  expect_s3_class(fit, "coherent4cast_causal")
  expect_true(fit$convergence$score_solved)
  expect_true(fit$convergence$valid)
  expect_equal(fit$n, nrow(d))
  for (arm in c("0", "1")) {
    a <- fit$arms[[arm]]
    direct <- rowMeans(a$values[, 1, 1, fit$group, drop = FALSE], dims = 3)
    expect_equal(unname(a$estimate), as.vector(direct), tolerance = 1e-12)
    expect_lt(max(abs(a$EIF_mean)), 1e-7 / sqrt(nrow(d)))
    expect_true(all(a$rates >= 0))
    expect_true(all(a$values[1, , , ] >= -1e-10 & a$values[1, , , ] <= 1 + 1e-10))
    rate <- (.55 + .25 * c(0, 1)) * exp(-.4 * as.numeric(arm))
    oracle <- c(mean(-expm1(-rate)), mean(-expm1(-rate) / rate))
    expect_equal(unname(a$estimate), oracle, tolerance = .015)
  }
  expect_equal(fit$joint_covariance, cov(fit$joint_EIF) / fit$n, tolerance = 1e-14)
  expect_equal(vcov(fit), cov(fit$arms[["1"]]$EIF - fit$arms[["0"]]$EIF) / fit$n, tolerance = 1e-14)
  effect <- predict(fit)
  expect_equal(effect$estimate, unname(coef(fit)))
  expect_true(all(effect$lower < effect$estimate & effect$upper > effect$estimate))
  ratio <- rc_causal_contrast(fit, "ratio")
  expect_equal(ratio$estimate, effect$arm1 / effect$arm0, tolerance = 1e-14)
  expect_true(all(ratio$lower > 0))
  log_ratio <- rc_causal_contrast(fit, "log_ratio")
  expect_equal(log_ratio$estimate, log(ratio$estimate))
  sim <- rc_causal_contrast(fit, interval = "multiplier", draws = 300, seed = 92)
  expect_identical(before, .Random.seed)
  expect_true(all(is.finite(sim$lower)))
  expect_true(all(is.na(predict(fit, interval = "none")$lower)))
  expect_equal(tidy(fit)$estimate, effect$estimate)
  expect_equal(glance(fit)$n, nrow(d))
  expect_equal(summary(fit)$estimates$estimate, effect$estimate)
  expect_output(print(fit), "Baseline-treatment")
  expect_error(logLik(fit), class = "coherent4cast_error")
  for (f in seq_along(fit$training_subjects))
    expect_false(any(fit$folds[fit$training_subjects[[f]]] == f))
  audit <- rc_audit_report(fit, include_session = FALSE)
  expect_match(paste(audit$model$assumptions, collapse = " "), "exchangeability")
  expect_identical(audit$fit_context$data_fingerprint$algorithm, "SHA-256")
  expect_match(audit$model$call, "audit_data_argument_redacted")
})

test_that("multitype counts and absorption remain coherent in causal fits", {
  d <- rc_simulate(180, horizon = 1, rates = c(.9, .6, .4), censor_rate = .15,
      dispersion = 0, seed = 12)
  d$treatment <- as.integer(d$id %% 2)
  d$censored <- as.integer(d$event_type == 0 & d$terminal == 0 & d$stop < 1)
  fit <- rc_causal(~1, d, "treatment", 1, "censored", "years",
      propensity = .5, spec = list(grid = 32L, working = "pooled"))
  expect_equal(length(fit$estimate), 4L)
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(fit$arms[[1]]$estimate >= 0))
  expect_true(all(fit$arms[[2]]$estimate >= 0))
  for (arm in fit$arms) {
    expect_true(arm$estimate["risk"] <= 1)
    expect_true(arm$estimate["rmst"] <= 1)
    absorbing <- dim(arm$values)[2]
    expect_true(all(abs(arm$values[1:2, absorbing, , ]) < 1e-12))
    expect_true(all(abs(arm$values[3, absorbing, , ] - 1) < 1e-12))
    expect_true(all(abs(arm$values[4, absorbing, , ]) < 1e-12))
  }
})

test_that("overlap, baseline, observation and unsupported nuisance inputs fail explicitly", {
  d <- causal_demo_data(80)
  fit <- function(z = d, ...) rc_causal(~x, z, "treatment", 1, "censored", "year", ...)
  expect_error(fit(propensity = .001), "overlap range")
  z <- d; z$treatment <- 0; expect_error(fit(z), "Both baseline")
  z <- d; z$treatment <- z$x; expect_error(fit(z), "training fold/stratum")
  z <- d; z$start <- .01; expect_error(fit(z), "contiguous")
  z <- d; z$stop[1] <- .05; z$terminal[1] <- 0; expect_error(fit(z), "explicit censor")
  expect_error(fit(spec = list(unknown = 1)), "Unknown")
  expect_error(fit(propensity = rep(.5, 3)), "fold-by-stratum")
  expect_error(fit(spec = list(initial_rates = list(1, 2))), "named lists")
  z <- d; z$censored[1] <- 1; z$terminal[1] <- 1; expect_error(fit(z), "without an event")
  z <- d; z$x[1] <- NA; expect_error(fit(z), "observed and fixed")
  z <- d; z$start[1] <- 0; z$stop[1] <- .1; z$terminal[1] <- 0
  extra <- z[1, ]; extra$start <- .1; extra$stop <- 1; extra$treatment <- 1 - extra$treatment
  z <- rbind(z, extra); z <- z[order(z$id, z$start), ]
  expect_error(fit(z), "remain fixed")
  expect_error(rc_causal_contrast(list()), class = "coherent4cast_error")
})

test_that("zero-event arms retain results and withhold regular intervals", {
  d <- causal_demo_data(80); d$stop <- 1; d$terminal <- 0
  expect_warning(fit <- rc_causal(~x, d, "treatment", 1, "censored", "year",
      spec = list(complete_followup = TRUE, maxit = 0)), "withheld")
  ci <- predict(fit)
  expect_true(all(is.na(ci$lower) & is.na(ci$upper)))
  expect_true(all(ci$inference_status == "withheld_nonregular"))
  expect_false(fit$convergence$valid)
})
