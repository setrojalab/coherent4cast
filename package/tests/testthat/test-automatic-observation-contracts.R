test_that("automatic fibres enforce observed-event and observation-window contracts", {
  h <- data.frame(id = "A", event_time_lower = .2, event_time_upper = .6,
                  event_type = "1", terminal = 0, followup = 1, x = 0)
  fibre <- function(d = h, ...) rc_coarsen_fibre(d, event_types = "1", time_unit = "years", ...)
  expect_error(rc_coarsen_fibre(h), "time_unit", class = "coherent4cast_error")
  expect_error(fibre(h[FALSE, ]), "fibre_data", class = "coherent4cast_error")
  expect_error(fibre(h[, setdiff(names(h), "terminal")]), "fibre_columns", class = "coherent4cast_error")
  expect_error(fibre(rounding_unit = 0), "rounding_unit", class = "coherent4cast_error")
  expect_error(fibre(seed = -1), "fibre_control", class = "coherent4cast_error")
  expect_error(fibre(breaks = c(1, Inf)), "breaks", class = "coherent4cast_error")
  bad <- h; bad$id <- NA_character_
  expect_error(fibre(bad), "fibre_unit", class = "coherent4cast_error")
  bad <- h; bad$terminal <- 2
  expect_error(fibre(bad), "fibre_event", class = "coherent4cast_error")
  expect_error(fibre(nonrecurrent = "unknown"), "nonrecurrent", class = "coherent4cast_error")
  expect_error(fibre(covariates = "followup"), "fibre_covariates", class = "coherent4cast_error")
  bad <- h; bad$.rc_source_unit <- "reserved"
  expect_error(fibre(bad), "fibre_covariates", class = "coherent4cast_error")
  bad <- h; bad$event_time_lower <- "not a time"
  expect_error(fibre(bad), "fibre_time", class = "coherent4cast_error")
  bad <- h; bad$followup <- 0
  expect_error(fibre(bad), "fibre_followup", class = "coherent4cast_error")
  rounded <- h; rounded$event_time_lower <- rounded$event_time_upper
  window <- h; window$event_time_lower <- .5
  a <- fibre(rounded, rounding_unit = .1, nodes = 2)
  b <- fibre(window, nodes = 2)
  expect_equal(a$data$stop, b$data$stop, tolerance = 1e-14)
  expect_equal(a$coarsening$log_weights, b$coarsening$log_weights, tolerance = 1e-14)
  expect_error(rc_fibre_mcse(a, Inf), "fibre_mcse", class = "coherent4cast_error")
})

test_that("colon fibres preserve surveillance truncation and training-only imputation", {
  raw <- data.frame(id = c(1, 1, 2, 2), etype = c(1, 2, 1, 2),
                    time = c(10, 20, 15, 30), status = c(0, 1, 1, 0),
                    rx = c("Obs", "Obs", "Lev", "Lev"), sex = c(0, 0, 1, 1),
                    age = c(60, 60, 65, 65), node4 = c(0, 0, 1, 1),
                    nodes = c(NA_real_, NA_real_, 5, 5))
  z <- rc_colon_fibre(raw, nodes = 2, covariate = "nodes", node_median = 7)
  expect_identical(z$metadata$surveillance_truncated_ids, "1")
  expect_equal(z$metadata$source_deaths, 1)
  expect_equal(sum(z$data$terminal), 0)
  expect_equal(unique(z$data$nodes_per_5[z$data$.rc_source_unit == "1"]), 7/5)
  expect_equal(unique(z$data$nodes_missing[z$data$.rc_source_unit == "1"]), 1)
  estimated <- rc_colon_fibre(raw, nodes = 2, covariate = "nodes")
  expect_equal(estimated$metadata$node_median, 5)
  expect_error(rc_colon_fibre(raw, covariate = "nodes", node_median = -1), "node_median", class = "coherent4cast_error")
  bad <- raw; bad$etype[1] <- NA_integer_
  expect_error(rc_colon_fibre(bad), "adapter_id", class = "coherent4cast_error")
  expect_error(rc_colon_fibre(raw[-1, ]), "adapter_pair", class = "coherent4cast_error")
  bad <- raw; bad$time[1] <- -1
  expect_error(rc_colon_fibre(bad), "adapter_time", class = "coherent4cast_error")
  bad <- raw; bad$age[1] <- 61
  expect_error(rc_colon_fibre(bad), "adapter_covariates", class = "coherent4cast_error")
})

test_that("history-only prediction rejects invalid measurement provenance", {
  training <- rc_simulate(40, horizon = 2, rates = c(.5, .2),
                           beta = c(-.3, .2), history = c(.2, .2), seed = 918)
  fit <- suppressWarnings(coherent4cast(~x, training, time_unit = "years",
                                     control = coherent4cast_control(phi_fixed = .4),
                                     frailty_diagnostics = FALSE))
  h <- data.frame(id = "A", event_time_lower = .2, event_time_upper = .6,
                  event_type = "1", terminal = 0, followup = 1, x = 0)
  forecast <- function(d = h, ...) rc_predict_coarsened(fit, d, 1, .2, type = "all",
                   quadrature = 4, coarsening_control = list(nodes = 2), ...)
  expect_error(rc_predict_coarsened(NULL, h, 1, .2), "fit", class = "coherent4cast_error")
  expect_error(forecast(h[FALSE, ]), "coarsened_history", class = "coherent4cast_error")
  expect_error(rc_predict_coarsened(fit, h, 0, .2), "landmark", class = "coherent4cast_error")
  expect_error(rc_predict_coarsened(fit, h, 1, .2, coarsening_control = list(2)), "coarsening_control", class = "coherent4cast_error")
  expect_error(forecast(unknown_setting = TRUE), "unused_argument", class = "coherent4cast_error")
  expect_error(forecast(h[, setdiff(names(h), "terminal")]), "coarsened_history", class = "coherent4cast_error")
  expect_error(forecast(newdata_future = list(time = 1)), "covariate_history", class = "coherent4cast_error")
  expect_error(forecast(newdata_future = data.frame(id = "unknown", time = 1.1, x = 0)), "covariate_history", class = "coherent4cast_error")
  expect_error(forecast(covariate_history = data.frame(id = "A", time = -.1, x = 0)), "nonnegative", class = "coherent4cast_error")
  expect_error(forecast(covariate_history = data.frame(id = "A", time = 2, x = 0)), "actual measurement", class = "coherent4cast_error")
  exact <- data.frame(id = "A", start = 0, stop = 1, event_type = 0, terminal = 0, x = 0)
  a <- forecast(exact)
  b <- predict(fit, exact, landmark = 1, horizon = .2, quadrature = 4, type = "all")
  expect_equal(a$estimate, b$estimate[match(a$target, b$target)], tolerance = 1e-12)
  exact$stop <- .9
  expect_error(forecast(exact), "path_history", class = "coherent4cast_error")
  exact$stop <- 1.1
  expect_error(forecast(exact), "path_future", class = "coherent4cast_error")
  future <- data.frame(id = "A", time = 1.1, x = 1)
  evolving <- forecast(covariate_evolution = "locf", newdata_future = future)
  expect_true(any(grepl("Future covariate rule: locf", attr(evolving, "assumptions"), fixed = TRUE)))
  expect_false(any(grepl("fixed at their landmark", attr(evolving, "assumptions"), fixed = TRUE)))
})
