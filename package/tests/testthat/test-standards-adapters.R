standards_adam_fixture <- function() data.frame(
  USUBJID = rep(c("001", "002", "003"), 2),
  PARAMCD = rep(c("OS", "OTHER"), each = 3), AVAL = c(5, 10, 7, 2, 3, 4),
  CNSR = c(0, 1, 2, 1, 0, 0), AVALU = "DAYS", stringsAsFactors = FALSE)

standards_adam <- function(d = standards_adam_fixture(), ...) cs_from_adam(d,
  paramcd = "OS", event_codes = 0, censor_codes = c(1, 2), time_unit = "days",
  duration_definition = "Synthetic elapsed days from process origin",
  ascertainment = "Complete synthetic single-endpoint follow-up",
  provenance = "Synthetic test fixture", ...)

standards_omop_fixture <- function() {
  origin <- rep(as.Date("2020-01-01"), 3)
  list(person = data.frame(person_id = c("001", "002", "003")),
    observation_period = data.frame(person_id = c("001", "002", "003"),
      observation_period_start_date = origin, observation_period_end_date = origin + c(12, 15, 10)),
    events = data.frame(person_id = c("002", "001", "001"),
      condition_concept_id = c("EXAMPLE_A", "EXAMPLE_B", "EXAMPLE_A"),
      condition_start_date = origin + c(4, 8, 3)),
    concept_map = c(EXAMPLE_A = "A", EXAMPLE_B = "B"), origin = origin,
    time_unit = "days", timing = "exact_day_boundaries",
    ascertainment = "Complete synthetic recurrent and death ascertainment from origin; one row per episode",
    provenance = "Synthetic reviewed test map, not clinical concept IDs",
    death = data.frame(person_id = "001", death_date = origin[1] + 10))
}

test_that("ADaM parameter selection preserves durations, IDs and nonbinary censor codes", {
  d <- standards_adam_fixture(); before <- serialize(d, NULL)
  z <- standards_adam(d)
  expect_s3_class(z, "cs_trajectory")
  expect_identical(z$id, c("001", "002", "003"))
  expect_equal(z$stop, c(5, 10, 7))
  expect_identical(z$terminal, c(1L, 0L, 0L))
  expect_identical(z$source_censor_code, c("0", "1", "2"))
  expect_identical(z$source_row, 1:3)
  expect_identical(attr(z, "cs_mapping")$paramcd, "OS")
  expect_identical(attr(cs_as_trajectory(z, event_types = character(), time_unit = "days"), "cs_mapping"), attr(z, "cs_mapping"))
  expect_identical(serialize(d, NULL), before)
  expect_true(attr(z, "observation_process_details")$continuous_observation_from_zero)
})

test_that("ADaM joins selected ADSL covariates by exact subject key", {
  a <- data.frame(USUBJID = c("003", "001", "002", "004"), age = c(30, 10, NA, 40))
  before <- serialize(a, NULL)
  z <- standards_adam(adsl = a, covariates = "age")
  expect_equal(z$age, c(10, NA, 30))
  expect_identical(serialize(a, NULL), before)
  a$USUBJID[4] <- "001"
  expect_error(standards_adam(adsl = a), "mapping_subject")
  expect_error(standards_adam(adsl = a[1, , drop = FALSE]), "mapping_subject")
  expect_error(standards_adam(covariates = "age"), "mapping_covariates")
  expect_error(standards_adam(adsl = data.frame(USUBJID = c("001", "002", "003"), start = 1), covariates = "start"), "mapping_covariates")
})

test_that("ADaM refuses ambiguous parameter, duration and censor mappings", {
  d <- standards_adam_fixture()
  expect_error(standards_adam(rbind(d, d[1, ])), "mapping_subject")
  bad <- d; bad$CNSR[3] <- 9
  expect_error(standards_adam(bad), "mapping_censor")
  args <- list(adtte = d, paramcd = "OS", event_codes = 0, censor_codes = c(0, 2),
    time_unit = "days", duration_definition = "elapsed days", ascertainment = "synthetic", provenance = "synthetic")
  expect_error(do.call(cs_from_adam, args), "mapping_censor")
  args$censor_codes <- c(1, 2); args$paramcd <- "ABSENT"
  expect_error(do.call(cs_from_adam, args), "mapping_parameter")
  for (value in c(0, -1, Inf, NA_real_)) {
    bad <- d; bad$AVAL[1] <- value
    expect_error(standards_adam(bad), "mapping_duration")
  }
  bad <- d; bad$AVALU[1] <- "years"
  expect_error(standards_adam(bad), "mapping_duration")
  bad <- d; bad$USUBJID[1] <- NA_character_
  expect_error(standards_adam(bad), "mapping_key")
})

test_that("OMOP complete histories preserve identifiers, event order and source inputs", {
  a <- standards_omop_fixture(); before <- serialize(a, NULL)
  z <- do.call(cs_from_omop, a)
  expect_s3_class(z, "cs_trajectory")
  expect_identical(unique(z$id), c("001", "002", "003"))
  expect_equal(z$start, c(0, 3, 8, 0, 4, 0))
  expect_equal(z$stop, c(3, 8, 10, 4, 15, 10))
  expect_identical(z$event_type, c("A", "B", "0", "A", "0", "0"))
  expect_identical(z$terminal, c(0L, 0L, 1L, 0L, 0L, 0L))
  expect_identical(z$source_row, c(3L, 2L, NA_integer_, 1L, NA_integer_, NA_integer_))
  expect_equal(sum(z$stop - z$start), 35)
  expect_true(attr(z, "observation_process_details")$continuous_observation_from_zero)
  expect_identical(attr(z, "cs_mapping")$concept_map, a$concept_map)
  expect_identical(attr(cs_as_trajectory(z, event_types = c("A", "B"), time_unit = "days"), "cs_mapping"), attr(z, "cs_mapping"))
  expect_identical(serialize(a, NULL), before)
  a$time_unit <- "years"
  years <- do.call(cs_from_omop, a)
  expect_equal(years$stop * 365.25, z$stop, tolerance = 1e-13)
})

test_that("OMOP rejects unknown concepts, people and unreviewed duplicate episodes", {
  a <- standards_omop_fixture(); a$events$condition_concept_id[1] <- "UNKNOWN"
  expect_error(do.call(cs_from_omop, a), "mapping_concepts")
  a <- standards_omop_fixture(); a$events$person_id[1] <- "999"
  expect_error(do.call(cs_from_omop, a), "mapping_subject")
  a <- standards_omop_fixture(); a$events <- rbind(a$events, a$events[1, ])
  expect_error(do.call(cs_from_omop, a), "mapping_duplicate")
  a <- standards_omop_fixture(); names(a$concept_map) <- c("EXAMPLE_A", "EXAMPLE_A")
  expect_error(do.call(cs_from_omop, a), "mapping_concepts")
  a <- standards_omop_fixture(); a$concept_map[1] <- "0"
  expect_error(do.call(cs_from_omop, a), "mapping_concepts")
})

test_that("OMOP does not invent timing for ties or violate terminal absorption", {
  a <- standards_omop_fixture(); a$events$condition_start_date[2] <- a$events$condition_start_date[3]
  expect_error(do.call(cs_from_omop, a), "mapping_tie")
  for (day in c(10, 11)) {
    a <- standards_omop_fixture(); a$events$condition_start_date[2] <- a$origin[1] + day
    expect_error(do.call(cs_from_omop, a), "mapping_absorption")
  }
  a <- standards_omop_fixture(); a$events$condition_start_date[1] <- a$origin[2] + 15
  expect_error(do.call(cs_from_omop, a), "mapping_absorption")
  a <- standards_omop_fixture(); a$timing <- "rounded_dates"
  expect_error(do.call(cs_from_omop, a), "mapping_timing")
})

test_that("OMOP requires complete periods from the declared process origin", {
  for (day in c(-1, 0, 16)) {
    a <- standards_omop_fixture(); a$events$condition_start_date[1] <- a$origin[2] + day
    expect_error(do.call(cs_from_omop, a), "mapping_containment")
  }
  a <- standards_omop_fixture(); a$observation_period <- rbind(a$observation_period, a$observation_period[1, ])
  expect_error(do.call(cs_from_omop, a), "mapping_period")
  a <- standards_omop_fixture(); a$observation_period <- a$observation_period[-1, ]
  expect_error(do.call(cs_from_omop, a), "mapping_period")
  a <- standards_omop_fixture(); a$origin[1] <- a$origin[1] - 1
  expect_error(do.call(cs_from_omop, a), "mapping_origin")
  a <- standards_omop_fixture(); a$origin <- a$origin[1]
  expect_error(do.call(cs_from_omop, a), "mapping_origin")
  a <- standards_omop_fixture(); a$observation_period$observation_period_end_date[1] <- a$origin[1]
  expect_error(do.call(cs_from_omop, a), "mapping_origin")
  for (day in c(0, 13)) {
    a <- standards_omop_fixture(); a$death$death_date <- a$origin[1] + day
    expect_error(do.call(cs_from_omop, a), "mapping_containment")
  }
  a <- standards_omop_fixture(); a$death <- rbind(a$death, a$death)
  expect_error(do.call(cs_from_omop, a), "mapping_subject")
})

test_that("OMOP handles declared zero-event extracts without inventing events", {
  a <- standards_omop_fixture(); a$events <- a$events[FALSE, ]; a$death <- NULL
  z <- do.call(cs_from_omop, a)
  expect_equal(nrow(z), 3)
  expect_equal(z$stop, c(12, 15, 10))
  expect_identical(z$event_type, rep("0", 3))
  expect_identical(z$terminal, rep(0L, 3))
})

test_that("standards mappings require declarations and reject unparsed dates", {
  a <- standards_omop_fixture(); a$ascertainment <- " "
  expect_error(do.call(cs_from_omop, a), "mapping_declaration")
  a <- standards_omop_fixture(); a$provenance <- NA_character_
  expect_error(do.call(cs_from_omop, a), "mapping_declaration")
  a <- standards_omop_fixture(); a$events$condition_start_date <- as.character(a$events$condition_start_date)
  expect_error(do.call(cs_from_omop, a), "mapping_date")
  a <- standards_omop_fixture(); a$person$person_id <- c(2^54, 2, 3)
  expect_error(do.call(cs_from_omop, a), "mapping_key")
})
