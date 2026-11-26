test_that("environment diagnostics preserve analytic references and unknown provenance", {
  report <- cs_validate_env()
  expect_s3_class(report, "coherentSurv_environment")
  expect_identical(report$report_type, "coherentSurv_environment_diagnostic")
  expect_identical(report$operational_status, "passed")
  expect_true(all(report$checks$status %in% c("passed", "failed", "unavailable")))
  expect_equal(sum(report$checks$category == "operational"), 4)
  numerical <- report$checks[report$checks$category == "operational", ]
  expect_true(all(numerical$maximum_absolute_error <= 1e-12))
  expect_true(all(is.na(report$checks$elapsed_seconds)))
  expect_identical(report$installation$build_compiler$status, "unavailable")
  expect_null(report$installation$build_compiler$compiler)
  expect_null(report$installation$build_compiler$flags)
  expect_identical(report$status, "unavailable")
  expect_equal(report$installation$package$version, as.character(getNamespaceVersion("coherentSurv")))
  expect_identical(report$installation$native$dynamic_lookup, FALSE)
  expect_match(report$installation$native$file_sha256, "^[0-9a-f]{64}$")
  expect_equal(report$fixtures$dense_absorption$observed[1:2],
               c(probability1 = exp(-.6), probability2 = -expm1(-.6)), tolerance = 1e-12)
  expect_equal(unname(report$fixtures$louis_mixture$observed["information"]), 4.96, tolerance = 1e-12)
  expect_match(paste(report$interpretation, collapse = " "), "unsigned and mutable")
})

test_that("environment fixtures do not use RNG and optional timings remain descriptive", {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (existed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(47); before <- .Random.seed
  one <- cs_validate_env(timings = TRUE)
  expect_identical(.Random.seed, before)
  two <- cs_validate_env()
  for (name in names(one$fixtures)) {
    expect_equal(one$fixtures[[name]]$observed, two$fixtures[[name]]$observed, tolerance = 0)
    expect_gt(one$fixtures[[name]]$result_object_bytes, 0)
    expect_null(two$fixtures[[name]]$result_object_bytes)
  }
  observed <- one$checks$elapsed_seconds[one$checks$category == "operational"]
  expect_true(all(is.finite(observed) & observed >= 0))
  expect_match(paste(one$interpretation, collapse = " "), "not peak resident memory")
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(cs_validate_env())
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("numerical errors warnings malformed results and disagreement cannot pass", {
  run <- get(".cs_env_fixture", envir = environment(cs_validate_env))
  failure <- run("failure", function() stop("controlled fixture failure"), c(value = 1))
  expect_identical(failure$check$status, "failed")
  expect_match(failure$check$detail, "controlled fixture failure")
  warning <- run("warning", function() { warning("controlled warning"); c(value = 1) }, c(value = 1))
  expect_identical(warning$check$status, "failed")
  expect_match(warning$check$detail, "controlled warning")
  for (value in list(c(value = NaN), c(value = 2), c(other = 1), c(value = 1, extra = 2))) {
    result <- run("invalid", function() value, c(value = 1))
    expect_identical(result$check$status, "failed")
  }
  good <- run("good", function() c(value = 1), c(value = 1))
  expect_identical(good$check$status, "passed")
})

test_that("an installation failure takes precedence over unavailable provenance", {
  original <- get(".cs_env_metadata", envir = environment(cs_validate_env))
  isolated <- new.env(parent = environment(cs_validate_env))
  isolated$.cs_env_metadata <- function() {
    x <- original(); x$native$registered_calls <- character(); x
  }
  diagnose <- cs_validate_env; environment(diagnose) <- isolated
  report <- diagnose()
  expect_identical(report$status, "failed")
  expect_identical(report$operational_status, "passed")
  expect_identical(report$checks$status[report$checks$check == "registered_native_entry_points"], "failed")
})

test_that("environment exports are portable and invalid options fail explicitly", {
  json_file <- tempfile(fileext = ".json"); text_file <- tempfile(fileext = ".txt")
  rds_file <- tempfile(fileext = ".rds")
  on.exit(unlink(c(json_file, text_file, rds_file)), add = TRUE)
  text <- cs_validate_env(format = "json", file = json_file)
  parsed <- jsonlite::fromJSON(text)
  expect_identical(parsed$operational_status, "passed")
  expect_identical(parsed$installation$build_compiler$status, "unavailable")
  expect_null(parsed$installation$build_compiler$compiler)
  expect_true(all(is.na(parsed$checks$elapsed_seconds)))
  expect_equal(jsonlite::fromJSON(json_file), parsed)
  lines <- cs_validate_env(format = "text", file = text_file)
  expect_identical(readLines(text_file, encoding = "UTF-8"), lines)
  expect_match(paste(lines, collapse = " "), "not IQ/OQ/PQ certification", fixed = TRUE)
  report <- cs_validate_env(file = rds_file)
  expect_identical(readRDS(rds_file), report)
  printer <- get("print.coherentSurv_environment", envir = environment(cs_validate_env))
  expect_output(expect_identical(printer(report), report), "environment diagnostic")
  for (bad in list(NA, 1, c(TRUE, FALSE), NULL))
    expect_error(cs_validate_env(timings = bad), "timings must", class = "coherentSurv_error")
  for (bad in list(NA_character_, "", c("a", "b"), 1))
    expect_error(cs_validate_env(file = bad), "output path", class = "coherentSurv_error")
  expect_error(cs_validate_env(format = "xml"), "arg")
})
