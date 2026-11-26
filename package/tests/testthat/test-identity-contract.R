test_that("the public namespace and native registration use the current identity", {
  exported <- getNamespaceExports("coherentSurv")
  required <- c("cs_fit", "cs_control", "cs_as_trajectory", "cs_build_history",
                "cs_surv", "cs_competing", "cs_causal", "cs_causal_contrast",
                "cs_bounds", "cs_target", "cs_audit_report")
  expect_true(all(required %in% exported))
  former_prefix <- paste0("r", "c", "_")
  expect_false(any(startsWith(exported, former_prefix)))
  registered <- names(getDLLRegisteredRoutines("coherentSurv")[[".Call"]])
  expect_gt(length(registered), 0L)
  expect_true(all(startsWith(registered, "_coherentSurv_")))
  expect_false(any(grepl(former_prefix, registered, fixed = TRUE)))
  for (generic in c("predict", "summary", "coef", "vcov", "logLik")) {
    expect_true(is.function(getS3method(generic, "coherentSurv", optional = TRUE)))
  }
})
