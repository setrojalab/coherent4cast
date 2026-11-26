test_that("native entry points are registered and dynamic lookup is disabled", {
  dll <- getLoadedDLLs()[["coherentSurv"]]
  expect_false(dll[["dynamicLookup"]])
  calls <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expected <- c(`_coherentSurv_cpp_cs_loglik` = 7L,
                `_coherentSurv_cpp_cs_coarsen` = 5L,
                `_coherentSurv_cpp_cs_propagate` = 6L,
                `_coherentSurv_cs_cpp_target` = 7L,
                `_coherentSurv_cs_cpp_bounds` = 8L,
                `_coherentSurv_cs_cpp_conformal` = 5L,
                `_coherentSurv_cs_cpp_pivot` = 3L,
                `_coherentSurv_cs_cpp_history` = 7L,
                `_coherentSurv_cs_cpp_cf` = 5L,
                `_coherentSurv_cs_cpp_continuous_target` = 9L,
                `_coherentSurv_cs_cpp_causal_target` = 10L,
                `_coherentSurv_cpp_cs_coarsen_fibre` = 12L,
                `_coherentSurv_cpp_cs_mixture_lik` = 12L,
                `_coherentSurv_cpp_cs_propagate_sparse` = 9L,
                `_coherentSurv_cpp_cs_bounds_sparse` = 10L)
  expect_equal(sort(names(calls)), sort(names(expected)))
  actual <- vapply(calls, function(x) x$numParameters, integer(1))
  expect_equal(actual[names(expected)], expected)
})

