test_that("native entry points are registered and dynamic lookup is disabled", {
  dll <- getLoadedDLLs()[["coherent4cast"]]
  expect_false(dll[["dynamicLookup"]])
  calls <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expected <- c(`_coherent4cast_cpp_rc_loglik` = 7L,
                `_coherent4cast_cpp_rc_coarsen` = 5L,
                `_coherent4cast_cpp_rc_propagate` = 6L,
                `_coherent4cast_rc_cpp_target` = 7L,
                `_coherent4cast_rc_cpp_bounds` = 8L,
                `_coherent4cast_rc_cpp_conformal` = 5L,
                `_coherent4cast_rc_cpp_pivot` = 3L,
                `_coherent4cast_rc_cpp_history` = 7L,
                `_coherent4cast_rc_cpp_cf` = 5L,
                `_coherent4cast_rc_cpp_continuous_target` = 9L,
                `_coherent4cast_rc_cpp_causal_target` = 10L,
                `_coherent4cast_cpp_rc_coarsen_fibre` = 12L,
                `_coherent4cast_cpp_rc_mixture_lik` = 12L,
                `_coherent4cast_cpp_rc_propagate_sparse` = 9L,
                `_coherent4cast_cpp_rc_bounds_sparse` = 10L)
  expect_equal(sort(names(calls)), sort(names(expected)))
  actual <- vapply(calls, function(x) x$numParameters, integer(1))
  expect_equal(actual[names(expected)], expected)
})

