test_that("native entry points are registered and dynamic lookup is disabled", {
  dll <- getLoadedDLLs()[["recurcast"]]
  expect_false(dll[["dynamicLookup"]])
  calls <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expected <- c(`_recurcast_cpp_rc_loglik` = 7L,
                `_recurcast_cpp_rc_coarsen` = 5L,
                `_recurcast_cpp_rc_propagate` = 6L,
                `_recurcast_rc_cpp_target` = 7L,
                `_recurcast_rc_cpp_bounds` = 8L,
                `_recurcast_rc_cpp_conformal` = 5L,
                `_recurcast_rc_cpp_pivot` = 3L,
                `_recurcast_rc_cpp_history` = 7L,
                `_recurcast_rc_cpp_cf` = 5L,
                `_recurcast_rc_cpp_continuous_target` = 9L,
                `_recurcast_cpp_rc_coarsen_fibre` = 12L,
                `_recurcast_cpp_rc_mixture_lik` = 12L,
                `_recurcast_cpp_rc_propagate_sparse` = 9L,
                `_recurcast_cpp_rc_bounds_sparse` = 10L)
  expect_equal(sort(names(calls)), sort(names(expected)))
  actual <- vapply(calls, function(x) x$numParameters, integer(1))
  expect_equal(actual[names(expected)], expected)
})
