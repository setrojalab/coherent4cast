test_that("lognormal refinement is audited at the returned optimized parameter", {
  d<-cs_simulate(60,seed=835)
  f<-suppressWarnings(cs_fit(~1,d,time_unit="years",frailty="lognormal",frailty_diagnostics=FALSE,
    control=cs_control(history_cap=0,phi_fixed=.5,quadrature=4,quadrature_max=128)))
  audit<-tail(f$audit$quadrature_refinement,1);expect_gt(nrow(audit),0)
  dd<-f$design
  run<-function(q)cs_frailty_loglik(coef(f),dd$unit,dd$exposure,dd$design,dd$event,
    "lognormal",phi_fixed=.5,quadrature=q,information=FALSE)
  low<-run(audit$nodes_low);high<-run(audit$nodes_high)
  expect_lt(abs(audit$max_subject_loglik_change-max(abs(high$subject_loglik-low$subject_loglik))),1e-12)
  expect_lt(abs(audit$max_mean_score_change-max(abs(high$gradient-low$gradient))),1e-12)
  expect_equal(f$control$quadrature,audit$nodes_high)
  # Starting at the maximum still compares with a coarser rule at the final fit.
  g<-suppressWarnings(cs_fit(~1,d,time_unit="years",frailty="lognormal",frailty_diagnostics=FALSE,
    control=cs_control(history_cap=0,phi_fixed=.5,quadrature=32,quadrature_max=32)))
  expect_equal(g$audit$quadrature_refinement$nodes_low,16L)
  expect_equal(g$audit$quadrature_refinement$nodes_high,32L)
  expect_identical(g$audit$quadrature_stable,
    g$audit$quadrature_refinement$max_subject_loglik_change<=1e-7 && g$audit$quadrature_refinement$max_mean_score_change<=1e-7)
})
