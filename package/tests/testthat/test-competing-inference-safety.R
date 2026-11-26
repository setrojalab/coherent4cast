competing_safety_data <- function() data.frame(id=1:12,
  time=c(.2,.3,.4,.5,.6,.8,1,1.2,1.3,1.5,1.7,2),
  status=c(1,2,0,1,0,2,1,0,2,0,1,0))

test_that("regular competing summaries preserve coefficient columns and limits", {
  fit<-cs_competing(~1,competing_safety_data(),time_unit="years",frailty_diagnostics=FALSE)
  table<-tidy(fit);se<-sqrt(diag(vcov(fit)))
  expect_true(fit$convergence$valid)
  expect_identical(names(table)[1:5],c("term","estimate","std.error","conf.low","conf.high"))
  expect_equal(table$estimate,unname(coef(fit)))
  expect_equal(unname(table$std.error),unname(se))
  expect_equal(unname(table$conf.low),unname(coef(fit)-stats::qnorm(.975)*se))
  expect_equal(unname(table$conf.high),unname(coef(fit)+stats::qnorm(.975)*se))
  expect_true(all(table$inference_status=="regular_asymptotic"))
  expect_identical(summary(fit)$inference_status,"regular_asymptotic")
  expect_equal(summary(fit)$coefficients,table)
})

test_that("each competing component can withhold joint Wald limits", {
  fit<-cs_competing(~1,competing_safety_data(),time_unit="years",frailty_diagnostics=FALSE)
  regular<-tidy(fit)
  for(k in seq_along(fit$fits)) {
    failed<-fit
    # A retained wrapper flag must not override a failed component diagnostic.
    failed$fits[[k]]$convergence<-list(valid=FALSE,code="parameter_boundary")
    expect_true(failed$convergence$valid)
    table<-tidy(failed);result<-summary(failed)
    expect_equal(table[,1:3],regular[,1:3])
    expect_true(all(is.na(table$conf.low)))
    expect_true(all(is.na(table$conf.high)))
    expect_true(all(table$inference_status=="unavailable_nonregular_fit"))
    expect_identical(result$inference_status,"unavailable_nonregular_fit")
    expect_equal(result$coefficients,table)
    expect_match(result$interval_scope,"withheld")
  }
  failed<-fit;failed$convergence$valid<-FALSE
  expect_true(all(is.na(tidy(failed)$conf.low)))
  expect_identical(summary(failed)$inference_status,"unavailable_nonregular_fit")
})

test_that("an actually constrained competing fit withholds regular limits", {
  fit<-suppressWarnings(cs_competing(~1,competing_safety_data(),time_unit="years",
    control=cs_control(coefficient_bounds=c(-.1,.1)),frailty_diagnostics=FALSE))
  expect_false(fit$convergence$valid)
  expect_true(any(vapply(fit$fits,function(f)isTRUE(f$convergence$at_boundary),logical(1))))
  table<-tidy(fit)
  expect_true(all(is.finite(table$estimate)))
  expect_true(all(is.na(table$conf.low)))
  expect_true(all(is.na(table$conf.high)))
  expect_identical(summary(fit)$inference_status,"unavailable_nonregular_fit")
})
