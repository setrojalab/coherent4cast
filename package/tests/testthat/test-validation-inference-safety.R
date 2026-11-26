test_that("validation rejects ignored landmarks and evaluates elapsed follow-up", {
  expect_error(validate(c(.8,.2),c(6,10),c(1,0),horizon=2,landmark=5),
    "Unsupported validation argument",class="coherentSurv_error")
  expect_error(validate(c(.8,.2),c(1,5),c(1,0),horizon=2,threshhold=.5),
    "Unsupported validation argument",class="coherentSurv_error")
  valid<-validate(c(.8,.2),time=c(6,10)-5,status=c(1,0),horizon=2)
  expect_equal(valid$metrics$estimate[valid$metrics$metric=="Brier"],.04)
  model<-structure(list(),class="coherentSurv")
  expect_error(validate(model,c(6,10),c(1,0),horizon=2,predicted=c(.8,.2),landmark=5),
    "Unsupported validation argument",class="coherentSurv_error")
})

test_that("separated calibration slopes are withheld in either direction", {
  outcome<-c(0,0,0,1,1,1)
  for(p in list(c(.1,.2,.3,.7,.8,.9),c(.9,.8,.7,.3,.2,.1),
                c(.1,.2,.5,.5,.8,.9))) {
    result<-validate(p,ifelse(outcome==1,1,3),outcome,horizon=2,
      censor_survival=list(case=rep(.8,6),horizon=rep(.7,6)))
    slope<-result$metrics[result$metrics$metric=="calibration_slope",]
    expect_true(is.na(slope$estimate));expect_false(slope$defined)
    expect_identical(result$calibration_status,"separated")
    expect_true(result$metrics$defined[result$metrics$metric=="CITL"])
    expect_true(result$metrics$defined[result$metrics$metric=="Brier"])
  }
  constant<-validate(rep(.4,6),ifelse(outcome==1,1,3),outcome,horizon=2)
  expect_identical(constant$calibration_status,"constant_predictions")
  expect_false(constant$metrics$defined[constant$metrics$metric=="calibration_slope"])
})

test_that("overlapping calibration outcomes retain the ordinary logistic estimate", {
  p<-c(.1,.2,.3,.4,.6,.7,.8,.9);y<-c(0,1,0,1,0,1,0,1)
  result<-validate(p,ifelse(y==1,1,3),y,horizon=2)
  oracle<-stats::glm(y~stats::qlogis(p),family=stats::quasibinomial())
  expect_equal(result$metrics$estimate[result$metrics$metric=="calibration_slope"],
    unname(stats::coef(oracle)[2]))
  expect_identical(result$calibration_status,"estimated")
})

test_that("ordinary summary and tidy retain nonregular inference status", {
  fit<-structure(list(coefficients=c(beta=1),covariance=matrix(.04),
    convergence=list(valid=TRUE,code="converged")),class="coherentSurv")
  regular<-tidy(fit)
  expect_equal(regular$conf.low,1-stats::qnorm(.975)*.2)
  expect_identical(regular$inference_status,"regular_asymptotic")
  for(code in c("parameter_boundary","nonregular_fit","quadrature_unresolved")) {
    fit$convergence<-list(valid=FALSE,code=code)
    result<-summary(fit);table<-tidy(fit)
    expect_equal(table$estimate,1);expect_equal(table$std.error,.2)
    expect_true(is.na(table$conf.low));expect_true(is.na(table$conf.high))
    expect_identical(table$inference_status,"unavailable_nonregular_fit")
    expect_identical(result$inference_status,"unavailable_nonregular_fit")
    expect_match(result$interval_scope,"withheld")
  }
})
