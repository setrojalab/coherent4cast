test_that("automatic fibres cannot manufacture exposure before delayed entry", {
  d<-data.frame(id=1,start=5,stop=6,event_type="1",terminal=0,x=0)
  expect_error(cs_coarsen_fibre(d,time_unit="years",event_types="1"),
    "delayed entry",class="coherentSurv_error")
  expect_error(cs_fit(~x,d,event_types="1",time_unit="years",observation="coarsened",frailty_diagnostics=FALSE),
    "delayed entry",class="coherentSurv_error")
})

test_that("exact forecasts reject an unobserved pre-entry history", {
  d<-cs_simulate(n=50,seed=711,time_unit="years")
  fit<-suppressWarnings(cs_fit(~x,d,event_types=c("1","2"),time_unit="years",frailty_diagnostics=FALSE,
    control=cs_control(phi_fixed=.4,history_cap=1L)))
  new<-data.frame(id=100,start=.2,stop=1,event_type="0",terminal=0,x=0)
  expect_error(predict(fit,newdata=new,landmark=.3,horizon=.2,type="risk"),
    "begin at time zero",class="coherentSurv_error")
})

test_that("conditional targeted prediction enforces its baseline experiment", {
  d<-read.csv(test_path("fixtures","continuous_target_subjects.csv"))
  fit<-suppressWarnings(cs_fit(~x,d,time_unit="year",targeting="continuous",
    targeting_spec=list(horizon=1,censor="censored",strata="x",folds="fold",grid=16,working="pooled")))
  subject<-d[1,,drop=FALSE];subject$start<-.2;subject$stop<-1
  subject$event_type<-"none";subject$terminal<-0;subject$censored<-0
  expect_error(predict(fit,newdata=subject,landmark=.3,horizon=.2,type="risk",interval="none"),
    "begin at time zero",class="coherentSurv_error")
  subject$start<-0
  expect_error(predict(fit,newdata=subject,landmark=.3,horizon=.2,type="unrecognized",interval="none"),
    "fitted target",class="coherentSurv_error")
  second<-subject;subject$stop<-.2;second$start<-.2;second$x<-1-subject$x
  expect_error(predict(fit,newdata=rbind(subject,second),landmark=.3,horizon=.2,type="risk",interval="none"),
    "baseline strata",class="coherentSurv_error")
  # A change after the landmark is unobserved at prediction and must not enter.
  second$start<-.5;subject$stop<-.5
  future<-rbind(subject,second);unchanged<-future;unchanged$x<-subject$x
  expect_equal(predict(fit,newdata=future,landmark=.3,horizon=.2,type="risk",interval="none")$estimate,
    predict(fit,newdata=unchanged,landmark=.3,horizon=.2,type="risk",interval="none")$estimate,tolerance=1e-13)
})

test_that("frailty-family comparisons retain identical subject case weights", {
  d<-cs_simulate(60,rates=c(.5,.2),beta=c(.1,.1),history=c(.2,.2),seed=45,time_unit="years")
  weights<-rep(c(1,3),30)
  fit<-suppressWarnings(cs_fit(~x,d,event_types="1",time_unit="years",weights=weights,frailty_diagnostics=FALSE,
    control=cs_control(history_cap=1L,quadrature=16L,quadrature_max=128L,quadrature_tolerance=1e-5,maxit=100L)))
  diagnostic<-cs_frailty_diagnose(fit,alternatives=TRUE)
  expect_match(diagnostic$comparison_scope,"descriptive weighted criterion")
  expect_match(diagnostic$null_variance_test$calibration,"fixed or exogenous")
  nd<-diagnostic$null_variance_test;dd<-fit$design;b<-nd$null_coefficients
  z<-cs_loglik(b,dd$unit,dd$exposure,dd$design,dd$event,phi_fixed=0)
  A1<-rowsum(dd$design*(dd$exposure*exp(drop(dd$design%*%b))),dd$unit,reorder=FALSE)
  analytic_B<-colSums(-weights*(z$counts-z$cumulative_hazard)*A1)
  expect_equal(nd$nuisance_sensitivity,unname(analytic_B),tolerance=1e-7)
  weighted_H<-apply(sweep(z$information_individual,3,weights,"*"),c(1,2),sum)
  manual_u<-weights*.5*((z$counts-z$cumulative_hazard)^2-z$counts)
  manual_S<-z$scores*weights
  manual_residual<-drop(manual_u+manual_S%*%solve(weighted_H,nd$nuisance_sensitivity))
  expect_equal(nd$effective_score_variance,sum(manual_residual^2),tolerance=1e-10)
  expect_equal(nd$z,sum(manual_u)/sqrt(sum(manual_residual^2)),tolerance=1e-10)
  expect_setequal(names(diagnostic$alternative_fits),c("lognormal","discrete"))
  for(alternative in diagnostic$alternative_fits) {
    expect_equal(alternative$case_weights,weights)
    dd<-alternative$design
    loglik<-cs_frailty_loglik(coef(alternative),dd$unit,dd$exposure,dd$design,dd$event,
      frailty=alternative$frailty,points=alternative$frailty_spec$points,quadrature=alternative$control$quadrature,
      information=FALSE)$subject_loglik
    expect_equal(alternative$loglik,sum(weights*loglik),tolerance=1e-10)
  }
})

test_that("a nonregular null nuisance fit cannot produce a boundary p value", {
  d<-cs_simulate(40,rates=c(.5,.2),beta=c(0,0),history=c(0,0),seed=92,time_unit="years")
  fit<-suppressWarnings(cs_fit(~1,d,event_types=c("1","unobserved"),time_unit="years",frailty_diagnostics=FALSE,
    control=cs_control(phi_fixed=.4,history_cap=1L,maxit=100L)))
  diagnostic<-cs_frailty_diagnose(fit)$null_variance_test
  expect_false(diagnostic$null_regular)
  expect_true(is.na(diagnostic$p.value))
  expect_match(diagnostic$calibration,"p value is withheld")
})
