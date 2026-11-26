test_that("ordinary confint supports parm and level without bypassing regularity", {
  object<-structure(list(coefficients=c(a=1,b=2),covariance=diag(c(.04,.09)),
    convergence=list(valid=TRUE)),class="coherentSurv")
  result<-confint(object,parm=c("b","a"),level=.9)
  expect_identical(rownames(result),c("b","a"))
  expect_equal(as.numeric(result[,1]),c(2,1)-qnorm(.95)*c(.3,.2))
  expect_equal(as.numeric(result[,2]),c(2,1)+qnorm(.95)*c(.3,.2))
  expect_identical(colnames(result),c("5 %","95 %"))
  expect_equal(confint(object,parm=2),confint(object,parm="b"))
  expect_equal(confint(object,parm=-1),confint(object,parm="b"))
  expect_identical(dim(confint(object,parm=integer())),c(0L,2L))
  object$convergence$valid<-FALSE
  expect_true(all(is.na(confint(object))))
  expect_true(all(attr(confint(object),"inference_status")=="unavailable_nonregular_fit"))
  object$convergence$valid<-TRUE;object$covariance[1,1]<-0
  expect_true(all(is.na(confint(object,parm="a"))))
  expect_true(all(is.finite(confint(object,parm="b"))))
  for(parm in list("unknown",NA_character_,3,1.5,c(-1,2),TRUE))
    expect_error(confint(object,parm=parm),class="coherentSurv_error")
  expect_error(confint(object,level=1),class="coherentSurv_error")
  expect_error(confint(object,unknown=TRUE),class="coherentSurv_error")
})

test_that("survival and competing confint respect component inference guards", {
  d<-data.frame(time=c(.2,.3,.4,.5,.6,.8,1,1.2,1.3,1.5,1.7,2),
    status=c(1,2,0,1,0,2,1,0,2,0,1,0))
  survival<-cs_surv(~1,transform(d,status=as.integer(status>0)),time_unit="years",frailty_diagnostics=FALSE)
  expect_equal(as.numeric(confint(survival)[,1]),summary(survival)$coefficients$conf.low)
  competing<-cs_competing(~1,d,time_unit="years",frailty_diagnostics=FALSE)
  expect_equal(as.numeric(confint(competing)[,2]),tidy(competing)$conf.high)
  for(k in seq_along(competing$fits)) {
    invalid<-competing;invalid$fits[[k]]$convergence$valid<-FALSE
    expect_true(all(is.na(confint(invalid))))
    expect_true(all(attr(confint(invalid),"inference_status")=="unavailable_nonregular_fit"))
  }
  expect_error(confint(competing,interval="wald"),class="coherentSurv_error")
})

test_that("targeted confint preserves Gaussian guards and projection scope", {
  d<-data.frame(id=1:12,path=rep(c(1,1,1,2,2,2),2),recorded=1)
  fit<-cs_fit(~1,d,time_unit="year",targeting="finite_path",targeting_spec=list(
    path="path",recorded="recorded",features=matrix(c(0,1),ncol=1,dimnames=list(NULL,"risk")),
    folds=rep(1:2,each=6),horizon=1,complete_followup=TRUE))
  set.seed(19);before<-.Random.seed
  p<-predict(fit,draws=100,seed=3,level=.9)
  result<-confint(fit,parm="risk",draws=100,seed=3,level=.9)
  expect_identical(.Random.seed,before)
  expect_equal(as.numeric(result),c(p$lower,p$upper))
  expect_identical(unname(attr(result,"inference_status")),p$inference_status)
  fit$convergence<-list(valid=FALSE,code="zero_event_boundary")
  explicit<-confint(fit,interval="logit_gaussian",draws=100)
  expect_true(all(is.na(explicit)))
  expect_match(attr(explicit,"inference_status"),"zero_event_boundary")
  automatic<-confint(fit,draws=100)
  expect_true(all(is.finite(automatic)))
  expect_match(attr(automatic,"inference_scope"),"does not inherit")
  projection<-confint(fit,interval="projection",draws=100)
  expect_equal(as.numeric(projection),as.numeric(automatic))
  expect_match(attr(projection,"inference_scope"),"used throughout")
  expect_error(confint(fit,landmark=.5),class="coherentSurv_error")
})

test_that("causal confint retains coordinate-specific guards and contrast methods", {
  # Known EIF arrays isolate the generic contract from nuisance fitting.
  e0<-cbind(risk=c(-.1,.1,-.1,.1),rmst=c(-.2,0,0,.2))
  e1<-cbind(risk=c(-.2,.2,-.2,.2),rmst=c(-.1,0,0,.1))
  fit<-structure(list(n=4,target_names=c("risk","rmst"),
    arms=list("0"=list(estimate=c(risk=.3,rmst=.8),EIF=e0,regular=c(TRUE,FALSE)),
      "1"=list(estimate=c(risk=.2,rmst=.9),EIF=e1,regular=c(TRUE,TRUE)))),
    class=c("coherentSurv_causal","coherentSurv"))
  for(contrast in c("difference","ratio","log_ratio")) {
    expected<-cs_causal_contrast(fit,contrast=contrast,level=.9)
    actual<-confint(fit,contrast=contrast,level=.9)
    expect_equal(as.numeric(actual),c(expected$lower,expected$upper))
    expect_true(all(is.finite(actual["risk",])))
    expect_true(all(is.na(actual["rmst",])))
    expect_identical(unname(attr(actual,"inference_status")),expected$inference_status)
  }
  expected<-cs_causal_contrast(fit,interval="multiplier",draws=100,seed=31)
  actual<-confint(fit,parm=1,interval="multiplier",draws=100,seed=31)
  expect_equal(as.numeric(actual),c(expected$lower[1],expected$upper[1]))
  expect_match(attr(actual,"inference_scope"),"regular fitted contrast coordinates")
  expect_error(confint(fit,unused=TRUE),class="coherentSurv_error")
})
