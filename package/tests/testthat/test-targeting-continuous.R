ct_fixture_fit <- local({
  saved <- NULL
  function() {
    if(is.null(saved)) {
      d<-read.csv(test_path("fixtures","continuous_target_subjects.csv"))
      saved<<-recurcast(~x,d,time_unit="year",event_types=c("recurrent1","recurrent2"),
        targeting="continuous",targeting_spec=list(horizon=1,censor="censored",strata="x",folds="fold",
          grid=64,working="pooled",censor_rates=rep(c(.35,.55),5)))
    }
    saved
  }
})
test_that("continuous primary API reproduces the V5 F2 martingale update", {
  fit<-ct_fixture_fit()
  expected<-read.csv(test_path("fixtures","continuous_target_expected.csv"))
  expect_s3_class(fit,"recurcast_targeted")
  expect_true(fit$convergence$valid)
  expect_equal(unname(fit$estimate),expected$estimate,tolerance=1e-11)
  expect_equal(unname(fit$se),expected$se,tolerance=1e-11)
  expect_equal(unname(fit$EIF),unname(as.matrix(read.csv(test_path("fixtures","continuous_target_eif.csv")))),tolerance=1e-10)
  expect_equal(as.numeric(fit$rates),read.csv(test_path("fixtures","continuous_target_rates.csv"))$rate,tolerance=1e-10)
  expect_lt(max(abs(fit$EIF_mean)),1e-7/sqrt(fit$n))
  expect_lt(fit$interpolation_error,1e-6)
  expect_gt(max(abs(fit$estimate-fit$initial)),1e-4)
  for(f in seq_along(fit$training_subjects))expect_false(any(fit$folds[fit$training_subjects[[f]]]==f))
  expect_true(all(fit$rates>=0))
  expect_equal(fit$covariance,cov(fit$EIF)/fit$n,tolerance=1e-14)
  expect_match(paste(fit$assumptions,collapse=" "),"right censoring")
})
test_that("targeted primary prediction automatically supplies finite-vector intervals", {
  fit<-ct_fixture_fit();set.seed(77);prior<-.Random.seed
  p<-predict(fit,draws=200,seed=8)
  expect_identical(.Random.seed,prior)
  expect_equal(p$estimate,unname(fit$estimate))
  expect_true(all(p$lower<p$estimate&p$upper>p$estimate))
  expect_identical(unique(p$inference_status),"regular_asymptotic_approximation")
  expect_gt(p$lower[p$target=="risk"],0);expect_lt(p$upper[p$target=="risk"],1)
  expect_gt(p$lower[p$target=="rmst"],0);expect_lt(p$upper[p$target=="rmst"],fit$horizon)
  raw<-predict(fit,interval="multiplier",draws=200,seed=8)
  expect_equal(raw$lower,pmax(0,raw$estimate-attr(raw,"critical_value")*raw$std_error),tolerance=1e-12)
  no<-predict(fit,interval="none");expect_true(all(is.na(no$lower)))
  s<-predict(fit,type="survival",draws=200,seed=8)
  expect_equal(s$estimate,1-fit$estimate["risk"],ignore_attr=TRUE)
  expect_equal(nrow(tidy(fit,draws=200)),4L)
  expect_equal(glance(fit)$n,500L)
  expect_equal(coef(fit),fit$estimate)
  expect_equal(summary(fit)$estimates$EIF_mean,fit$EIF_mean,ignore_attr=TRUE)
})
test_that("conditional forecasts update one fold-mixture law without using future marks", {
  fit<-ct_fixture_fit();d<-fit$data
  subject<-d[d$id==1,,drop=FALSE]
  at_zero<-predict(fit,newdata=subject,landmark=0,horizon=1,interval="none")
  x<-fit$subject_stratum[1];groups<-(seq_len(fit$nfolds)-1)*length(fit$strata)+x
  expected<-as.vector(fit$values[,1,1,groups]%*%fit$fold_prior[x,])
  expect_equal(at_zero$estimate,expected,tolerance=1e-11)
  candidate<-unique(d$id)[vapply(split(d,d$id),function(z)max(z$stop)>=.2&&!any(z$terminal==1&z$stop<=.2),logical(1))][1]
  subject<-d[d$id==candidate,,drop=FALSE]
  p<-predict(fit,newdata=subject,landmark=.2,horizon=.4,interval="none")
  w<-attr(p,"fold_posterior")[[1]]
  expect_equal(sum(w),1,tolerance=1e-14)
  expect_true(all(w>=0))
  expect_true(all(is.na(p$std_error)))
  later<-subject$stop>.2;subject$event_type[later]<-"none";subject$terminal[later]<-0
  changed<-predict(fit,newdata=subject,landmark=.2,horizon=.4,interval="none")
  expect_equal(p$estimate,changed$estimate,tolerance=1e-12)
  risk<-predict(fit,newdata=subject,landmark=.2,horizon=.4,type="risk",interval="none")
  surv<-predict(fit,newdata=subject,landmark=.2,horizon=.4,type="survival",interval="none")
  expect_equal(risk$estimate+surv$estimate,1,tolerance=1e-13)
  expect_error(predict(fit,newdata=subject,landmark=.2,horizon=.4,interval="multiplier"),class="recurcast_error")
  expect_error(predict(fit,newdata=subject,landmark=.2,horizon=1),class="recurcast_error")
})
test_that("finite-path targeting is exposed through the same fit and predict interface", {
  d<-read.csv(test_path("fixtures","target_subjects.csv"));d$id<-seq_len(nrow(d))
  features<-as.matrix(read.csv(test_path("fixtures","target_features.csv")))
  z<-recurcast(~strata,d,time_unit="window",targeting="finite_path",
     targeting_spec=list(path="path",recorded="recorded",strata="strata",features=features,folds="fold",horizon=2))
  old<-rc_target(d$path,d$recorded,d$strata,features,d$fold)
  expect_s3_class(z,"recurcast_targeted")
  expect_equal(z$estimate,old$estimate,tolerance=1e-12)
  expect_equal(z$EIF,old$EIF,tolerance=1e-12)
  p<-predict(z,draws=200)
  expect_equal(nrow(p),ncol(features));expect_true(all(is.finite(p$lower)))
  expect_error(predict(z,landmark=1,horizon=1),class="recurcast_error")
})
test_that("continuous targeting requires explicit censor and baseline observation assumptions", {
  d<-read.csv(test_path("fixtures","continuous_target_subjects.csv"))
  spec<-list(horizon=1,censor="censored",strata="x",folds="fold",grid=16,working="pooled")
  expect_error(recurcast(~x,d,time_unit="year",targeting="continuous"),class="recurcast_error")
  missing<-spec;missing$censor<-"missing"
  expect_error(recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=missing),class="recurcast_error")
  bad<-d;bad$censored[]<-0
  expect_error(recurcast(~x,bad,time_unit="year",targeting="continuous",targeting_spec=spec),class="recurcast_error")
  bad<-d;ids<-which(duplicated(bad$id));bad$x[ids[1]]<-1-bad$x[ids[1]]
  expect_error(recurcast(~x,bad,time_unit="year",targeting="continuous",targeting_spec=spec),class="recurcast_error")
  singular<-spec;singular$folds<-rep(1,500)
  expect_error(recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=singular),class="recurcast_error")
  short<-spec;short$maxit<-0L
  expect_warning(z<-recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=short),"Targeting score")
  expect_false(z$convergence$valid)
  expect_true(all(is.na(predict(z,interval="logit_gaussian",draws=200)$lower)))
})

test_that("zero-event Gaussian precision is withheld and bounded uncensored fallback is explicit", {
 d<-data.frame(id=1:40,start=0,stop=1,event_type=rep(c("recurrent1","none"),20),terminal=0,censored=0,x=0)
 spec<-list(horizon=1,censor="censored",strata="x",nfolds=2,grid=16,working="pooled")
 expect_warning(fit<-recurcast(~x,d,time_unit="year",targeting="continuous",event_types=c("recurrent1","recurrent2"),targeting_spec=spec),"boundary")
 expect_true(any(fit$zero_event_boundary))
 expect_false(fit$convergence$valid)
 p<-predict(fit,draws=100)
 expect_true(all(p$lower[p$target!="rmst"]==0))
 expect_lte(p$lower[p$target=="rmst"],1)
 expect_match(p$inference_status[1],"unbounded_burden_outer_support")
 expect_true(all(is.na(predict(fit,interval="logit_gaussian",draws=100)$lower)))
 spec$complete_followup<-TRUE
 expect_warning(complete<-recurcast(~x,d,time_unit="year",targeting="continuous",event_types=c("recurrent1","recurrent2"),targeting_spec=spec),"boundary")
 bounded<-predict(complete,draws=100)
 expect_identical(unique(bounded$interval),"bounded_coarsening_projection")
 risk<-bounded[bounded$target=="risk",]
 expect_equal(risk$lower,0)
 expect_gt(risk$upper,0)
 expect_equal(risk$upper,sqrt(log(4/.05)/(2*40)),tolerance=1e-13)
 expect_true(all(is.infinite(bounded$upper[grepl("^burden_",bounded$target)])))
 expect_true(all(is.finite(bounded$upper[bounded$target%in%c("risk","rmst")])))
 expect_equal(complete$censor_rates,rep(0,length(complete$censor_rates)))
})
