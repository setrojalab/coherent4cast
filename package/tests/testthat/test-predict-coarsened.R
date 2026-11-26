test_that("automatic coarsened forecasts preserve subjects and match explicit mixtures", {
  training<-cs_simulate(n=60,horizon=2,rates=c(.5,.2),beta=c(-.3,.2),history=c(.2,.2),seed=912)
  fit<-suppressWarnings(cs_fit(~x,training,time_unit="years",control=cs_control(phi_fixed=.4),frailty_diagnostics=FALSE))
  h<-data.frame(id=c("A-17","B-02"),event_time_lower=c(.2,.4),event_time_upper=c(.6,.8),
    event_type="1",terminal=0,followup=1,x=c(0,1))
  f<-cs_coarsen_fibre(h,time_unit="years",event_types="1",nodes=2)
  explicit<-cs_predict_paths(fit,f$data,f$coarsening$owner,f$coarsening$log_weights,1,.3,type="all",quadrature=8)
  automatic<-cs_predict_coarsened(fit,h,1,.3,type="all",quadrature=8,coarsening_control=list(nodes=2))
  expect_setequal(automatic$unit,c("A-17","B-02"))
  expect_equal(automatic$estimate,explicit$estimate,tolerance=1e-12)
  expect_true(all(attr(automatic,"path_weights")$source_unit%in%h$id))
  for(id in h$id)expect_equal(sum(automatic$estimate[automatic$unit==id&automatic$target%in%c("risk","survival")]),1,tolerance=1e-12)
  if("newdata_observation"%in%names(formals(getS3method("predict","coherentSurv")))) {
    routed<-predict(fit,h,landmark=1,horizon=.3,type="all",quadrature=8,newdata_observation="coarsened",coarsening_control=list(nodes=2))
    expect_equal(routed$estimate,automatic$estimate,tolerance=1e-12)
  }
  future<-h;future$event_time_upper[1]<-1.2
  expect_error(cs_predict_coarsened(fit,future,1,.3),class="coherentSurv_error",regexp="path_future")
  future<-h;future$followup<-1.2
  expect_error(cs_predict_coarsened(fit,future,1,.3),class="coherentSurv_error",regexp="path_future")
  dead<-h;dead$terminal[1]<-1
  expect_error(cs_predict_coarsened(fit,dead,1,.3),class="coherentSurv_error",regexp="landmark_terminal")
  expect_error(cs_predict_coarsened(fit,h,1,.3,coarsening_control=list(time_unit="days")),class="coherentSurv_error",regexp="coarsening_control")
  future_schedule<-data.frame(id=h$id,time=1.1,x=c(1,0))
  map<-unique(f$data[,c("id",".cs_source_unit")])
  mapped<-do.call(rbind,lapply(seq_len(nrow(map)),function(i) {
    z<-future_schedule[future_schedule$id==map$.cs_source_unit[i],];z$id<-map$id[i];z
  }))
  manual_evolving<-cs_predict_paths(fit,f$data,f$coarsening$owner,f$coarsening$log_weights,1,.3,type="all",quadrature=8,
    covariate_evolution="locf",newdata_future=mapped)
  auto_evolving<-cs_predict_coarsened(fit,h,1,.3,type="all",quadrature=8,coarsening_control=list(nodes=2),
    covariate_evolution="locf",newdata_future=future_schedule)
  expect_equal(auto_evolving$estimate,manual_evolving$estimate,tolerance=1e-12)
  expect_true(all(auto_evolving$covariate_evolution=="locf"))
})

test_that("coarsened joint covariate forecasts use measured times and exclude future measurements", {
  if(!"covariate_history"%in%names(formals(getS3method("predict","coherentSurv"))))skip("measured-history integration not yet installed")
  training<-cs_simulate(n=50,horizon=2,rates=c(.5,.2),beta=c(-.3,.2),history=c(.2,.2),seed=913)
  fit<-suppressWarnings(cs_fit(~x,training,time_unit="years",control=cs_control(phi_fixed=.4),frailty_diagnostics=FALSE))
  h<-data.frame(id="A-17",event_time_lower=.2,event_time_upper=.6,event_type="1",terminal=0,followup=1,x=0)
  model<-structure(list(models=list(x=list(beta=c(0,.2),D=matrix(.05,1,1),sigma2=.1)),
    variables="x",unit="id",time="time",center=0,scale=1,random_slope=FALSE,
    assumptions="fixed Gaussian working submodel for this numerical test"),class="cs_covariate_model")
  measured<-data.frame(id="A-17",time=c(0,.5),x=c(0,.15))
  a<-cs_predict_coarsened(fit,h,1,.3,type="all",quadrature=4,coarsening_control=list(nodes=2),
    covariate_evolution="joint_submodel",covariate_model=model,covariate_history=measured,
    covariate_draws=8,covariate_grid=2,seed=271)
  future<-rbind(measured,data.frame(id="A-17",time=2,x=999))
  b<-cs_predict_coarsened(fit,h,1,.3,type="all",quadrature=4,coarsening_control=list(nodes=2),
    covariate_evolution="joint_submodel",covariate_model=model,covariate_history=future,
    covariate_draws=8,covariate_grid=2,seed=271)
  expect_equal(a$estimate,b$estimate,tolerance=1e-12)
  expect_true(all(a$covariate_evolution=="joint_submodel"))
  expect_error(cs_predict_coarsened(fit,h,1,.3,covariate_evolution="joint_submodel",covariate_model=model),
    class="coherentSurv_error",regexp="covariate_history")
  duplicate<-rbind(measured,measured[1,])
  expect_error(cs_predict_coarsened(fit,h,1,.3,covariate_evolution="joint_submodel",covariate_model=model,covariate_history=duplicate),
    class="coherentSurv_error",regexp="one measured")
})
