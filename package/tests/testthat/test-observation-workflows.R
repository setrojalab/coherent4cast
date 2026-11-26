test_that("duplicate legal paths integrate to the exact-history likelihood", {
  data(sim_perinatal)
  d<-sim_perinatal;other<-d;other$id<-other$id+160
  fit<-cs_fit(~x,d,time_unit="months",breaks=c(0,1,Inf),event_types=c("1","2"),control=cs_control(phi_fixed=.4))
  mix<-cs_fit(~x,rbind(d,other),time_unit="months",breaks=c(0,1,Inf),event_types=c("1","2"),
    observation="coarsened",coarsening=list(owner=rep(1:160,2),log_weights=rep(log(.5),320)),control=cs_control(phi_fixed=.4))
  expect_equal(mix$n,160);expect_equal(mix$npaths,320)
  expect_equal(mix$loglik,fit$loglik,tolerance=1e-7)
  expect_equal(coef(mix),coef(fit),tolerance=1e-5)
  expect_equal(vcov(mix),vcov(fit),tolerance=1e-5)
  expect_error(predict(mix,landmark=0,horizon=1),class="coherentSurv_error")
  weight<-cs_fit(~x,d,time_unit="months",weights=rep(2,160),control=cs_control(phi_fixed=.4))
  normal<-cs_fit(~x,d,time_unit="months",control=cs_control(phi_fixed=.4))
  expect_equal(coef(weight),coef(normal),tolerance=1e-9)
  expect_equal(vcov(weight),vcov(normal),tolerance=1e-9)
})

test_that("whole-unit selection and omission preserve complete histories", {
  data(sim_perinatal)
  chosen<-suppressWarnings(cs_fit(~x,sim_perinatal,time_unit="months",subset=id<=80,control=cs_control(phi_fixed=.4)))
  expect_equal(chosen$n,80)
  missing<-sim_perinatal;missing$x[1]<-NA
  clean<-cs_fit(~x,missing,time_unit="months",na.action=stats::na.omit,control=cs_control(phi_fixed=.4))
  expect_equal(clean$n,159);expect_identical(clean$audit$omitted_units,"1")
  expect_false(any(clean$data$id==1))
  raw<-cs_fit(~x,sim_perinatal,time_unit="months",control=cs_control(phi_fixed=.4))
  expect_equal(validate(raw,c(1,4),c(1,0),3,predicted=c(.8,.2))$n,2)
  for(t in c("risk","survival","rmst"))expect_equal(nrow(predict(raw,sim_perinatal[sim_perinatal$id==1,],landmark=0,horizon=.2,type=t)),1)
})

test_that("nonrepeatable event eligibility removes inactive coefficients", {
  d<-cs_simulate(100,rates=c(.6,.15),beta=c(0,0),history=c(0,0),seed=500)
  for(rows in split(seq_len(nrow(d)),d$id)) {
    jumps<-rows[d$event_type[rows]==1]
    if(length(jumps)>1)d$event_type[jumps[-1]]<-0
  }
  fit<-cs_fit(~1,d,time_unit="years",nonrecurrent="1",control=cs_control(phi_fixed=0))
  expect_true("event1:count_1" %in% fit$audit$structurally_inactive_coefficients)
  expect_false("event1:count_1" %in% names(coef(fit)))
  pred<-predict(fit,d[d$id==1,],landmark=0,horizon=5,type="burden")
  expect_true(all(pred$estimate>=0&pred$estimate<=1+1e-10))
})

test_that("invalid model specifications and nonregular fits remain visible", {
  d<-data.frame(id=1:8,start=0,stop=1:8,event_type=0,terminal=c(1,1,0,1,0,1,0,0))
  fit<-function(...)cs_fit(~1,d,time_unit="years",...)
  expect_error(fit(targeting="continuous"),class="coherentSurv_error")
  expect_error(cs_fit(stop~id,d,time_unit="years"),class="coherentSurv_error")
  expect_error(fit(breaks=c(1,Inf)),class="coherentSurv_error")
  expect_error(fit(subset=TRUE),class="coherentSurv_error")
  expect_error(cs_fit(~unknown,d,time_unit="years"),class="coherentSurv_error")
  expect_error(fit(nonrecurrent="1"),class="coherentSurv_error")
  auto<-suppressWarnings(fit(observation="coarsened"))
  exact<-suppressWarnings(fit())
  expect_equal(auto$loglik,exact$loglik,tolerance=1e-8)
  expect_false(is.null(auto$audit$automatic_fibre))
  expect_error(fit(observation="coarsened",coarsening=list(owner=1,log_weights=0)),class="coherentSurv_error")
  expect_error(fit(observation="coarsened",coarsening=list(owner=1:8,log_weights=rep(0,8)),weights=1),class="coherentSurv_error")
  expect_warning(bad<-fit(control=cs_control(maxit=1)),"nonregular")
  expect_false(bad$convergence$valid)
  expect_equal(coherentSurv:::.cs_inverse(matrix(0,2,2)),matrix(0,2,2))
  expect_error(cs_regularize(bad,maxit=NA),class="coherentSurv_error")
  expect_error(cs_simulate(n=0),class="coherentSurv_error")
  expect_error(cs_simulate(n=1,history=c(Inf,0,0)),class="coherentSurv_error")
  expect_error(cs_build_history(d),class="coherentSurv_error")
  tr<-cs_as_trajectory(d,time_unit="years")
  expect_error(cs_build_history(tr,cap=-1),class="coherentSurv_error")
  expect_error(cs_build_history(tr,decay=-1),class="coherentSurv_error")
  expect_error(cs_as_trajectory(d),class="coherentSurv_error")
  expect_error(cs_as_trajectory(d[FALSE,],time_unit="years"),class="coherentSurv_error")
  expect_error(cs_control(phi_fixed=-.1),class="coherentSurv_error")
})

test_that("prediction rejects unsupported conditioning and numerical overflow", {
  data(sim_perinatal)
  fit<-cs_fit(~x,sim_perinatal,time_unit="months",event_types=c("1","2"))
  d<-sim_perinatal[sim_perinatal$id==1,]
  expect_error(predict(fit,horizon=1),class="coherentSurv_error")
  expect_error(predict(fit,landmark=0,horizon=0),class="coherentSurv_error")
  expect_error(predict(fit,landmark=0,horizon=1,quadrature=1),class="coherentSurv_error")
  short<-fit;short$control$max_states<-2
  expect_error(predict(short,d,landmark=0,horizon=1),class="coherentSurv_error")
  overflow<-fit;overflow$coefficients[1]<-1000
  expect_error(predict(overflow,d,landmark=0,horizon=1),class="coherentSurv_error")
  expect_error(coherentSurv:::.cs_gamma_nodes(0,1,5),class="coherentSurv_error")
  gap<-data.frame(id=c(1,1),start=c(0,2),stop=c(1,3),event_type=0,terminal=0,x=0)
  expect_error(predict(fit,gap,landmark=.5,horizon=1),class="coherentSurv_error")
  paths<-data.frame(id=1,start=0,stop=1,event_type=0,terminal=0,x=0)
  expect_error(cs_predict_paths(NULL,paths,1,0,1,1),class="coherentSurv_error")
  expect_error(cs_predict_paths(fit,paths,1,0,0,1),class="coherentSurv_error")
  expect_error(cs_predict_paths(fit,paths,2,0,1,1),class="coherentSurv_error")
  expect_error(cs_predict_paths(fit,paths,1,NA,1,1),class="coherentSurv_error")
  delayed<-paths;delayed$start<-.1
  expect_error(cs_predict_paths(fit,delayed,1,0,1,1),class="coherentSurv_error")
})

test_that("clinical source errors and validation degeneracy are explicit", {
  expect_error(cs_readmission(data.frame()),class="coherentSurv_error")
  expect_error(cs_readmission(data.frame(id=1)),class="coherentSurv_error")
  expect_error(coherentSurv:::.cs_source_data("x","not_an_installed_package_coherentSurv",NULL),class="coherentSurv_error")
  small<-data.frame(id=1,t.start=0,t.stop=2,event=0,death=0,sex="Female",chemo="Treated",dukes="C",charlson="0")
  bad<-small;bad$id<-NA
  expect_error(cs_readmission(bad),class="coherentSurv_error")
  gap<-rbind(small,small);gap$t.start<-c(0,3);gap$t.stop<-c(2,4)
  expect_error(cs_readmission(gap),class="coherentSurv_error")
  expect_equal(coherentSurv:::.cs_legendre_unit(1),list(x=.5,w=1))
  expect_error(validate(numeric(),numeric(),numeric(),1),class="coherentSurv_error")
  expect_error(validate(c(.1,.2),1:2,c(1,0),1,threshold=NA),class="coherentSurv_error")
  expect_error(validate(c(.1,.2),1:2,c(1,0),1,strata=1),class="coherentSurv_error")
  expect_error(validate(c(.1,.2),1:2,c(1,0),1,censor_survival=list(case=1,horizon=1)),class="coherentSurv_error")
  grouped<-validate(c(.1,.2,.7,.9),1:4,c(1,0,1,0),2,strata=c(1,1,2,2),censor_survival=list(case=rep(.9,4),horizon=rep(.8,4)))
  expect_equal(length(grouped),2)
  expect_equal(bounds(c(.4,.6),c(0,1),kappa=1,recording_fraction=.8)$lower,.6,tolerance=1e-12)
  expect_equal(length(conformal(c(.1,.2),c(.15,.3))$candidates$included),2)
})
