fibre_records <- function(lower,upper,marks=seq_along(lower),terminal=rep(0,length(lower)),end=max(upper))
  data.frame(id=1,event_time_lower=lower,event_time_upper=upper,event_type=as.character(marks),terminal=terminal,followup=end)
fibre_mass <- function(f) sum(exp(f$coarsening$log_weights))

test_that("same-day recurrence and death integrate their legal triangle", {
  f<-cs_coarsen_fibre(fibre_records(c(1,1),c(2,2),c("a","0"),c(0,1)),time_unit="years")
  expect_equal(fibre_mass(f),.5,tolerance=1e-13)
  expect_equal(f$tie_audit$canonical_orders_enumerated,1)
  expect_true(all(vapply(split(f$data,f$data$id),function(d)all(diff(d$stop)>0)&&tail(d$terminal,1)==1,logical(1))))
  expect_identical(f$tie_audit$method,"quadrature")
})
test_that("distinct marks with overlapping windows retain both orders", {
  f<-cs_coarsen_fibre(fibre_records(c(0,0),c(1,1),c("a","b")),time_unit="years")
  expect_equal(fibre_mass(f),1,tolerance=1e-13)
  expect_equal(f$tie_audit$canonical_orders_enumerated,2)
  first<-vapply(split(f$data,f$data$id),function(d)d$event_type[1],character(1))
  expect_setequal(first,c("a","b"))
})
test_that("identical same-mark windows are not duplicated factorial times", {
  f<-cs_coarsen_fibre(fibre_records(rep(0,3),rep(1,3),rep("a",3)),time_unit="years")
  expect_equal(fibre_mass(f),1/6,tolerance=1e-13)
  expect_equal(f$tie_audit$canonical_orders_enumerated,1)
})
test_that("same-mark overlapping windows use their union not assignment multiplicity", {
  f<-cs_coarsen_fibre(fibre_records(c(0,1),c(2,3),c("a","a")),time_unit="years")
  expect_equal(fibre_mass(f),3.5,tolerance=1e-13)
  # Independent midpoint union integration over the ordered time square.
  grid<-(seq_len(600)-.5)*3/600
  t1<-matrix(grid,600,600);t2<-t(t1)
  in_union<-t1<t2 & ((t1<=2 & t2>=1) | (t1>=1 & t2<=2))
  oracle<-sum(in_union)*(3/600)^2
  expect_lt(abs(fibre_mass(f)-oracle),.01)
})
test_that("canonical matching allows a wide window before a later narrow window", {
  f<-cs_coarsen_fibre(fibre_records(c(1,0),c(2,3),c("a","a")),time_unit="years")
  expect_equal(fibre_mass(f),2.5,tolerance=1e-13)
  expect_equal(f$tie_audit$canonical_orders_enumerated,2)
})
test_that("mixed point and interval observations preserve their different masses", {
  f<-cs_coarsen_fibre(fibre_records(c(.5,1),c(.5,2),c("a","b")),time_unit="years")
  expect_equal(fibre_mass(f),1,tolerance=1e-13)
  expect_true(all(f$data$stop[f$data$event_type=="a"]==.5))
})
test_that("exact ties require a declared coarsening resolution", {
  d<-data.frame(id=1,start=c(0,1),stop=c(1,1),event_type=c("a","0"),terminal=c(0,1))
  expect_error(cs_coarsen_fibre(d,time_unit="days"),class="coherentSurv_error")
  f<-cs_coarsen_fibre(d,time_unit="days",rounding_unit=1)
  expect_equal(fibre_mass(f),.5,tolerance=1e-13)
  expect_error(cs_coarsen_fibre(d),class="coherentSurv_error",regexp="time_unit")
})
test_that("terminal-before-recurrence observations are not repaired by deleting events", {
  expect_error(cs_coarsen_fibre(fibre_records(c(2,0),c(3,1),c("a","0"),c(0,1)),time_unit="years"),class="coherentSurv_error")
})
test_that("nonrepeatable duplicate events are rejected", {
  expect_error(cs_coarsen_fibre(fibre_records(c(0,1),c(1,2),c("a","a")),
    time_unit="years",nonrecurrent="a"),class="coherentSurv_error",regexp="eligibility")
})
test_that("censor-only subjects retain a unit-mass empty-event path", {
  d<-fibre_records(NA_real_,NA_real_,"0",0,end=3)
  f<-cs_coarsen_fibre(d,time_unit="years",event_types="a")
  expect_equal(fibre_mass(f),1)
  expect_equal(f$data$start,0);expect_equal(f$data$stop,3)
  expect_equal(f$data$event_type,"0")
})
test_that("event at final censoring time produces no zero-length extra row", {
  f<-cs_coarsen_fibre(fibre_records(1,1,"a"),time_unit="years")
  expect_equal(nrow(f$data),1);expect_equal(f$data$stop,1)
})
test_that("quadrature respects baseline discontinuity splitting", {
  f<-cs_coarsen_fibre(fibre_records(0,2,"a"),time_unit="years",nodes=2,breaks=c(0,1,Inf))
  times<-f$data$stop[f$data$event_type=="a"]
  expect_length(times,4);expect_true(any(times<1)&&any(times>1))
  expect_equal(fibre_mass(f),2,tolerance=1e-13)
})
test_that("importance fallback is explicit and reproducible without touching R RNG", {
  d<-fibre_records(rep(0,6),rep(1,6),letters[1:6])
  set.seed(77);before<-.Random.seed
  a<-cs_coarsen_fibre(d,time_unit="years",max_perms=5,draws=600,seed=4242)
  b<-cs_coarsen_fibre(d,time_unit="years",max_perms=5,draws=600,seed=4242)
  expect_identical(.Random.seed,before)
  expect_identical(a$data,b$data);expect_identical(a$coarsening,b$coarsening)
  expect_identical(a$tie_audit$method,"importance_sampling")
  expect_equal(a$tie_audit$attempted_draws,600)
  expect_true(a$metadata$sampled)
  expect_lt(abs(fibre_mass(a)-1),5*a$tie_audit$volume_mcse)
  se<-cs_fibre_mcse(a,rep(-10000,length(a$coarsening$owner)))
  expect_equal(se$relative_likelihood_mcse,a$tie_audit$volume_mcse/fibre_mass(a),tolerance=1e-11)
})
test_that("quadrature path budget also invokes a disclosed sampling fallback", {
  d<-fibre_records(rep(0,5),rep(1,5),rep("a",5))
  f<-cs_coarsen_fibre(d,time_unit="years",nodes=4,max_paths=100,draws=100)
  expect_identical(f$tie_audit$method,"importance_sampling")
  expect_match(f$tie_audit$reason,"max_paths")
})
test_that("quadrature does not fabricate a Monte Carlo standard error", {
  f<-cs_coarsen_fibre(fibre_records(0,1,"a"),time_unit="years")
  e<-cs_fibre_mcse(f,rep(0,length(f$coarsening$owner)))
  expect_true(is.na(e$relative_likelihood_mcse))
  expect_equal(e$log_likelihood,0,tolerance=1e-13)
})
test_that("unknown counts in an unrecorded counting-process gap are rejected", {
  d<-data.frame(id=1,start=c(0,2),stop=c(1,3),event_type=c("a","0"),terminal=0)
  expect_error(cs_coarsen_fibre(d,time_unit="days",rounding_unit=1),class="coherentSurv_error",regexp="fibre_gap")
})
test_that("time-varying covariates are not silently reordered with hidden events", {
  d<-fibre_records(c(0,0),c(1,1),c("a","b"));d$x<-c(0,1)
  expect_error(cs_coarsen_fibre(d,time_unit="days"),class="coherentSurv_error",regexp="fibre_covariates")
})
test_that("custom column names and source identifiers remain traceable", {
  d<-data.frame(person="patient-A",lo=0,hi=1,kind="visit",dead=0,end=2,x=1)
  f<-cs_coarsen_fibre(d,unit="person",lower="lo",upper="hi",event_type="kind",terminal="dead",followup="end",time_unit="days")
  expect_true(all(c("person","start","stop","kind","dead","x")%in%names(f$data)))
  expect_true(all(f$data$.cs_source_unit=="patient-A"))
  expect_true(all(f$data$x==1));expect_equal(unique(f$coarsening$owner),1)
})
test_that("the automatic colon generator matches the validated day-coarsening adapter", {
  skip_if_not_installed("survival")
  raw<-getExportedValue("survival","colon")
  selected<-raw[raw$id%in%c(1,2,125,277,324,365,670),]
  a<-cs_colon(selected,nodes=4)
  b<-cs_colon_fibre(selected,nodes=4)
  expect_equal(a$data$start,b$data$start,tolerance=1e-13)
  expect_equal(a$data$stop,b$data$stop,tolerance=1e-13)
  expect_equal(a$data$event_type,b$data$event_type)
  expect_equal(a$data$terminal,b$data$terminal)
  expect_equal(a$coarsening$owner,b$coarsening$owner)
  expect_equal(a$coarsening$log_weights,b$coarsening$log_weights,tolerance=1e-11)
  expect_equal(b$metadata$same_day_ids,c(125,277,324,365,670))
})
test_that("the main coarsened fit automatically constructs the same observation likelihood", {
  if(!"coarsening_control"%in%names(formals(cs_fit)))skip("main API integration is not installed in the isolated reference library")
  d<-cs_simulate(n=35,horizon=2,rates=c(.45,.2),beta=c(-.1,.1),history=c(.2,.2),seed=7701)
  d$start<-ceiling(d$start*4)/4;d$stop<-ceiling(d$stop*4)/4
  fibre<-cs_coarsen_fibre(d,time_unit="years",rounding_unit=.25,nodes=2,covariates="x")
  ctl<-cs_control(phi_fixed=.4,history_cap=1,maxit=200)
  manual<-suppressWarnings(cs_fit(~x,fibre$data,time_unit="years",observation="coarsened",
    coarsening=fibre$coarsening,control=ctl,frailty_diagnostics=FALSE))
  automatic<-suppressWarnings(cs_fit(~x,d,time_unit="years",observation="coarsened",
    coarsening_control=list(rounding_unit=.25,nodes=2),control=ctl,frailty_diagnostics=FALSE))
  expect_equal(automatic$n,35)
  expect_equal(automatic$loglik,manual$loglik,tolerance=1e-11)
  expect_equal(automatic$coefficients,manual$coefficients,tolerance=1e-11)
  expect_equal(automatic$information,manual$information,tolerance=1e-10)
  expect_equal(automatic$audit$automatic_fibre$metadata$subjects,35)
})
