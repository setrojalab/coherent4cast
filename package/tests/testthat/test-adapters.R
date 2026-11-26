p4_colon_rows <- function(day=2,death=2,recurrence=1,dead=1,id=1) data.frame(
  id=id,etype=c(1,2),time=c(day,death),status=c(recurrence,dead),
  rx=factor("Obs",levels=c("Obs","Lev","Lev+5FU")),sex=1,age=60,node4=0,nodes=2)

test_that("canonical colon retains all endpoint records and five same-day pairs", {
  skip_if_not_installed("survival")
  a <- cs_colon()
  expect_equal(a$metadata$patients,929)
  expect_equal(a$metadata$recurrences,468)
  expect_equal(a$metadata$source_deaths,452)
  expect_equal(a$metadata$same_day_ids,c(125,277,324,365,670))
  expect_equal(a$metadata$paths,7447)
  expect_length(unique(a$coarsening$owner),929)
  expect_identical(a$nonrecurrent,"1")
})
test_that("same-day triangle has its analytic area and centroid", {
  a <- cs_colon(p4_colon_rows(),nodes=4)
  w <- exp(a$coarsening$log_weights)
  r <- a$data$stop[a$data$event_type=="1"]*365.25
  d <- a$data$stop[a$data$terminal==1]*365.25
  expect_true(all(r<d))
  expect_equal(sum(w),.5/365.25^2,tolerance=1e-14)
  expect_equal(sum(w*r)/sum(w),4/3,tolerance=1e-13)
  expect_equal(sum(w*d)/sum(w),5/3,tolerance=1e-13)
})
test_that("quadrature splits event cells at baseline discontinuities", {
  a <- cs_colon(p4_colon_rows(day=366,death=400),nodes=2)
  r <- a$data$stop[a$data$event_type=="1"]*365.25
  expect_true(any(r<365.25));expect_true(any(r>365.25))
  expect_equal(a$metadata$paths,8)
  expect_equal(sum(exp(a$coarsening$log_weights)),1/365.25^2,tolerance=1e-14)
})
test_that("censoring endpoints retain unit integration mass", {
  a <- cs_colon(p4_colon_rows(day=3,death=3,recurrence=0,dead=0))
  expect_equal(a$coarsening$log_weights,0)
  expect_equal(a$data$stop,3/365.25)
  expect_equal(a$data$terminal,0)
})
test_that("early recurrence surveillance truncation is disclosed", {
  a <- cs_colon(p4_colon_rows(day=3,death=5,recurrence=0,dead=1))
  expect_equal(a$metadata$surveillance_truncated_ids,1)
  expect_equal(a$data$stop,3/365.25); expect_equal(a$data$terminal,0)
})
test_that("recurrence after death is rejected", {
  expect_error(cs_colon(p4_colon_rows(day=5,death=3)),class="coherentSurv_error",regexp="adapter_order")
})
test_that("duplicate endpoint rows are rejected", {
  d <- p4_colon_rows();d$etype <- c(1,1)
  expect_error(cs_colon(d),class="coherentSurv_error",regexp="adapter_pair")
})
test_that("inconsistent paired covariates are not silently reconciled", {
  d <- p4_colon_rows();d$age[2] <- 61
  expect_error(cs_colon(d),class="coherentSurv_error",regexp="adapter_covariates")
})
test_that("continuous nodes use the declared training median and missing indicator", {
  d <- p4_colon_rows();d$nodes <- NA_real_
  a <- cs_colon(d,covariate="nodes",node_median=7)
  expect_true(all(a$data$nodes_per_5==7/5));expect_true(all(a$data$nodes_missing==1))
  expect_error(cs_colon(d,covariate="nodes"),class="coherentSurv_error",regexp="node_median")
})
test_that("node4 is retained even when count coding differs", {
  d <- p4_colon_rows();d$nodes <- 10;d$node4 <- 0
  expect_true(all(cs_colon(d)$data$node4==0))
})
test_that("invalid quadrature sizes cannot silently change the observation model", {
  for (n in list(NA_real_,0,1.5,c(2,3),65))
    expect_error(cs_colon(p4_colon_rows(),nodes=n),class="coherentSurv_error",regexp="quadrature")
})
test_that("canonical readmission distinguishes intervals from recurrences", {
  skip_if_not_installed("frailtypack")
  a <- cs_readmission()
  expect_equal(a$metadata$patients,403);expect_equal(nrow(a$data),861)
  expect_equal(a$metadata$recurrences,458);expect_equal(a$metadata$terminal_events,109)
  b <- cs_readmission(time_unit="days")
  expect_equal(b$data$stop/365.25,a$data$stop,tolerance=1e-14)
})
test_that("later Charlson changes do not leak into baseline covariates", {
  d <- data.frame(id=1,t.start=c(0,1),t.stop=c(1,2),event=c(1,0),death=0,
    sex="Female",chemo="Treated",dukes="C",charlson=c("0","3"))
  a <- cs_readmission(d)
  expect_equal(a$data$charlson_3,c(0,0));expect_equal(a$data$recorded_charlson,c("0","3"))
})
test_that("unknown clinical categories are not silently coded as the reference level", {
  d <- data.frame(id=1,t.start=0,t.stop=2,event=0,death=0,
    sex="Female",chemo="Treated",dukes="C",charlson="0")
  for (v in c("sex","chemo","dukes","charlson")) {
    bad <- d;bad[[v]] <- "unknown"
    expect_error(cs_readmission(bad),class="coherentSurv_error",regexp="adapter_covariates")
  }
})
test_that("path forecasts equal direct averaging and reject future data", {
  train <- cs_simulate(n=45,horizon=2,seed=177)
  fit <- suppressWarnings(cs_fit(~x,train,time_unit="years",control=cs_control(phi_fixed=.4,maxit=100)))
  paths <- data.frame(id=c(1,2,2),start=c(0,0,.4),stop=c(1,.4,1),event_type=c("0","1","0"),terminal=0,x=0)
  a <- cs_predict_paths(fit,paths,owner=c(1,1),log_weights=log(c(.4,.6)),landmark=1,horizon=.25,type="all")
  q <- attr(a,"path_weights")$posterior_weight
  exact <- predict(fit,paths,landmark=1,horizon=.25,type="all")
  for (target in a$target) {
    z <- exact[exact$target==target,]
    expect_equal(a$estimate[a$target==target],sum(q*z$estimate),tolerance=1e-12)
  }
  expect_equal(sum(q),1,tolerance=1e-15)
  future <- paths; future$stop[c(1,3)] <- 1.1
  expect_error(cs_predict_paths(fit,future,c(1,1),c(0,0),1,.25),class="coherentSurv_error",regexp="path_future")
  short <- paths;short$stop[c(1,3)] <- .9
  expect_error(cs_predict_paths(fit,short,c(1,1),c(0,0),1,.25),class="coherentSurv_error",regexp="path_history")
})
