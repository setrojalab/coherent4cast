p4_history <- function() data.frame(id = c(1,1,2,2), start = c(0,1,0,1),
  stop = c(1,2,1,2), event_type = c("1","0","1","0"),
  terminal = c(0,1,0,0), x = c(0,0,1,1))
p4_trajectory <- function(d, ...) cs_as_trajectory(d, time_unit = "years", event_types = "1", ...)
p4_bad <- function(expr, code) expect_error(expr, regexp = code, class = "coherentSurv_error")

test_that("missing unit identifiers are rejected without dropping rows", {
  d <- p4_history(); d$id[2] <- NA; p4_bad(p4_trajectory(d), "unit")
})
test_that("empty unit identifiers cannot join otherwise unrelated histories", {
  d <- p4_history(); d$id <- as.character(d$id); d$id[1] <- ""; p4_bad(p4_trajectory(d), "unit")
})
test_that("a missing endpoint column is rejected", {
  d <- p4_history(); d$stop <- NULL; p4_bad(p4_trajectory(d), "columns")
})
test_that("numeric times stored as strings require an explicit conversion", {
  d <- p4_history(); d$start <- as.character(d$start); p4_bad(p4_trajectory(d), "time")
})
test_that("missing times do not silently disappear", {
  d <- p4_history(); d$stop[1] <- NA_real_; p4_bad(p4_trajectory(d), "time")
})
test_that("infinite observed followup is rejected", {
  d <- p4_history(); d$stop[2] <- Inf; p4_bad(p4_trajectory(d), "time")
})
test_that("negative starts are invalid", {
  d <- p4_history(); d$start[1] <- -1; p4_bad(p4_trajectory(d), "interval")
})
test_that("zero duration events require a coarsening model", {
  d <- p4_history(); d$stop[1] <- 0; p4_bad(p4_trajectory(d), "interval")
})
test_that("reversed intervals are rejected", {
  d <- p4_history(); d$stop[2] <- 0.5; p4_bad(p4_trajectory(d), "interval")
})
test_that("terminal coding must be binary", {
  d <- p4_history(); d$terminal[2] <- 2; p4_bad(p4_trajectory(d), "terminal")
})
test_that("missing terminal indicators are rejected", {
  d <- p4_history(); d$terminal[2] <- NA; p4_bad(p4_trajectory(d), "terminal")
})
test_that("a terminal jump cannot share an exact row with recurrence", {
  d <- p4_history(); d$terminal[1] <- 1; p4_bad(p4_trajectory(d), "tie")
})
test_that("missing event labels are rejected", {
  d <- p4_history(); d$event_type[1] <- NA; p4_bad(p4_trajectory(d), "event_type")
})
test_that("an undeclared event label is not silently pooled", {
  d <- p4_history(); d$event_type[1] <- "sepsis"; p4_bad(p4_trajectory(d), "event_type")
})
test_that("duplicate declared marks are invalid", {
  p4_bad(cs_as_trajectory(p4_history(), time_unit="years", event_types=c("1","1")), "event_types")
})
test_that("overlapping intervals cannot double-count time at risk", {
  d <- p4_history(); d$start[2] <- .5; p4_bad(p4_trajectory(d), "overlap")
})
test_that("out-of-order histories must be sorted explicitly", {
  d <- p4_history()[c(2,1,3,4),]; p4_bad(p4_trajectory(d), "ordering")
})
test_that("interleaved patients are not silently reordered", {
  d <- p4_history()[c(1,3,2,4),]; p4_bad(p4_trajectory(d), "ordering")
})
test_that("records after a terminal event are rejected", {
  d <- p4_history(); d$event_type[1] <- "0"; d$terminal <- c(1,0,0,0)
  p4_bad(p4_trajectory(d), "post_terminal")
})
test_that("nonbinary observation masks cannot become TRUE", {
  d <- p4_history(); d$recorded <- c(1,2,1,1)
  p4_bad(p4_trajectory(d, observed="recorded"), "observation")
})
test_that("exact fitting rejects an unobserved interval", {
  d <- p4_history(); d$recorded <- c(1,0,1,1)
  p4_bad(cs_fit(~x,d,time_unit="years",observed="recorded"), "unrecorded_history")
})
test_that("unrecorded gaps are recorded and rejected by exact fitting", {
  d <- p4_history(); d$start[2] <- 1.1
  expect_identical(attr(p4_trajectory(d),"cs_schema")$gaps, 2L)
  p4_bad(cs_fit(~x,d,time_unit="years"), "unrecorded_history")
})
test_that("delayed entry is not treated as survival from origin", {
  d <- p4_history(); d$start[1] <- .2
  p4_bad(cs_fit(~x,d,time_unit="years"), "left_truncation")
})
test_that("future endpoint columns cannot enter the covariate formula", {
  p4_bad(cs_fit(~stop,p4_history(),time_unit="years"), "formula_leakage")
})
test_that("a row subset cannot truncate a patient's history silently", {
  p4_bad(cs_fit(~x,p4_history(),time_unit="years",subset=start==0), "subset")
})
test_that("default fitting reports missing covariates", {
  d <- p4_history(); d$x[1] <- NA
  p4_bad(cs_fit(~x,d,time_unit="years"), "missing")
})
test_that("zero events are not fitted as an identified likelihood", {
  d <- p4_history(); d$event_type <- "0"; d$terminal <- 0
  p4_bad(cs_fit(~x,d,time_unit="years",event_types="1"), "no_events")
})
test_that("case weights cannot be negative", {
  p4_bad(cs_fit(~x,p4_history(),time_unit="years",weights=c(1,-1)), "weights")
})
test_that("case weights are indexed by patients rather than rows", {
  p4_bad(cs_fit(~x,p4_history(),time_unit="years",weights=rep(1,4)), "weights")
})
test_that("causal interventions are not silently interpreted as prognosis", {
  p4_bad(cs_fit(~x,p4_history(),time_unit="years",intervention=list(x=0)), "causal_scope")
})
test_that("recurrent records cannot enter a nonrepeatable mark model", {
  d <- p4_history(); d$event_type[4] <- "1"
  p4_bad(cs_fit(~x,d,time_unit="years",nonrecurrent="1"), "eligibility")
})
test_that("a missing control value produces a structured error", {
  p4_bad(cs_control(maxit=NA_integer_), "control")
})
test_that("vector controls are rejected rather than partially used", {
  p4_bad(cs_control(tolerance=c(1e-8,1e-6)), "control")
})
test_that("fractional state caps are not silently rounded", {
  p4_bad(cs_control(history_cap=1.5), "control")
})
test_that("reversed optimizer bounds are rejected", {
  p4_bad(cs_control(coefficient_bounds=c(8,-12)), "control")
})
test_that("negative path event counts cannot define gamma posteriors", {
  p4_bad(cs_filter_paths(0,0,-1,1,.5), "path_filter")
})
test_that("negative integrated hazards cannot define gamma posteriors", {
  p4_bad(cs_filter_paths(0,0,1,-1,.5), "path_filter")
})
test_that("noninteger path owners are rejected", {
  p4_bad(cs_filter_paths(0,0,1,1,.5,owner=.5), "path_filter")
})
test_that("path weights remain stable when likelihoods underflow directly", {
  z <- cs_filter_paths(c(-10000,-10001),c(0,0),c(0,1),c(1,2),.5)
  expect_equal(sum(z$weight),1,tolerance=1e-15)
  expect_equal(z$weight[1],1/(1+exp(-1)),tolerance=1e-15)
  expect_equal(z$shape,c(2,3))
})
test_that("predictable counts do not contain the current or future jump", {
  h <- cs_build_history(p4_trajectory(p4_history()),decay=1)
  expect_equal(unname(h$counts[,1]),c(0,1,0,1))
  expect_true(all(h$recency[c(1,3),1]==0))
  expect_true(all(is.na(h$gap[c(1,3),1])))
})
test_that("factor event labels preserve the explicit event dictionary", {
  d <- p4_history(); d$event_type <- factor(d$event_type,levels=c("0","1","unused"))
  z <- p4_trajectory(d); expect_equal(z$.cs_mark,c(1,0,1,0)); expect_equal(nrow(z),4)
})
test_that("unknown outcomes beyond a landmark cannot alter its forecast", {
  train <- cs_simulate(n=35,horizon=2,seed=4242)
  fit <- suppressWarnings(cs_fit(~x,train,time_unit="years",control=cs_control(phi_fixed=.4,maxit=100)))
  a <- data.frame(id=c(1,1),start=c(0,1),stop=c(1,2),event_type=c("0","1"),terminal=0,x=c(0,1))
  b <- a; b$x[2] <- 100; b$event_type[2] <- "0"; b$terminal[2] <- 1
  pa <- predict(fit,a,landmark=.5,horizon=.2,type="all")
  pb <- predict(fit,b,landmark=.5,horizon=.2,type="all")
  expect_equal(pa$estimate,pb$estimate,tolerance=1e-12)
})
