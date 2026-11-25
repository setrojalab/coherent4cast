test_that("alternative-frailty analytic scores and Louis information match differences", {
  X<-cbind(1,c(-.3,.2,.5,-.4,.1,.8));id<-rep(1:3,each=2);ex<-c(.4,.5,.2,.7,.3,.6);ev<-c(1,0,0,1,1,1)
  for(case in c("lognormal","discrete2","discrete3")) {
    fam<-if(case=="lognormal")"lognormal"else "discrete";points<-if(case=="discrete3")3L else 2L
    th<-c(-.7,.2,if(fam=="lognormal")log(.4)else if(points==2)c(.3,-1.1)else c(.2,-.3,-1.2,.4))
    f<-function(t)rc_frailty_loglik(t,id,ex,X,ev,frailty=fam,quadrature=40,points=points)
    z<-f(th);eps<-1e-5
    grad<-vapply(seq_along(th),function(j){up<-dn<-th;up[j]<-up[j]+eps;dn[j]<-dn[j]-eps;(sum(f(up)$subject_loglik)-sum(f(dn)$subject_loglik))/(2*eps)},numeric(1))
    hess<-sapply(seq_along(th),function(j){up<-dn<-th;up[j]<-up[j]+eps;dn[j]<-dn[j]-eps;(colSums(f(up)$scores)-colSums(f(dn)$scores))/(2*eps)})
    expect_equal(colSums(z$scores),grad,tolerance=2e-7)
    expect_equal(z$information,-hess,tolerance=2e-6)
    expect_equal(rowSums(z$posterior_weights),rep(1,3),tolerance=1e-13)
    expect_equal(apply(z$information_individual,c(1,2),sum),z$information,tolerance=1e-13)
    mixed<-function(t)rc_coarsen(f(t),owner=c(1,1,2),log_weights=log(c(.3,.7,1)))
    hess2<-sapply(seq_along(th),function(j){up<-dn<-th;up[j]<-up[j]+eps;dn[j]<-dn[j]-eps;(colSums(mixed(up)$scores)-colSums(mixed(dn)$scores))/(2*eps)})
    expect_equal(mixed(th)$information,-hess2,tolerance=2e-6)
  }
})
test_that("lognormal quadrature matches direct integration and discrete scale is identified", {
  th<-c(log(.4),log(.7));z<-rc_frailty_loglik(th,1,2,matrix(1,1,1),1,"lognormal",quadrature=80)
  v<-log1p(.7);truth<-log(.4)+log(stats::integrate(function(t)exp(t)*exp(-.8*exp(t))*stats::dnorm(t,-v/2,sqrt(v)),-12,12,rel.tol=1e-12)$value)
  expect_equal(z$subject_loglik,truth,tolerance=1e-10)
  sp<-coherent4cast:::.rc_frailty_spec("discrete",coherent4cast_control(),3)
  nd<-coherent4cast:::.rc_frailty_nodes(c(.3,-.1,-1,.4),sp)
  expect_equal(sum(exp(nd$h+nd$g)),1,tolerance=1e-13)
  expect_equal(sum(exp(nd$g)),1,tolerance=1e-13)
  expect_error(rc_frailty_loglik(th,1,-1,matrix(1,1,1),1,"lognormal"),class="coherent4cast_error")
})
test_that("main fit and prediction support alternative frailty and ordered covariates", {
  d<-rc_simulate(n=90,seed=120)
  # Use baseline x supplied by the simulator when available; otherwise a fixed synthetic covariate.
  d$x<-as.numeric(factor(d$id))%%2
  for(fam in c("gamma","lognormal","discrete")) {
    fit<-suppressWarnings(coherent4cast(~x,d,time_unit="years",frailty=fam,frailty_diagnostics=FALSE,control=coherent4cast_control(history_cap=0,maxit=100,quadrature=12)))
    expect_equal(fit$frailty,fam)
    one<-d[d$id==d$id[1],,drop=FALSE]
    p<-predict(fit,one,landmark=0,horizon=.4,type="all",quadrature=12)
    expect_equal(sum(p$estimate[p$target %in% c("risk","survival")]),1,tolerance=1e-10)
    expect_true(all(is.finite(p$estimate)))
  }
  fit<-suppressWarnings(coherent4cast(~x,d,time_unit="years",frailty_diagnostics=FALSE,control=coherent4cast_control(history_cap=0,phi_fixed=0)))
  fit$coefficients[grepl(":x$",names(fit$coefficients))]<-1
  one<-d[d$id==d$id[1],,drop=FALSE];one$x<-0
  frozen<-predict(fit,one,0,.5,type="risk")
  schedule<-data.frame(id=one$id[1],time=.25,x=2)
  evolving<-predict(fit,one,0,.5,type="risk",covariate_evolution="locf",newdata_future=schedule)
  term<-fit$coefficients[grepl("terminal:baseline",names(fit$coefficients))]
  expect_equal(evolving$estimate,unname(1-exp(-exp(term)*(.25+.25*exp(2)))),tolerance=1e-10)
  expect_gt(evolving$estimate,frozen$estimate)
  expect_error(predict(fit,one,0,.5,newdata_future=schedule),class="coherent4cast_error")
})
test_that("Gaussian future covariates condition on past only and preserve RNG", {
  skip_if_not_installed("nlme")
  set.seed(44);train<-expand.grid(time=0:4,id=1:45);b0<-rnorm(45);b1<-rnorm(45,0,.2)
  train$x<-.3*train$time+b0[train$id]+b1[train$id]*train$time+rnorm(nrow(train),0,.3)
  m<-rc_covariate_model(train,variables="x")
  past<-data.frame(start=c(0,1),x=c(.3,.5))
  before<-.Random.seed
  draws<-coherent4cast:::.rc_covariate_draws(m,past,c(1,2,3),"start",100,8)
  expect_identical(.Random.seed,before)
  expect_equal(draws,coherent4cast:::.rc_covariate_draws(m,past,c(1,2,3),"start",100,8))
  expect_equal(dim(draws[[1]]),c(3L,1L))
  expect_gt(stats::sd(vapply(draws,function(x)x[3,1],numeric(1))),0)
  expect_error(rc_covariate_model(rbind(train,train[1,]),variables="x"),class="coherent4cast_error")
})
test_that("conditional random intercept draws agree with Gaussian conjugacy", {
  m<-structure(list(variables="x",center=0,scale=1,random_slope=FALSE,
    models=list(x=list(beta=c(0,0),D=matrix(.5,1,1),sigma2=.25))),class="rc_covariate_model")
  a<-coherent4cast:::.rc_covariate_draws(m,data.frame(time=c(0,1),x=c(.5,1)),c(1,2),"time",10000,72)
  vals<-vapply(a,function(z)z[1,1],numeric(1))
  expect_lt(abs(mean(vals)-.6),.015)
  expect_lt(abs(stats::var(vals)-.1),.006)
  expect_true(all(vapply(a,function(z)z[1,1]==z[2,1],logical(1))))
})
test_that("lognormal forecast compression respects independent integral", {
  d<-rc_simulate(40,seed=82)
  f<-suppressWarnings(coherent4cast(~1,d,time_unit="years",frailty_diagnostics=FALSE,control=coherent4cast_control(phi_fixed=0,history_cap=0)))
  f$control$phi_fixed<-NULL;f$control$quadrature<-320L
  f$frailty<-"lognormal";f$frailty_spec<-coherent4cast:::.rc_frailty_spec("lognormal",f$control)
  f$coefficients<-c(log(c(.3,.2,.1)),log(1.5));f$phi<-1.5
  one<-d[d$id==d$id[1],,drop=FALSE]
  p<-predict(f,one,0,1,type="all")
  v<-log1p(1.5);oracle<-stats::integrate(function(x)exp(-.1*exp(x))*stats::dnorm(x,-v/2,sqrt(v)),-12,12,rel.tol=1e-12)$value
  expect_equal(p$estimate[p$target=="survival"],oracle,tolerance=1e-10)
  expect_equal(p$estimate[p$target=="burden_1"],3*(1-oracle),tolerance=1e-10)
  expect_true(all(p$frailty_compression_bound>=0))
  expect_lt(max(p$frailty_compression_bound),1e-10)
})
test_that("alternative frailty validates parameters and information contracts", {
  X<-matrix(1,2,1);theta<-c(-1,log(.3));u<-1:2;t<-c(1,2);event<-c(1,0)
  z<-rc_frailty_loglik(theta,u,t,X,event,"lognormal",information=FALSE)
  expect_null(z$information)
  expect_equal(length(rc_coarsen(z,owner=1:2,log_weights=c(0,0))$subject_loglik),2L)
  expect_error(rc_frailty_loglik(theta,c(1,3),t,X,event,"lognormal"),class="coherent4cast_error")
  expect_error(rc_frailty_loglik(theta,u,t,X,event,"discrete",phi_fixed=.3),class="coherent4cast_error")
  expect_error(rc_frailty_loglik(theta,u,t,X,event,"lognormal",quadrature=1),class="coherent4cast_error")
  expect_error(rc_frailty_loglik(theta[-1],u,t,X,event,"lognormal"),class="coherent4cast_error")
  expect_error(rc_frailty_loglik(theta,u,t,X,c(0,2),"lognormal"),class="coherent4cast_error")
  expect_error(rc_frailty_loglik(theta,u,t,X,event,"discrete",points=1),class="coherent4cast_error")
  expect_error(rc_covariate_model(data.frame(id=1:3,time=1,x=1),variables="x"),class="coherent4cast_error")
  fixed<-rc_frailty_loglik(-1,u,t,X,event,"lognormal",phi_fixed=0)
  expect_equal(fixed$subject_loglik,rc_loglik(-1,u,t,X,event,phi_fixed=0)$subject_loglik,tolerance=1e-13)
})
