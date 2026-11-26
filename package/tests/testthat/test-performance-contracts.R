test_that("the kernel has a fixed-dimensional replicated-workload contract", {
  x <- matrix(c(1,1,1,1,-1,1,-1,1),ncol=2)
  id <- c(1L,1L,2L,2L);dt<-c(.5,1,.8,.4);event<-c(1,0,0,1)
  theta<-c(-.4,.2,log(.7))
  one<-cs_loglik(theta,id,dt,x,event)
  copied<-cs_loglik(theta,c(id,id+2L),rep(dt,2),rbind(x,x),rep(event,2))
  expect_equal(copied$objective,one$objective,tolerance=1e-13)
  expect_equal(copied$gradient,one$gradient,tolerance=1e-13)
  expect_equal(copied$information,2*one$information,tolerance=1e-12)
  expect_equal(copied$scores,rbind(one$scores,one$scores),tolerance=1e-13)
  # Permuting records preserves the likelihood because rows carry their unit ID.
  permutation<-c(4,1,3,2)
  shuffled<-cs_loglik(theta,id[permutation],dt[permutation],x[permutation,,drop=FALSE],event[permutation])
  expect_equal(shuffled$loglik,one$loglik,tolerance=1e-13)
  expect_equal(shuffled$scores,one$scores,tolerance=1e-13)
})

test_that("gamma integration matches a direct independent R calculation", {
  id<-c(1,1,2,2,3);dt<-c(.2,.4,.7,.3,.9)
  x<-cbind(1,c(0,1,0,0,1));event<-c(1,1,0,1,0)
  theta<-c(-.2,.3,log(.4));a<-exp(-theta[3])
  eta<-drop(x%*%theta[1:2]);v<-dt*exp(eta)
  A<-drop(rowsum(v,id,reorder=FALSE));d<-drop(rowsum(event,id,reorder=FALSE))
  ll<-drop(rowsum(event*eta,id,reorder=FALSE))+lgamma(a+d)-lgamma(a)+a*log(a)-(a+d)*log(a+A)
  A1<-rowsum(x*v,id,reorder=FALSE);ev<-rowsum(x*event,id,reorder=FALSE)
  sb<-ev-A1*((a+d)/(a+A))
  sp<- -a*(digamma(a+d)-digamma(a)+log(a)+1-log(a+A)-(a+d)/(a+A))
  out<-cs_loglik(theta,id,dt,x,event)
  expect_equal(out$subject_loglik,unname(ll),tolerance=1e-12)
  expect_equal(unname(out$scores),unname(cbind(sb,sp)),tolerance=1e-12)
  individual<-lapply(unique(id),function(i){
    take<-id==i
    cs_loglik(theta,rep(1,sum(take)),dt[take],x[take,,drop=FALSE],event[take])
  })
  expect_equal(sum(vapply(individual,`[[`,numeric(1),"loglik")),out$loglik,tolerance=1e-12)
  expect_equal(Reduce(`+`,lapply(individual,`[[`,"information")),out$information,tolerance=1e-12)
})

test_that("gamma quadrature converges for a known posterior transform", {
  # This tests quadrature independently of fitted coefficients or optimizer noise.
  shape<-2.5;rate<-3.2;hazard<-.7
  truth<-(rate/(rate+hazard))^shape
  estimate<-vapply(c(5L,15L,30L),function(k){
    q<-coherentSurv:::.cs_gamma_nodes(shape,rate,k)
    expect_equal(sum(q$weights),1,tolerance=1e-13)
    expect_true(all(q$nodes>0))
    sum(q$weights*exp(-hazard*q$nodes))
  },numeric(1))
  expect_lt(abs(estimate[3]-truth),1e-12)
  expect_lt(abs(estimate[2]-truth),1e-10)
  expect_lt(abs(estimate[3]-truth),abs(estimate[1]-truth))
})

test_that("forecast coordinates share an analytic constant-rate law", {
  d<-cs_simulate(30,seed=24,time_unit="years")
  fitted<-suppressWarnings(cs_fit(~1,d,event_types=c("1","2"),time_unit="years",
    control=cs_control(phi_fixed=0,history_cap=1L,maxit=100L)))
  # Set a declared constant-rate law to test propagation independently of fitting.
  rates<-c(.3,.2,.1);fitted$coefficients[]<-0
  for(k in 1:3)fitted$coefficients[(k-1L)*fitted$design$block_size+1L]<-log(rates[k])
  subject<-d[d$id==1,,drop=FALSE]
  output<-predict(fitted,subject,landmark=0,horizon=c(.5,2),type="all")
  for(h in c(.5,2)) {
    z<-output[output$horizon==h,];v<-setNames(z$estimate,z$target)
    alive<-exp(-.1*h);restricted<-(1-alive)/.1
    expect_equal(unname(v["risk"]+v["survival"]),1,tolerance=1e-12)
    expect_equal(unname(v[c("risk","survival","rmst","burden_1","burden_2")]),
      c(1-alive,alive,restricted,.3*restricted,.2*restricted),tolerance=1e-11)
    expect_true(all(z$series_error_bound>=0))
  }
  fitted$phi<-.4
  fitted$control$phi_fixed<-.4
  mixed<-predict(fitted,subject,landmark=0,horizon=2,type="all",quadrature=30)
  values<-setNames(mixed$estimate,mixed$target)
  alive<-(1+.4*.1*2)^(-1/.4)
  expect_equal(unname(values["survival"]),alive,tolerance=1e-11)
  expect_equal(unname(values["burden_1"]),.3/.1*(1-alive),tolerance=1e-11)
})
