# Reproducible revision validation; no target is imposed on fitted estimates.
library(recurcast)
out<-Sys.getenv("RECURCAST_BENCHMARK_OUT",unset="frailty_future_results")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
set.seed(1204)
trace<-list();fits<-list();a<-rc_readmission()
for(fam in c("gamma","lognormal","discrete")) {
  elapsed<-system.time(fit<-suppressWarnings(recurcast(a$formula,a$data,time_unit="years",event_types="1",frailty=fam,
    frailty_diagnostics=TRUE,control=recurcast_control(history_cap=1,quadrature=40,maxit=500))))[["elapsed"]]
  fits[[fam]]<-fit
  if(nrow(fit$audit$quadrature_refinement))write.csv(fit$audit$quadrature_refinement,file.path(out,paste0(fam,"_fit_quadrature_audit.csv")),row.names=FALSE)
  diag<-fit$frailty_diagnostics$null_variance_test
  trace[[fam]]<-data.frame(family=fam,n=fit$n,logLik=fit$loglik,AIC=-2*fit$loglik+2*length(coef(fit)),variance=fit$phi,
    parameters=length(coef(fit)),converged=fit$convergence$valid,gradient=fit$convergence$gradient,
    minimum_information=min(fit$convergence$information_eigenvalues),boundary_score=diag$statistic,boundary_p=diag$p.value,elapsed_seconds=elapsed)
}
write.csv(do.call(rbind,trace),file.path(out,"readmission_frailty_families.csv"),row.names=FALSE)
saveRDS(fits,file.path(out,"readmission_frailty_fits.rds"))
# Direct integration and independent finite differences, retained as numerical evidence.
th<-c(log(.4),log(.7));v<-log1p(.7)
truth<-log(.4)+log(integrate(function(t)exp(t)*exp(-.8*exp(t))*dnorm(t,-v/2,sqrt(v)),-12,12,rel.tol=1e-12)$value)
quad<-do.call(rbind,lapply(c(10,20,40,80),function(q){z<-rc_frailty_loglik(th,1,2,matrix(1,1,1),1,"lognormal",quadrature=q);data.frame(nodes=q,loglik=z$subject_loglik,reference=truth,absolute_error=abs(z$subject_loglik-truth))}))
write.csv(quad,file.path(out,"lognormal_quadrature.csv"),row.names=FALSE)
stopifnot(tail(quad$absolute_error,1)<1e-10)
# A separate synthetic model gives an independent closed-form forecast check.
d<-rc_simulate(n=150,seed=321,history_cap=0,dispersion=0)
f<-suppressWarnings(recurcast(~x,d,time_unit="years",frailty_diagnostics=FALSE,control=recurcast_control(phi_fixed=0,history_cap=0)))
f$coefficients[grepl(":x$",names(coef(f)))]<-.7
past<-data.frame(id=999,start=0,stop=1,event_type=0,terminal=0,x=0)
schedule<-data.frame(id=999,time=c(1,1.5),x=c(0,1))
p<-predict(f,past,landmark=1,horizon=1,type="all",covariate_evolution="locf",newdata_future=schedule)
lambda<-exp(coef(f)[grepl("terminal:baseline",names(coef(f)))])
oracle<-1-exp(-lambda*(.5+.5*exp(.7)))
stopifnot(abs(p$estimate[p$target=="risk"]-oracle)<1e-10)
write.csv(data.frame(estimate=p$estimate[p$target=="risk"],oracle=unname(oracle),absolute_error=abs(p$estimate[p$target=="risk"]-oracle)),file.path(out,"ordered_covariate_oracle.csv"),row.names=FALSE)
# Fit a longitudinal covariate model and integrate future random coefficients.
train<-expand.grid(time=0:4,id=1:60);b0<-rnorm(60,0,.5);b1<-rnorm(60,0,.15)
train$x<-.2*train$time+b0[train$id]+b1[train$id]*train$time+rnorm(nrow(train),0,.2)
model<-rc_covariate_model(train,variables="x")
measured<-data.frame(id=999,time=c(0,1),x=c(0,.2))
future<-do.call(rbind,lapply(c(8L,16L,32L),function(grid){
  pp<-predict(f,past,1,1,type="all",covariate_evolution="joint_submodel",covariate_model=model,
    covariate_history=measured,covariate_draws=500,covariate_grid=grid,seed=91)
  pp$grid<-grid;pp
}))
write.csv(future,file.path(out,"future_covariate_grid.csv"),row.names=FALSE)
stopifnot(all(is.finite(future$estimate)),all(future$covariate_mcse>=0))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("Frailty comparisons, integration oracles and future-covariate grids completed.\n")
