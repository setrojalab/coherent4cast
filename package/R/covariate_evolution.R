#' Estimate a Gaussian longitudinal model for future covariates
#' @export
rc_covariate_model <- function(data,unit="id",time="time",variables,random_slope=TRUE,method=c("REML","ML")) {
  method<-match.arg(method)
  if(!requireNamespace("nlme",quietly=TRUE)) .rc_error("optional_dependency","Install nlme to estimate the longitudinal covariate submodel.")
  needed<-c(unit,time,variables)
  if(!is.data.frame(data)||!length(variables)||anyDuplicated(variables)||!all(needed %in% names(data))||anyNA(data[,needed,drop=FALSE])) .rc_error("covariate_model","Supply complete unit, time and named numeric covariate columns.")
  if(!is.numeric(data[[time]])||any(!is.finite(data[[time]]))||any(!vapply(data[variables],function(x)is.numeric(x)&&all(is.finite(x)),logical(1)))) .rc_error("covariate_model","Longitudinal times and responses must be finite numeric values.")
  if(anyDuplicated(data[,c(unit,time)])) .rc_error("covariate_model","One measurement per unit and time is required.")
  center<-mean(data[[time]]);scale<-stats::sd(data[[time]])
  if(!is.finite(scale)||scale<=0||length(unique(data[[unit]]))<3L) .rc_error("covariate_model","At least three units and varying measurement times are required.")
  models<-lapply(variables,function(v){
    frame<-data.frame(.id=factor(data[[unit]]),.time=(data[[time]]-center)/scale,.y=data[[v]])
    f<-nlme::lme(.y~.time,data=frame,random=if(random_slope)~1+.time|.id else ~1|.id,method=method,
      control=nlme::lmeControl(maxIter=100,msMaxIter=200,returnObject=FALSE))
    D<-as.matrix(nlme::getVarCov(f,type="random.effects"));list(beta=as.numeric(nlme::fixef(f)),D=D,sigma2=f$sigma^2,fit=f)
  });names(models)<-variables
  structure(list(models=models,variables=variables,unit=unit,time=time,center=center,scale=scale,random_slope=random_slope,
    method=method,training_units=length(unique(data[[unit]])),assumptions=c("Gaussian linear random-coefficient longitudinal submodels",
    "variables conditionally independent across submodels", "future covariate evolution has no feedback from future events conditional on supplied past history",
    "event intensities condition on observed predictable covariates; this is not a jointly estimated latent measurement-error survival model")),class="rc_covariate_model")
}
.rc_with_seed <- function(seed,code) {
  if(length(seed)!=1L||!is.finite(seed)||seed<0||seed>.Machine$integer.max||seed!=floor(seed)) .rc_error("seed","Use one nonnegative integer seed.")
  had<-exists(".Random.seed",.GlobalEnv,inherits=FALSE);old<-if(had)get(".Random.seed",.GlobalEnv)else NULL
  on.exit(if(had)assign(".Random.seed",old,.GlobalEnv)else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed);force(code)
}
.rc_covariate_draws <- function(model,past,times,time_column,draws,seed) {
  if(!inherits(model,"rc_covariate_model")) .rc_error("covariate_model","Supply an rc_covariate_model object for joint_submodel prediction.")
  if(!all(model$variables %in% names(past))) .rc_error("covariate_model","Past history lacks a modeled covariate.")
  .rc_with_seed(seed,{
    out<-vector("list",draws);for(b in seq_len(draws))out[[b]]<-matrix(0,length(times),length(model$variables),dimnames=list(NULL,model$variables))
    for(v in model$variables){
      m<-model$models[[v]];d<-past[!duplicated(past[[time_column]]),,drop=FALSE]
      t<-(d[[time_column]]-model$center)/model$scale;X<-cbind(1,t);Z<-if(model$random_slope)X else matrix(1,nrow(d),1)
      # Gaussian conditioning, expressed in measurement space; invertible because sigma2 > 0.
      V<-Z%*%m$D%*%t(Z)+diag(m$sigma2,nrow(d));C<-m$D%*%t(Z)%*%solve(V)
      mu<-as.vector(C%*%(d[[v]]-X%*%m$beta));Sigma<-m$D-C%*%Z%*%m$D
      e<-eigen((Sigma+t(Sigma))/2,symmetric=TRUE)
      if(min(e$values) < -1e-8) .rc_error("covariate_covariance","Conditional random-effect covariance is not positive semidefinite.")
      L<-sweep(e$vectors,2,sqrt(pmax(0,e$values)),"*");Xf<-cbind(1,(times-model$center)/model$scale);Zf<-if(model$random_slope)Xf else matrix(1,length(times),1)
      for(b in seq_len(draws))out[[b]][,v]<-as.vector(Xf%*%m$beta+Zf%*%(mu+L%*%stats::rnorm(length(mu))))
    };out
  })
}
.rc_future_scenarios <- function(object,d,rows,current,ident,landmark,horizon,mode,future,model,draws,grid,seed,measurement_history=NULL) {
  sc<-object$schema;last<-d[current,,drop=FALSE];end<-landmark+horizon
  times<-sort(unique(c(landmark,object$breaks[object$breaks>landmark & object$breaks<end])))
  frames<-list()
  if(mode=="freeze") {
    if(!is.null(future)) .rc_error("future_covariates","newdata_future requires covariate_evolution='locf'.")
    frames[[1]]<-last[rep(1,length(times)),,drop=FALSE]
  } else if(mode=="locf") {
    if(!is.null(future)) {
      vars<-all.vars(object$formula)
      if(!is.data.frame(future)||!all(c(sc$unit,"time",vars) %in% names(future))||anyNA(future[,c(sc$unit,"time",vars),drop=FALSE])||!is.numeric(future$time)||any(!is.finite(future$time))||any(future$time<landmark)) .rc_error("future_covariates","Future schedules need unit, time >= landmark, and all formula covariates, without missing values.")
      ff<-future[as.character(future[[sc$unit]])==ident & future$time<end,,drop=FALSE]
      if(anyDuplicated(ff$time)) .rc_error("future_covariates","A unit cannot have duplicate future update times.")
      ff<-ff[order(ff$time),,drop=FALSE];times<-sort(unique(c(times,ff$time)))
      f<-last[rep(1,length(times)),,drop=FALSE]
      for(j in seq_along(times)){ix<-which(ff$time<=times[j]);if(length(ix))f[j,vars]<-ff[max(ix),vars,drop=FALSE]}
      frames[[1]]<-f
    }else frames[[1]]<-last[rep(1,length(times)),,drop=FALSE]
  }else {
    if(!is.null(future)) .rc_error("future_covariates","A deterministic schedule and a stochastic joint_submodel are separate forecast scenarios.")
    if(!inherits(model,"rc_covariate_model")) .rc_error("covariate_model","Supply a fitted rc_covariate_model.")
    if(!all(model$variables %in% all.vars(object$formula))) .rc_error("covariate_model","Modeled variables must occur in the event-model formula.")
    if(length(draws)!=1L||!is.finite(draws)||draws<2||draws!=floor(draws)||length(grid)!=1L||!is.finite(grid)||grid<1||grid!=floor(grid)) .rc_error("covariate_model","Use at least two Monte Carlo draws and one time-grid segment.")
    times<-sort(unique(c(times,head(seq(landmark,end,length.out=grid+1L),-1L))))
    needed<-c(sc$unit,"time",model$variables)
    if(!is.data.frame(measurement_history)||!all(needed %in% names(measurement_history))||anyNA(measurement_history[,needed,drop=FALSE])||!is.numeric(measurement_history$time)||any(!is.finite(measurement_history$time))) .rc_error("covariate_history","Supply covariate_history with unit, time and measured variables; event intervals are not independent biomarker measurements.")
    past<-measurement_history[as.character(measurement_history[[sc$unit]])==ident & measurement_history$time<=landmark,,drop=FALSE]
    if(!nrow(past)||anyDuplicated(past$time)) .rc_error("covariate_history","Each forecast unit needs unique measurement times at or before landmark.")
    sampled<-.rc_covariate_draws(model,past,times,"time",draws,seed)
    for(b in seq_len(draws)){frames[[b]]<-last[rep(1,length(times)),,drop=FALSE];frames[[b]][,model$variables]<-sampled[[b]]}
  }
  mm<-lapply(frames,function(f){mf<-stats::model.frame(object$design$terms,f,na.action=stats::na.fail,xlev=object$design$xlevels);stats::model.matrix(object$design$terms,mf,contrasts.arg=object$design$contrasts)[,object$design$covariate_columns,drop=FALSE]})
  list(cuts=c(times,end),covariates=mm,mode=mode)
}
