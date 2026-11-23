#' Validate predictions at a fixed horizon
#' @export
validate <- function(object,...) UseMethod("validate")

#' @export
validate.default <- function(object,time,status,horizon,censor_survival=NULL,threshold=0.2,strata=NULL,...) {
  predicted<-as.numeric(object);n<-length(predicted)
  if(!n||length(time)!=n||length(status)!=n||any(!is.finite(c(predicted,time)))||any(time<0)||any(!status %in% c(0,1))||anyNA(status)||any(predicted<0|predicted>1)) .rc_error("validation","Provide aligned finite risks, nonnegative times and binary status.")
  if(length(horizon)!=1L||!is.finite(horizon)||horizon<=0||length(threshold)!=1L||!is.finite(threshold)||threshold<=0||threshold>=1) .rc_error("validation","Require a positive horizon and a threshold strictly between zero and one.")
  if(!is.null(strata)) {
    if(length(strata)!=n||anyNA(strata)) .rc_error("strata","Supply one nonmissing stratum per observation.")
    return(lapply(split(seq_len(n),strata),function(i)validate(predicted[i],time[i],status[i],horizon,
      if(is.null(censor_survival))NULL else list(case=censor_survival$case[i],horizon=censor_survival$horizon[i]),threshold)))
  }
  case<-status==1 & time<=horizon;control<-time>horizon
  if(is.null(censor_survival)) {
    times<-sort(unique(time));at<-vapply(times,function(t)sum(time>=t),numeric(1));dc<-vapply(times,function(t)sum(time==t & status==0),numeric(1));surv<-cumprod(1-dc/at)
    gcase<-vapply(time,function(t){i<-which(times<t);if(length(i))surv[max(i)]else 1},numeric(1))
    i<-which(times<=horizon);gh<-if(length(i))surv[max(i)]else 1
    ghorizon<-rep(gh,n);assumption<-"Marginal independent censoring for the Kaplan-Meier censoring estimate. Use supplied conditional probabilities when needed."
  } else {
    gcase<-censor_survival$case;ghorizon<-censor_survival$horizon
    if(length(gcase)!=n||length(ghorizon)!=n) .rc_error("censoring","Supply case and horizon censoring-survival vectors of length n.")
    assumption<-"Supplied censoring probabilities must be valid for the target population and conditioning history."
  }
  if(any(!is.finite(c(gcase,ghorizon)))||any(gcase<0|gcase>1)||any(ghorizon<0|ghorizon>1)||any(gcase[case]<=0)||any(ghorizon[control]<=0)) .rc_error("positivity","Needed censoring-survival probabilities must be strictly positive and at most one.")
  w<-numeric(n);w[case]<-1/gcase[case];w[control]<-1/ghorizon[control];y<-as.numeric(case);valid<-w>0
  if(!any(valid)) .rc_error("validation_support","No observed cases or horizon survivors are available; use a supported evaluation horizon or dataset.")
  brier<-sum(w*(y-predicted)^2)/n
  wc<-sum(w[case]);wn<-sum(w[control]);auc<-NA_real_
  if(wc>0&&wn>0) {
    ranked<-order(predicted);ordered<-split(ranked,cumsum(c(TRUE,diff(predicted[ranked])!=0)));cum_control<-0;num<-0
    for(i in ordered){cw<-sum(w[i][case[i]]);nw<-sum(w[i][control[i]]);num<-num+cw*(cum_control+.5*nw);cum_control<-cum_control+nw}
    auc<-num/(wc*wn)
  }
  positive<-predicted>=threshold;tp<-sum(w[case&positive]);fn<-sum(w[case&!positive]);fp<-sum(w[control&positive]);tn<-sum(w[control&!positive])
  ratio<-function(a,b)if(b>0)a/b else NA_real_
  citl<-slope<-ici<-NA_real_;cal_status<-"insufficient_cases_or_controls"
  if(wc>0&&wn>0&&sum(valid)>3) {
    z<-stats::qlogis(pmin(1-1e-8,pmax(1e-8,predicted[valid])));yy<-y[valid];ww<-w[valid]
    fits<-tryCatch(suppressWarnings(list(stats::glm(yy~1+offset(z),weights=ww,family=stats::quasibinomial()),stats::glm(yy~z,weights=ww,family=stats::quasibinomial()))),error=function(e)NULL)
    if(!is.null(fits)){citl<-unname(stats::coef(fits[[1]])[1]);slope<-unname(stats::coef(fits[[2]])[2]);cal_status<-if(all(vapply(fits,function(f)isTRUE(f$converged),logical(1))))"estimated" else "nonconverged"}
    if(length(unique(predicted[valid]))>=4) {
      local<-tryCatch(suppressWarnings(stats::loess(yy~pp,data=data.frame(yy=yy,pp=predicted[valid]),weights=ww,span=.75,degree=1,control=stats::loess.control(surface="direct"))),error=function(e)NULL)
      if(!is.null(local)){cal<-as.numeric(stats::predict(local));if(all(is.finite(cal)))ici<-sum(ww*abs(pmin(1,pmax(0,cal))-predicted[valid]))/sum(ww)}
    }
  }
  values<-c(Brier=brier,AUROC=auc,CITL=citl,calibration_slope=slope,ICI=ici,sensitivity=ratio(tp,tp+fn),specificity=ratio(tn,tn+fp),PPV=ratio(tp,tp+fp),NPV=ratio(tn,tn+fn),net_benefit=(tp-fp*threshold/(1-threshold))/n)
  list(metrics=data.frame(metric=names(values),estimate=unname(values),defined=is.finite(values)),
       confusion=matrix(c(tp,fn,fp,tn),2,2,byrow=TRUE,dimnames=list(outcome=c("event","no_event"),prediction=c("positive","negative"))),
       horizon=horizon,threshold=threshold,n=n,usable=sum(valid),cases=sum(case),controls=sum(control),
       calibration_status=cal_status,assumptions=assumption,weights=w)
}

#' @export
validate.recurcast <- function(object,time,status,horizon,predicted,censor_survival=NULL,...) {
  validate.default(predicted,time,status,horizon,censor_survival,...)
}

#' @export
print.recurcast <- function(x,...) {
  if(identical(x$frailty_requested,"none"))cat("RecurCast no-frailty proportional-hazards model\n") else cat("RecurCast shared",if(is.null(x$frailty))"gamma"else x$frailty,"model\n")
  cat(x$n,"independent units;",x$design$K,"nonterminal types; time:",x$schema$time_unit,"\n")
  cat("Frailty variance:",format(x$phi)," Convergence:",x$convergence$code,"\n")
  cat("Observation:",x$schema$observation,"; data fingerprint:",x$audit$data_hash,"\n")
  invisible(x)
}

#' @export
summary.recurcast <- function(object,...) {
  se<-sqrt(pmax(0,diag(object$covariance)));estimate<-object$coefficients
  out<-list(coefficients=data.frame(term=names(estimate),estimate=unname(estimate),std.error=se,
    conf.low=estimate-stats::qnorm(.975)*se,conf.high=estimate+stats::qnorm(.975)*se),
    phi=object$phi,loglik=object$loglik,n=object$n,convergence=object$convergence,
    interval_scope="Regular asymptotic sandwich intervals; not a finite-sample guarantee.",frailty_diagnostics=object$frailty_diagnostics,audit=object$audit)
  class(out)<-"summary.recurcast";out
}
#' @export
print.summary.recurcast <- function(x,...) {print(x$coefficients,row.names=FALSE);cat("\n",x$interval_scope,"\nConvergence:",x$convergence$code,"\n");if(!is.null(x$frailty_diagnostics)){cat("\nNo-frailty boundary score test:\n");print(x$frailty_diagnostics$null_variance_test);print(x$frailty_diagnostics$comparison,row.names=FALSE)};invisible(x)}
#' @export
coef.recurcast <- function(object,...)object$coefficients
#' @export
vcov.recurcast <- function(object,...)object$covariance
#' @export
logLik.recurcast <- function(object,...)structure(object$loglik,df=length(object$coefficients),nobs=object$n,class="logLik")
#' @export
tidy <- function(x,...)UseMethod("tidy")
#' @export
tidy.recurcast <- function(x,...)summary(x)$coefficients
#' @export
glance <- function(x,...)UseMethod("glance")
#' @export
glance.recurcast <- function(x,...)data.frame(n=x$n,logLik=x$loglik,AIC=-2*x$loglik+2*length(x$coefficients),frailty_variance=x$phi,converged=x$convergence$valid,gradient_max=x$convergence$gradient)
#' @export
augment <- function(x,...)UseMethod("augment")
#' @export
augment.recurcast <- function(x,newdata=NULL,landmark,horizon,...)predict(x,newdata,landmark,horizon,type="all",...)
#' @export
bounds <- function(object,...)UseMethod("bounds")
#' @export
bounds.default <- function(object,...)rc_bounds(p=object,...)
#' @export
bounds.recurcast <- function(object,...) .rc_error("bounds_schema","Recording-sensitivity bounds require an explicit finite-support recorded law. Use rc_bounds().")
#' @export
conformal <- function(object,...)rc_conformal(scores=object,...)
#' @export
plot.recurcast_prediction <- function(x,unit=unique(x$unit)[1],target="risk",...) {
  d<-x[x$unit==unit & x$target==target,,drop=FALSE]
  if(!nrow(d)) .rc_error("plot","No matching unit and target.")
  d<-d[order(d$horizon),];graphics::plot(d$horizon,d$estimate,type="b",xlab=paste("Horizon",d$time_unit[1]),ylab=target,...);invisible(x)
}
#' @export
autoplot <- function(object,...)UseMethod("autoplot")
#' @export
autoplot.recurcast_prediction <- function(object,unit=unique(object$unit)[1],target="risk",...) {
  if(!requireNamespace("ggplot2",quietly=TRUE)) .rc_error("optional_dependency","Install ggplot2 for autoplot(); base plot() is also available.")
  d<-object[object$unit==unit & object$target==target,,drop=FALSE]
  ggplot2::ggplot(d,ggplot2::aes(x=horizon,y=estimate))+ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::labs(x=paste("Horizon",d$time_unit[1]),y=target)
}
utils::globalVariables(c("horizon","estimate"))
