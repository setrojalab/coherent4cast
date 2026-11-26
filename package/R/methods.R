#' Validate predictions at a fixed horizon
#' @export
validate <- function(object,...) UseMethod("validate")

#' @export
validate.default <- function(object,time,status,horizon,censor_survival=NULL,threshold=0.2,strata=NULL,...) {
  if(length(list(...))) .cs_error("unused_argument","Unsupported validation argument. Supply follow-up durations measured from the assessment landmark; landmark is not a validation argument.")
  predicted<-as.numeric(object);n<-length(predicted)
  if(!n||length(time)!=n||length(status)!=n||any(!is.finite(c(predicted,time)))||any(time<0)||any(!status %in% c(0,1))||anyNA(status)||any(predicted<0|predicted>1)) .cs_error("validation","Provide aligned finite risks, nonnegative times and binary status.")
  if(length(horizon)!=1L||!is.finite(horizon)||horizon<=0||length(threshold)!=1L||!is.finite(threshold)||threshold<=0||threshold>=1) .cs_error("validation","Require a positive horizon and a threshold strictly between zero and one.")
  if(!is.null(strata)) {
    if(length(strata)!=n||anyNA(strata)) .cs_error("strata","Supply one nonmissing stratum per observation.")
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
    if(length(gcase)!=n||length(ghorizon)!=n) .cs_error("censoring","Supply case and horizon censoring-survival vectors of length n.")
    assumption<-"Supplied censoring probabilities must be valid for the target population and conditioning history."
  }
  if(any(!is.finite(c(gcase,ghorizon)))||any(gcase<0|gcase>1)||any(ghorizon<0|ghorizon>1)||any(gcase[case]<=0)||any(ghorizon[control]<=0)) .cs_error("positivity","Needed censoring-survival probabilities must be strictly positive and at most one.")
  w<-numeric(n);w[case]<-1/gcase[case];w[control]<-1/ghorizon[control];y<-as.numeric(case);valid<-w>0
  if(!any(valid)) .cs_error("validation_support","No observed cases or horizon survivors are available; use a supported evaluation horizon or dataset.")
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
    intercept_fit<-tryCatch(suppressWarnings(stats::glm(yy~1+offset(z),weights=ww,family=stats::quasibinomial())),error=function(e)NULL)
    intercept_valid<-!is.null(intercept_fit)&&isTRUE(intercept_fit$converged)&&is.finite(stats::coef(intercept_fit)[1])
    if(intercept_valid)citl<-unname(stats::coef(intercept_fit)[1])
    # For an intercept and one predictor, disjoint outcome ranges (including a
    # shared boundary) characterize complete or quasi-complete separation.
    # Positive IPCW weights do not change that criterion. glm() can report
    # convergence at a large finite iterate even though the slope has no MLE.
    if(length(unique(z))<2L)cal_status<-"constant_predictions" else if(max(z[yy==0])<=min(z[yy==1])||max(z[yy==1])<=min(z[yy==0]))cal_status<-"separated" else {
      slope_fit<-tryCatch(suppressWarnings(stats::glm(yy~z,weights=ww,family=stats::quasibinomial())),error=function(e)NULL)
      slope_valid<-!is.null(slope_fit)&&isTRUE(slope_fit$converged)&&is.finite(stats::coef(slope_fit)[2])
      if(slope_valid)slope<-unname(stats::coef(slope_fit)[2])
      cal_status<-if(intercept_valid&&slope_valid)"estimated" else "nonconverged"
    }
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
validate.coherentSurv <- function(object,time,status,horizon,predicted,censor_survival=NULL,...) {
  validate.default(predicted,time,status,horizon,censor_survival,...)
}

#' @export
print.coherentSurv <- function(x,...) {
  if(identical(x$frailty_requested,"none"))cat("coherentSurv no-frailty proportional-hazards model\n") else cat("coherentSurv shared",if(is.null(x$frailty))"gamma"else x$frailty,"model\n")
  cat(x$n,"independent units;",x$design$K,"nonterminal types; time:",x$schema$time_unit,"\n")
  cat("Frailty variance:",format(x$phi)," Convergence:",x$convergence$code,"\n")
  cat("Observation:",x$schema$observation,"; data fingerprint:",x$audit$data_hash,"\n")
  invisible(x)
}

#' @export
summary.coherentSurv <- function(object,...) {
  se<-sqrt(pmax(0,diag(object$covariance)));estimate<-object$coefficients
  regular<-isTRUE(object$convergence$valid)
  inference_status<-if(regular)"regular_asymptotic" else "unavailable_nonregular_fit"
  lower<-upper<-rep(NA_real_,length(estimate))
  if(regular){lower<-estimate-stats::qnorm(.975)*se;upper<-estimate+stats::qnorm(.975)*se}
  out<-list(coefficients=data.frame(term=names(estimate),estimate=unname(estimate),std.error=se,
    conf.low=lower,conf.high=upper,inference_status=inference_status),
    phi=object$phi,loglik=object$loglik,n=object$n,convergence=object$convergence,
    inference_status=inference_status,
    interval_scope=if(regular)"Regular asymptotic sandwich intervals; not a finite-sample guarantee." else "Regular sandwich intervals withheld: the fit did not meet convergence and regularity criteria. Standard errors are numerical diagnostics only.",frailty_diagnostics=object$frailty_diagnostics,audit=object$audit)
  class(out)<-"summary.coherentSurv";out
}
#' @export
print.summary.coherentSurv <- function(x,...) {print(x$coefficients,row.names=FALSE);cat("\n",x$interval_scope,"\nConvergence:",x$convergence$code,"\n");if(!is.null(x$audit$estimand)){cat("\nAnalyst-declared estimand (not validated by the software):\n");print(x$audit$estimand)};if(!is.null(x$frailty_diagnostics)){cat("\nNo-frailty boundary score test:\n");print(x$frailty_diagnostics$null_variance_test);print(x$frailty_diagnostics$comparison,row.names=FALSE)};invisible(x)}
#' @export
coef.coherentSurv <- function(object,...)object$coefficients
#' @export
vcov.coherentSurv <- function(object,...)object$covariance
#' @export
logLik.coherentSurv <- function(object,...)structure(object$loglik,df=length(object$coefficients),nobs=object$n,class="logLik")

.cs_confint_selection <- function(terms,parm,level,dots) {
  if(length(dots)) .cs_error("unused_argument","An unsupported confidence-interval argument was supplied.")
  if(!is.numeric(level)||length(level)!=1L||!is.finite(level)||level<=0||level>=1)
    .cs_error("confidence_level","level must be a single number strictly between zero and one.")
  if(!length(terms)||anyNA(terms)||anyDuplicated(terms))
    .cs_error("confidence_parameters","The fitted object must have distinct named parameters or targets.")
  if(is.null(parm))return(seq_along(terms))
  if(is.character(parm)) {
    if(anyNA(parm)||any(!parm%in%terms)) .cs_error("confidence_parameters","parm contains an unknown parameter or target name.")
    return(match(parm,terms))
  }
  if(!is.numeric(parm)||any(!is.finite(parm))||any(parm!=floor(parm))||
     any(abs(parm)>length(terms))||(any(parm>0)&&any(parm<0)))
    .cs_error("confidence_parameters","parm must contain valid names or integer parameter indices; positive and negative indices cannot be mixed.")
  seq_along(terms)[parm]
}

.cs_confint_matrix <- function(lower,upper,terms,selected,level,status,interval,scope) {
  probabilities<-c((1-level)/2,(1+level)/2)
  result<-cbind(lower,upper)[selected,,drop=FALSE]
  dimnames(result)<-list(terms[selected],paste0(format(100*probabilities,trim=TRUE,scientific=FALSE,digits=3)," %"))
  if(length(status)==1L)status<-rep(status,length(terms))
  attr(result,"inference_status")<-stats::setNames(status[selected],terms[selected])
  attr(result,"level")<-level;attr(result,"interval")<-interval
  attr(result,"inference_scope")<-scope
  result
}

.cs_confint_wald <- function(estimate,covariance,regular,parm,level,dots) {
  terms<-names(estimate);selected<-.cs_confint_selection(terms,parm,level,dots)
  if(!is.matrix(covariance)||!identical(dim(covariance),c(length(estimate),length(estimate))))
    .cs_error("confidence_covariance","The fitted covariance does not match the parameter vector.")
  variance<-diag(covariance);valid<-regular&is.finite(estimate)&is.finite(variance)&variance>0
  lower<-upper<-rep(NA_real_,length(estimate))
  lower[valid]<-estimate[valid]-stats::qnorm((1+level)/2)*sqrt(variance[valid])
  upper[valid]<-estimate[valid]+stats::qnorm((1+level)/2)*sqrt(variance[valid])
  status<-if(!regular)rep("unavailable_nonregular_fit",length(estimate)) else
    ifelse(valid,"regular_asymptotic","unavailable_variance")
  .cs_confint_matrix(lower,upper,terms,selected,level,status,"wald",
    "Pointwise regular asymptotic sandwich intervals; nonregular fits and invalid variances are withheld. No finite-sample guarantee.")
}

#' @export
confint.coherentSurv <- function(object,parm=NULL,level=.95,...) {
  .cs_confint_wald(object$coefficients,object$covariance,isTRUE(object$convergence$valid),parm,level,list(...))
}

#' @export
confint.coherentSurv_competing <- function(object,parm=NULL,level=.95,...) {
  regular<-isTRUE(object$convergence$valid)&&length(object$fits)>0L&&
    all(vapply(object$fits,function(f)isTRUE(f$convergence$valid),logical(1)))
  .cs_confint_wald(object$coefficients,object$covariance,regular,parm,level,list(...))
}

#' @export
confint.coherentSurv_targeted <- function(object,parm=NULL,level=.95,interval="auto",draws=1000L,seed=1L,...) {
  selected<-.cs_confint_selection(object$target_names,parm,level,list(...))
  result<-.cs_targeted_predict(object,interval=interval,level=level,draws=draws,seed=seed)
  out<-.cs_confint_matrix(result$lower,result$upper,result$target,selected,level,
    result$inference_status,unique(result$interval),attr(result,"inference_scope"))
  attr(out,"requested_interval")<-interval
  out
}

#' @export
confint.coherentSurv_causal <- function(object,parm=NULL,level=.95,
    contrast=c("difference","ratio","log_ratio"),interval=c("wald","multiplier","none"),draws=2000L,seed=1L,...) {
  selected<-.cs_confint_selection(object$target_names,parm,level,list(...))
  contrast<-match.arg(contrast);interval<-match.arg(interval)
  result<-cs_causal_contrast(object,contrast=contrast,level=level,interval=interval,draws=draws,seed=seed)
  scope<-if(interval=="multiplier")"Simultaneous asymptotic intervals over the regular fitted contrast coordinates; no finite-sample guarantee." else
    if(interval=="none")"Confidence intervals not requested." else "Pointwise asymptotic causal-contrast intervals under the stated identification and regularity assumptions; no finite-sample guarantee."
  out<-.cs_confint_matrix(result$lower,result$upper,result$target,selected,level,result$inference_status,interval,scope)
  attr(out,"contrast")<-contrast;out
}

#' @export
tidy <- function(x,...)UseMethod("tidy")
#' @export
tidy.coherentSurv <- function(x,...)summary(x)$coefficients
#' @export
glance <- function(x,...)UseMethod("glance")
#' @export
glance.coherentSurv <- function(x,...)data.frame(n=x$n,logLik=x$loglik,AIC=-2*x$loglik+2*length(x$coefficients),frailty_variance=x$phi,converged=x$convergence$valid,gradient_max=x$convergence$gradient)
#' @export
augment <- function(x,...)UseMethod("augment")
#' @export
augment.coherentSurv <- function(x,newdata=NULL,landmark,horizon,...)predict(x,newdata,landmark,horizon,type="all",...)
#' @export
bounds <- function(object,...)UseMethod("bounds")
#' @export
bounds.default <- function(object,...)cs_bounds(p=object,...)
#' @export
bounds.coherentSurv <- function(object,...) .cs_error("bounds_schema","Recording-sensitivity bounds require an explicit finite-support recorded law. Use cs_bounds().")
#' @export
conformal <- function(object,...)cs_conformal(scores=object,...)
#' @export
plot.coherentSurv_prediction <- function(x,unit=unique(x$unit)[1],target="risk",...) {
  d<-x[x$unit==unit & x$target==target,,drop=FALSE]
  if(!nrow(d)) .cs_error("plot","No matching unit and target.")
  d<-d[order(d$horizon),];graphics::plot(d$horizon,d$estimate,type="b",xlab=paste("Horizon",d$time_unit[1]),ylab=target,...);invisible(x)
}
#' @export
autoplot <- function(object,...)UseMethod("autoplot")
#' @export
autoplot.coherentSurv_prediction <- function(object,unit=unique(object$unit)[1],target="risk",...) {
  if(!requireNamespace("ggplot2",quietly=TRUE)) .cs_error("optional_dependency","Install ggplot2 for autoplot(); base plot() is also available.")
  d<-object[object$unit==unit & object$target==target,,drop=FALSE]
  ggplot2::ggplot(d,ggplot2::aes(x=horizon,y=estimate))+ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::labs(x=paste("Horizon",d$time_unit[1]),y=target)
}
utils::globalVariables(c("horizon","estimate"))
