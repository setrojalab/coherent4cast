.cs_wrapper_formula <- function(formula,time,status,start) {
  if(!inherits(formula,"formula")) .cs_error("formula","Supply a formula.")
  if(length(formula)==3L) {
    lhs<-formula[[2L]]
    name<-if(is.call(lhs))paste(deparse(lhs[[1L]]),collapse="")else ""
    if(!name %in% c("Surv","survival::Surv") || !length(lhs) %in% c(3L,4L) ||
       !all(vapply(as.list(lhs)[-1L],is.symbol,logical(1))))
      .cs_error("survival_formula","Use Surv(time,status) or Surv(start,stop,status) with column names, or a one-sided formula and named columns.")
    cols<-vapply(as.list(lhs)[-1L],as.character,character(1))
    start<-if(length(cols)==3L)cols[1L]else NULL
    time<-cols[length(cols)-1L];status<-cols[length(cols)]
    formula<-stats::as.formula(call("~",formula[[3L]]),env=environment(formula))
  }
  if(length(formula)!=2L) .cs_error("formula","Use a one-sided or supported Surv formula.")
  if(any(all.vars(formula) %in% c(time,status,start)))
    .cs_error("formula_leakage","Outcome times and event labels cannot be predictors.")
  list(formula=formula,time=time,status=status,start=start)
}

.cs_wrapper_data <- function(data,spec,competing=FALSE,profiles=FALSE,landmark=0,horizon=1) {
  if(!is.data.frame(data)||!nrow(data)) .cs_error("data","Supply a nonempty data frame.")
  internal<-c(".rcw_unit",".rcw_start",".rcw_stop",".rcw_recurrent",".rcw_terminal",".rcw_cause")
  if(any(internal %in% names(data))) .cs_error("reserved_columns","Input columns beginning .rcw_ are reserved for wrapper mappings.")
  if(!all(all.vars(spec$formula)%in%names(data))) .cs_error("columns","A predictor column is missing.")
  available<-c(spec$time,spec$status) %in% names(data)
  profile_only<-profiles&&!any(available)
  if(!profile_only&&!all(available)) .cs_error("columns","Supply both follow-up and status columns, or covariate-only prediction profiles.")
  if(!is.null(spec$unit)&&!spec$unit%in%names(data)) .cs_error("unit","The declared unit column is missing.")
  if(!is.null(spec$start)&&is.null(spec$unit)) .cs_error("unit","Start-stop data require an explicit unit column.")
  d<-data
  d$.rcw_unit<-if(is.null(spec$unit))as.character(seq_len(nrow(d)))else as.character(d[[spec$unit]])
  if(profile_only) {
    if(anyDuplicated(d$.rcw_unit)) .cs_error("profiles","Use one covariate-only profile per eligible unit.")
    d$.rcw_start<-0;d$.rcw_stop<-landmark+max(horizon)+1
    code<-rep(as.character(spec$censor_code),nrow(d))
  } else {
    if(!is.null(spec$start)&&!spec$start%in%names(d)) .cs_error("columns","The start column is absent.")
    d$.rcw_start<-if(is.null(spec$start))rep(0,nrow(d))else d[[spec$start]]
    d$.rcw_stop<-d[[spec$time]];code<-as.character(d[[spec$status]])
    if(anyNA(code)) .cs_error("status","Event status cannot be missing.")
    if(!competing&&!all(code%in%c("0","1"))) .cs_error("status","Survival status must be 0=censored or 1=event.")
    if(competing&&is.factor(d[[spec$status]])&&!as.character(spec$censor_code)%in%levels(d[[spec$status]]))
      .cs_error("censor_code","For a factor cause column, explicitly match censor_code to its censoring level.")
    if(competing&&!all(code%in%c(as.character(spec$censor_code),spec$causes)))
      .cs_error("cause","A cause label is outside the declared causes.")
  }
  d$.rcw_recurrent<-"none";d$.rcw_terminal<-as.integer(code!=as.character(spec$censor_code))
  d$.rcw_cause<-ifelse(d$.rcw_terminal==1L,code,"0")
  out<-cs_as_trajectory(d,unit=".rcw_unit",start=".rcw_start",stop=".rcw_stop",
      event_type=".rcw_recurrent",terminal=".rcw_terminal",event_types=character(),
      time_unit=spec$time_unit,terminal_cause=if(competing)".rcw_cause"else NULL)
  attr(out,"profile_only")<-profile_only
  out
}

.cs_wrapper_dots <- function(dots) {
  if(length(dots)&&(is.null(names(dots))||any(!nzchar(names(dots))))) .cs_error("wrapper_arguments","Additional fit arguments must be named.")
  banned<-c("frailty","unit","start","stop","event_type","terminal","event_types","time_unit","targeting","observation","coarsening","subset","nonrecurrent")
  if(any(names(dots)%in%banned)) .cs_error("wrapper_scope","Schema, no-frailty, exact observation and ordinary estimation are fixed by these wrappers; use cs_fit() for other experiments.")
}

#' Fit right-censored survival using the existing no-frailty likelihood
#' @export
cs_surv <- function(formula,data,time="time",status="status",unit=NULL,start=NULL,
                    time_unit,breaks=c(0,Inf),control=cs_control(),...) {
  .cs_wrapper_dots(list(...))
  spec<-.cs_wrapper_formula(formula,time,status,start)
  spec$unit<-unit;spec$time_unit<-time_unit;spec$censor_code<-0
  d<-.cs_wrapper_data(data,spec)
  fit<-cs_fit(spec$formula,d,unit=".rcw_unit",start=".rcw_start",stop=".rcw_stop",
    event_type=".rcw_recurrent",terminal=".rcw_terminal",event_types=character(),
    time_unit=time_unit,breaks=breaks,control=control,frailty="none",...)
  fit$wrapper<-spec;fit$audit$wrapper_call<-match.call()
  class(fit)<-c("coherentSurv_surv",class(fit));fit
}

.cs_wrapper_future <- function(future,spec) {
  if(is.null(future))return(NULL)
  if(!is.data.frame(future)) .cs_error("future_covariates","Supply a data frame of future updates.")
  column<-if(is.null(spec$unit))"unit"else spec$unit
  if(!column%in%names(future)) .cs_error("future_covariates",paste("Future updates require unit column",column))
  future$.rcw_unit<-as.character(future[[column]]);future
}

#' @export
predict.coherentSurv_surv <- function(object,newdata=NULL,landmark=0,horizon,
                                  newdata_future=NULL,covariate_history=NULL,...) {
  d<-if(is.null(newdata))NULL else .cs_wrapper_data(newdata,object$wrapper,profiles=TRUE,landmark=landmark,horizon=horizon)
  out<-predict.coherentSurv(object,newdata=d,landmark=landmark,horizon=horizon,
    newdata_future=.cs_wrapper_future(newdata_future,object$wrapper),
    covariate_history=.cs_wrapper_future(covariate_history,object$wrapper),...)
  if(!is.null(d)&&isTRUE(attr(d,"profile_only"))) attr(out,"profile_assumption")<-"Each supplied covariate profile is declared eligible and event-free at the landmark; synthetic mapping intervals are not observed follow-up."
  out
}

#' Fit cause-specific no-frailty hazards and coherent cumulative incidence
#' @export
cs_competing <- function(formula,data,time="time",status="status",unit=NULL,start=NULL,
                         causes=NULL,censor_code=0,time_unit,breaks=c(0,Inf),
                         control=cs_control(),...) {
  dots<-list(...);.cs_wrapper_dots(dots)
  spec<-.cs_wrapper_formula(formula,time,status,start)
  if(length(censor_code)!=1L||is.na(censor_code)||!nzchar(as.character(censor_code))) .cs_error("censor_code","Declare one nonmissing censoring code.")
  if(!is.data.frame(data)||!spec$status%in%names(data)) .cs_error("columns","The status/cause column is missing.")
  code<-as.character(data[[spec$status]])
  if(is.null(causes))causes<-sort(unique(code[!is.na(code)&code!=as.character(censor_code)]))
  causes<-as.character(causes)
  if(!length(causes)||anyNA(causes)||any(!nzchar(causes))||anyDuplicated(causes)||as.character(censor_code)%in%causes||"0"%in%causes)
    .cs_error("causes","Declare distinct nonzero absorbing cause labels, excluding the censoring code; cause '0' is reserved for no event.")
  spec$unit<-unit;spec$time_unit<-time_unit;spec$censor_code<-censor_code;spec$causes<-causes
  d<-.cs_wrapper_data(data,spec,competing=TRUE)
  if(length(attr(d,"cs_schema")$delayed_entry)) .cs_error("left_truncation","These wrappers currently require histories from the time origin.")
  if(length(attr(d,"cs_schema")$gaps)) .cs_error("unrecorded_history","Competing-risk fitting requires contiguous observed histories.")
  fits<-stats::setNames(vector("list",length(causes)),causes)
  for(k in seq_along(causes)) {
    dk<-d;class(dk)<-"data.frame";attr(dk,"cs_schema")<-NULL
    dk$.rcw_terminal<-as.integer(d$.rcw_cause==causes[k]);dk$.rcw_cause<-ifelse(dk$.rcw_terminal==1L,causes[k],"0")
    if(!sum(dk$.rcw_terminal)) .cs_error("no_cause_events",paste("No events support the declared cause",causes[k],"; its hazard is not estimated by this regular wrapper."))
    fits[[k]]<-cs_fit(spec$formula,dk,unit=".rcw_unit",start=".rcw_start",stop=".rcw_stop",
      event_type=".rcw_recurrent",terminal=".rcw_terminal",event_types=character(),
      time_unit=time_unit,breaks=breaks,control=control,frailty="none",...)
  }
  n<-fits[[1]]$n
  if(any(vapply(fits,function(f)f$n!=n,logical(1)))) .cs_error("competing_alignment","Cause-specific fits must use identical eligible units.")
  all_names<-unlist(lapply(seq_along(fits),function(k)paste0("cause[",causes[k],"]:",names(fits[[k]]$coefficients))))
  coefs<-stats::setNames(unlist(lapply(fits,function(f)unname(f$coefficients))),all_names)
  sizes<-vapply(fits,function(f)length(f$coefficients),integer(1));H<-matrix(0,sum(sizes),sum(sizes));pos<-0L
  for(k in seq_along(fits)){ix<-pos+seq_len(sizes[k]);H[ix,ix]<-fits[[k]]$information;pos<-pos+sizes[k]}
  scores<-do.call(cbind,lapply(fits,`[[`,"scores"));weighted_scores<-scores*fits[[1]]$case_weights
  inv<-.cs_inverse(H);V<-inv%*%crossprod(weighted_scores)%*%inv;dimnames(V)<-list(all_names,all_names)
  converged<-all(vapply(fits,function(f)isTRUE(f$convergence$valid),logical(1)))
  out<-list(causes=causes,fits=fits,coefficients=coefs,covariance=V,information=H,scores=scores,
    loglik=sum(vapply(fits,`[[`,numeric(1),"loglik")),n=n,data=d,wrapper=spec,
    schema=attr(d,"cs_schema"),breaks=breaks,control=fits[[1]]$control,
    convergence=list(valid=converged,code=if(converged)"converged"else "nonregular_cause_fit",
                     causes=lapply(fits,`[[`,"convergence")),
    audit=c(list(call=match.call(),assumptions=c("cause-specific piecewise-exponential proportional hazards",
      "independent units and valid censoring observation model","mutually exclusive absorbing first-event causes",
      "no frailty; cumulative incidence integrates every cause-specific hazard")),.cs_capture_audit_context(d,match.call())))
  class(out)<-"coherentSurv_competing";out
}

#' @export
predict.coherentSurv_competing <- function(object,newdata=NULL,landmark=0,horizon,
    type=c("cif","survival","risk","rmst","all"),covariate_evolution=c("freeze","locf"),
    newdata_future=NULL,...) {
  if(length(list(...))) .cs_error("unused_argument","Unsupported competing-risk prediction argument.")
  type<-match.arg(type);mode<-match.arg(covariate_evolution)
  if(length(landmark)!=1L||!is.finite(landmark)||landmark<0||!is.numeric(horizon)||!length(horizon)||any(!is.finite(horizon)|horizon<=0))
    .cs_error("horizon","Declare a nonnegative landmark and positive finite durations.")
  d<-if(is.null(newdata))object$data else .cs_wrapper_data(newdata,object$wrapper,TRUE,TRUE,landmark,horizon)
  sc<-attr(d,"cs_schema")
  if(length(sc$delayed_entry)||length(sc$gaps)||any(!d$.cs_observed)) .cs_error("observation","Prediction requires contiguous eligible histories from the origin.")
  future<-.cs_wrapper_future(newdata_future,object$wrapper)
  base<-object$fits[[1]];C<-length(object$causes);output<-list();excluded<-character()
  for(rows in split(seq_len(nrow(d)),d$.cs_unit)) {
    ident<-as.character(d$.rcw_unit[rows[1]])
    if(any(d$.rcw_terminal[rows]==1&d$.rcw_stop[rows]<=landmark)||max(d$.rcw_stop[rows])<landmark){excluded<-c(excluded,ident);next}
    current<-tail(rows[d$.rcw_start[rows]<=landmark],1)
    if(!length(current)){excluded<-c(excluded,ident);next}
    for(delta in horizon) {
      scenario<-.cs_future_scenarios(base,d,rows,current,ident,landmark,delta,mode,future,NULL,1L,1L,1L)
      generators<-vector("list",length(scenario$cuts)-1L)
      for(j in seq_along(generators)) {
        bin<-findInterval(scenario$cuts[j],object$breaks,all.inside=TRUE)
        baseline<-numeric(base$design$base_size);baseline[bin]<-1
        v<-c(baseline,scenario$covariates[[1]][j,])
        rates<-vapply(object$fits,function(f)exp(sum(v*.cs_linear_coef(f))),numeric(1))
        if(any(!is.finite(rates))) .cs_error("forecast_overflow","Cause-specific rates overflow; rescale predictors.")
        generators[[j]]<-cs_generator_csr(C+1L,rep(1L,C),seq_len(C)+1L,rates)
      }
      result<-cs_propagate_sparse(c(1,rep(0,C)),generators,diff(scenario$cuts),
        reward=matrix(c(1,rep(0,C)),ncol=1),tolerance=object$control$propagation_tolerance)
      values<-c(stats::setNames(result$probability[-1L],paste0("cif_",object$causes)),
        survival=result$probability[1],risk=sum(result$probability[-1L]),rmst=result$reward[1])
      wanted<-switch(type,cif=paste0("cif_",object$causes),all=names(values),type)
      error<-stats::setNames(rep(result$error_bound,C+3L),names(values));error["rmst"]<-result$reward_error_bound[1]
      error["risk"]<-C*result$error_bound
      output[[length(output)+1L]]<-data.frame(unit=ident,landmark=landmark,horizon=delta,target=wanted,
        estimate=unname(values[wanted]),time_unit=object$schema$time_unit,
        series_error_bound=unname(error[wanted]),probability_roundoff_bound=result$roundoff_bound,
        covariate_evolution=mode,stringsAsFactors=FALSE)
    }
  }
  if(!length(output)) .cs_error("eligibility","No unit is event-free and observed at the landmark.")
  out<-do.call(rbind,output);rownames(out)<-NULL;class(out)<-c("coherentSurv_prediction","cs_forecast","data.frame")
  attr(out,"excluded_units")<-excluded
  attr(out,"assumptions")<-c(object$audit$assumptions,"Forecast numerical bounds exclude sampling and model error.")
  if(isTRUE(attr(d,"profile_only")))attr(out,"profile_assumption")<-"Supplied profiles are declared event-free at the landmark; mapping intervals are not observed follow-up."
  out
}

#' @export
print.coherentSurv_competing <- function(x,...) {cat("coherentSurv cause-specific no-frailty competing-risk model\n",x$n,"independent units; causes:",paste(x$causes,collapse=", "),"\n");invisible(x)}
#' @export
coef.coherentSurv_competing <- function(object,...)object$coefficients
#' @export
vcov.coherentSurv_competing <- function(object,...)object$covariance
#' @export
logLik.coherentSurv_competing <- function(object,...)structure(object$loglik,df=length(object$coefficients),nobs=object$n,class="logLik")
#' @export
tidy.coherentSurv_competing <- function(x,...) {
  se<-sqrt(pmax(0,diag(x$covariance)))
  regular<-isTRUE(x$convergence$valid)&&length(x$fits)>0L&&
    all(vapply(x$fits,function(f)isTRUE(f$convergence$valid),logical(1)))
  lower<-upper<-rep(NA_real_,length(x$coefficients))
  if(regular){lower<-unname(x$coefficients)-stats::qnorm(.975)*se;upper<-unname(x$coefficients)+stats::qnorm(.975)*se}
  data.frame(term=names(x$coefficients),estimate=unname(x$coefficients),std.error=se,
    conf.low=lower,conf.high=upper,
    inference_status=if(regular)"regular_asymptotic" else "unavailable_nonregular_fit")
}
#' @export
summary.coherentSurv_competing <- function(object,...) {
  coefficients<-tidy(object)
  regular<-all(coefficients$inference_status=="regular_asymptotic")
  list(coefficients=coefficients,convergence=object$convergence,
    inference_status=if(regular)"regular_asymptotic" else "unavailable_nonregular_fit",
    interval_scope=if(regular)"Regular joint sandwich approximation; no automatic CIF confidence bands." else "Regular joint sandwich intervals withheld: the wrapper or a cause-specific fit did not meet convergence and regularity criteria; no automatic CIF confidence bands. Standard errors are numerical diagnostics only.",audit=object$audit)
}
#' @export
glance.coherentSurv_competing <- function(x,...)data.frame(n=x$n,causes=length(x$causes),logLik=x$loglik,
    AIC=-2*x$loglik+2*length(x$coefficients),converged=x$convergence$valid)

#' Validate cause-specific cumulative incidence without censoring competing events
#' @export
validate.coherentSurv_competing <- function(object,time,status,horizon,cause,predicted,
    censor_survival=NULL,threshold=.2,strata=NULL,...) {
  if(length(list(...))) .cs_error("unused_argument","Unsupported competing-risk validation argument.")
  if(length(cause)!=1L||!as.character(cause)%in%object$causes) .cs_error("cause","Select one fitted absorbing cause.")
  n<-length(time);code<-as.character(status);cc<-as.character(object$wrapper$censor_code)
  if(!is.numeric(time)||!n||length(code)!=n||length(predicted)!=n||any(!is.finite(time))||any(time<0)||anyNA(code)||!all(code%in%c(cc,object$causes)))
    .cs_error("validation","Provide aligned nonnegative follow-up, declared cause labels and risks.")
  if(length(horizon)!=1L||!is.finite(horizon)||horizon<=0) .cs_error("horizon","Supply one finite positive validation horizon.")
  if(is.null(censor_survival)) {
    times<-sort(unique(time));at<-vapply(times,function(t)sum(time>=t),numeric(1))
    dc<-vapply(times,function(t)sum(time==t&code==cc),numeric(1));surv<-cumprod(1-dc/at)
    ge<-vapply(time,function(t){j<-which(times<t);if(length(j))surv[max(j)]else 1},numeric(1))
    j<-which(times<=horizon);gh<-rep(if(length(j))surv[max(j)]else 1,n)
    assumption<-"Marginal independent right censoring; competing events are observed noncases for this cumulative-incidence endpoint."
  } else {
    ge<-censor_survival$event;gh<-censor_survival$horizon
    if(length(ge)!=n||length(gh)!=n) .cs_error("censoring","Supply event=G(T_i-) for every failure and horizon=G(h), each of length n.")
    assumption<-"Supplied conditional censoring probabilities must be valid and correctly aligned; competing events are not treated as censoring."
  }
  if(any(!is.finite(c(ge,gh)))||any(c(ge,gh)<0|c(ge,gh)>1)) .cs_error("censoring","Censoring-survival values must lie in [0,1].")
  other<-code!=cc&code!=as.character(cause)&time<=horizon
  mapped_time<-time;mapped_time[other]<-horizon+max(1,horizon)*.5
  if(any(!is.finite(mapped_time))) .cs_error("horizon","Rescale the validation time unit to avoid overflow.")
  gh[other]<-ge[other]
  result<-validate.default(predicted,mapped_time,as.integer(code==as.character(cause)),horizon,
    censor_survival=list(case=ge,horizon=gh),threshold=threshold,strata=strata)
  annotate<-function(x){x$assumptions<-assumption;x$cause<-as.character(cause);x$endpoint<-"Cause-specific cumulative incidence by horizon; earlier competing events are observed noncases.";x}
  if(is.null(strata)){result<-annotate(result);result$competing_events_by_horizon<-sum(other)}else result<-lapply(result,annotate)
  result
}
