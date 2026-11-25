.rc_target_rng <- function(seed, code) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed); force(code)
}
.rc_target_spec_value <- function(spec, name, default) if (is.null(spec[[name]])) default else spec[[name]]

.rc_targeted_fit <- function(formula, data, unit="id", start="start", stop="stop",
 event_type="event_type", terminal="terminal", observed=NULL, event_types=NULL,
 time_unit, targeting, control, targeting_spec, seed=1L, call=NULL) {
  spec <- targeting_spec
  if(!inherits(formula,"formula")||length(formula)!=2L) .rc_error("targeting_formula","Supply a one-sided baseline covariate formula.")
  if(missing(time_unit)||length(time_unit)!=1L||is.na(time_unit)||!nzchar(time_unit)) .rc_error("time_unit","Declare one nonempty time unit.")
  if (!is.list(spec)) .rc_error("targeting_spec", "Supply an explicit targeting_spec list with the estimand and observation schema.")
  if (targeting == "finite_path") {
    if(is.null(spec$strata))spec$strata<-all.vars(formula)
    if(!all(all.vars(formula)%in%spec$strata)) .rc_error("targeting_strata","Every formula predictor must enter the finite-path baseline strata.")
    required <- c("path", "recorded", "features")
    if (!all(required %in% names(spec))) .rc_error("targeting_spec", "Finite-path targeting requires path and recorded column names and a feature matrix.")
    columns <- c(spec$path, spec$recorded, spec$strata)
    if (!all(columns %in% names(data))) .rc_error("targeting_spec", "A finite-path schema column is absent.")
    if (unit %in% names(data) && anyDuplicated(data[[unit]])) .rc_error("targeting_subjects", "The finite-path experiment needs one row per independent subject.")
    strata <- if (!length(spec$strata)) rep(1L,nrow(data)) else do.call(interaction,c(data[,spec$strata,drop=FALSE],list(drop=TRUE,lex.order=TRUE)))
    folds <- spec$folds
    if (is.character(folds) && length(folds)==1L && folds %in% names(data)) folds <- data[[folds]]
    fit <- rc_target(data[[spec$path]],data[[spec$recorded]],strata,spec$features,folds,
                     spec$initial_laws,spec$recording,.rc_target_spec_value(spec,"control",list()))
    target_names <- colnames(as.matrix(spec$features));if(is.null(target_names))target_names<-paste0("target_",seq_along(fit$estimate))
    fit$target_names<-target_names;names(fit$estimate)<-target_names
    fit$target_limits<-rbind(apply(as.matrix(spec$features),2,min),apply(as.matrix(spec$features),2,max))
    feature_matrix<-as.matrix(spec$features);recorded_values<-feature_matrix[data[[spec$path]][data[[spec$recorded]]==1],,drop=FALSE]
    fit$zero_event_boundary<-vapply(seq_len(ncol(feature_matrix)),function(j)all(recorded_values[,j]==fit$target_limits[1,j])||all(recorded_values[,j]==fit$target_limits[2,j]),logical(1))
    lower_observed<-matrix(rep(fit$target_limits[1,],each=nrow(data)),nrow(data));upper_observed<-matrix(rep(fit$target_limits[2,],each=nrow(data)),nrow(data))
    known<-which(data[[spec$recorded]]==1);lower_observed[known,]<-upper_observed[known,]<-recorded_values
    fit$observation_bounds<-list(lower=lower_observed,upper=upper_observed)

    fit$convergence$score_solved<-isTRUE(fit$convergence$valid)
    if(any(fit$zero_event_boundary)){fit$convergence$valid<-FALSE;fit$convergence$code<-"zero_event_boundary"}

    fit$targeting<-targeting;fit$n<-nrow(data);fit$data<-data;fit$formula<-formula;fit$targeting_spec<-spec
    fit$horizon<-.rc_target_spec_value(spec,"horizon",1);fit$time_unit<-time_unit
    if(isTRUE(spec$complete_followup)){if(any(data[[spec$recorded]]!=1)) .rc_error("complete_followup","The complete-followup declaration conflicts with missing paths.");fit$bounded_data<-as.matrix(spec$features)[data[[spec$path]],,drop=FALSE]}
    fit$audit<-list(call=call,seed=seed,assumptions=fit$assumptions);class(fit)<-c("coherent4cast_targeted","coherent4cast")
    return(fit)
  }
  if (targeting != "continuous") .rc_error("targeting", "Unknown targeting route.")
  tau<-spec$horizon;.rc_scalar(tau,"target horizon",0,Inf,TRUE)
  if(is.null(spec$censor)||length(spec$censor)!=1L||!spec$censor %in% names(data)) .rc_error("censor_schema","Continuous targeting requires a named right-censor indicator column.")
  d<-as_trajectory(data,unit,start,stop,event_type,terminal,observed,event_types,time_unit,"exact")
  schema<-attr(d,"rc_schema")
  if(length(schema$gaps)||length(schema$delayed_entry)||any(!d$.rc_observed)) .rc_error("targeting_observation","Continuous targeting currently requires contiguous histories from time zero, with right censoring explicitly marked.")
  if(anyNA(d[[spec$censor]])||!all(d[[spec$censor]] %in% c(0,1))) .rc_error("censor_schema","Right-censor indicators must be binary.")
  strata_names<-.rc_target_spec_value(spec,"strata",all.vars(formula))
  if(!all(strata_names %in% names(d))) .rc_error("targeting_strata","A baseline stratum column is absent.")
  if(!all(all.vars(formula) %in% strata_names)) .rc_error("targeting_strata","Every formula predictor must enter the declared baseline strata; predictors are not silently discarded.")
  subjects<-split(seq_len(nrow(d)),d$.rc_unit);n<-length(subjects)
  if(n<4L) .rc_error("targeting_subjects","Continuous targeting requires at least four independent subjects.")
  for(ix in subjects){
    if(any(d[[spec$censor]][ix[-length(ix)]]!=0)||any(d[[spec$censor]][ix]==1 & (d$.rc_mark[ix]>0|d$.rc_terminal[ix]==1))) .rc_error("censor_schema","Right censoring must terminate a subject's record and cannot coincide with an event.")
    last<-tail(ix,1);if(d$.rc_stop[last]<tau && d$.rc_terminal[last]==0 && d[[spec$censor]][last]!=1) .rc_error("censor_schema","A history ending before the horizon needs an explicit terminal event or censoring indicator.")
    if(length(strata_names)&&any(vapply(d[ix,strata_names,drop=FALSE],function(z)length(unique(z))!=1L,logical(1)))) .rc_error("targeting_strata","This finite-state targeting implementation uses baseline discrete strata; changing predictors require an expanded declared state.")
  }
  if(isTRUE(spec$complete_followup)&&any(d[[spec$censor]]==1 & d$.rc_stop<tau)) .rc_error("complete_followup","The declared uncensored design contains right-censoring events.")
  first<-vapply(subjects,`[`,integer(1),1)
  if(length(strata_names)) {keys<-do.call(paste,c(d[first,strata_names,drop=FALSE],sep="\r"));levels<-sort(unique(keys));x<-match(keys,levels)} else {keys<-rep("all",n);levels<-"all";x<-rep(1L,n)}
  nx<-length(levels);if(nx>30L) .rc_error("targeting_strata","At most 30 prespecified baseline strata are supported; discretize scientific strata explicitly.")
  folds<-spec$folds;if(is.character(folds)&&length(folds)==1L&&folds %in% names(d))folds<-d[[folds]][first]
  nf<-.rc_target_spec_value(spec,"nfolds",5L);.rc_scalar(nf,"number of folds",2,20)
  if(is.null(folds))folds<-.rc_target_rng(seed,sample(rep(seq_len(min(nf,n)),length.out=n)))
  if(length(folds)!=n||anyNA(folds)) .rc_error("targeting_folds","One fold label is required per independent subject.")
  folds<-match(folds,sort(unique(folds)));nf<-max(folds);if(nf<2L) .rc_error("targeting_folds","At least two folds are required.")
  L<-.rc_target_spec_value(spec,"grid",64L);.rc_scalar(L,"targeting time grid",4,512);if(L!=floor(L)) .rc_error("targeting_grid","grid must be integer.")
  cap<-.rc_target_spec_value(spec,"history_cap",control$history_cap);.rc_scalar(cap,"retained history cap",0,6)
  K<-length(schema$event_types);M<-K+1L;J<-K+2L
  if(K<1L||K>6L||cap!=floor(cap)||(cap+1)^K>128) .rc_error("targeting_state","Declare one to six recurrent marks and at most 128 sufficient living states.")
  states<-as.matrix(expand.grid(rep(list(0:cap),K)))[,rev(seq_len(K)),drop=FALSE];S<-nrow(states)
  powers<-(cap+1)^rev(seq.int(0,K-1L));nxt<-matrix(S,S,M)
  for(q in seq_len(S))for(k in seq_len(K)){s<-states[q,];s[k]<-min(cap,s[k]+1L);nxt[q,k]<-sum(s*powers)}
  storage.mode(nxt)<-"integer"
  history<-build_history(d,cap=cap)$counts
  q<-as.integer(as.vector(history %*% powers))+1L
  # Restrict the experiment to the declared horizon without treating later events as observed.
  keep<-d$.rc_start<tau;dd<-d[keep,,drop=FALSE];q<-q[keep];original_end<-dd$.rc_stop
  dd$.rc_stop<-pmin(dd$.rc_stop,tau);after<-original_end>tau
  dd$.rc_mark[after]<-0L;dd$.rc_terminal[after]<-0L;dd[[spec$censor]][after]<-0L
  u<-dd$.rc_unit;group<-(folds-1L)*nx+x;G<-nf*nx
  event<-ifelse(dd$.rc_terminal==1,K,ifelse(dd$.rc_mark>0,dd$.rc_mark-1L,ifelse(dd[[spec$censor]]==1,M,-1L)))
  records<-cbind(u-1L,group[u]-1L,q-1L,dd$.rc_start,dd$.rc_stop,event)
  exposure<-dd$.rc_stop-dd$.rc_start
  nonrecurrent<-match(.rc_target_spec_value(spec,"nonrecurrent",character()),schema$event_types)
  if(anyNA(nonrecurrent)) .rc_error("targeting_eligibility","An ineligible mark label is absent from event_types.")
  for(k in nonrecurrent)if(any(dd$.rc_mark==k & history[keep,k]>0)) .rc_error("targeting_eligibility","A non-repeatable mark recurs in the supplied history.")
  rates<-array(0,c(M,S,L,G));censor<-numeric(G);training<-vector("list",nf)
  working<-.rc_target_spec_value(spec,"working","loglinear");if(!working %in% c("loglinear","pooled")) .rc_error("targeting_working","working must be loglinear or pooled.")
  offsets<-spec$rate_offset;if(!is.null(offsets)&&(!is.matrix(offsets)||!identical(dim(offsets),c(M,S))||any(!is.finite(offsets))||any(offsets<=0))) .rc_error("targeting_offset","rate_offset must be a positive mark-by-state matrix.")
  for(f in seq_len(nf)){
    training[[f]]<-which(folds!=f)
    for(a in seq_len(nx)){
      g<-(f-1L)*nx+a;ii<-which(folds[u]!=f & x[u]==a)
      if(!length(ii)||sum(exposure[ii])<=0) .rc_error("targeting_training","A training fold has no exposure in a declared stratum.")
      risk_time<-sum(exposure[ii]);censor[g]<-(sum(event[ii]==M)+.5)/(risk_time+1)
      for(k in seq_len(M)){
        if(!is.null(offsets))lambda<-(sum(event[ii]==k-1L)+.5)/(sum(exposure[ii]*offsets[k,q[ii]])+1)*offsets[k,] else if(working=="pooled")lambda<-rep((sum(event[ii]==k-1L)+.5)/(risk_time+1),S) else {
          # tabulate weights without silently dropping unvisited retained states.
          tt<-vapply(seq_len(S),function(s)sum(exposure[ii[q[ii]==s]]),numeric(1))
          cc<-tabulate(q[ii[event[ii]==k-1L]],S)+.5/S;tt<-tt+1/S
          X<-cbind(1,states);init<-c(log(sum(cc)/sum(tt)),rep(0,K))
          fn<-function(b){eta<-as.vector(X%*%b);sum(tt*exp(eta)-cc*eta)+1e-6*sum(b[-1]^2)}
          gr<-function(b)as.vector(crossprod(X,tt*exp(as.vector(X%*%b))-cc))+c(0,2e-6*b[-1])
          opt<-stats::optim(init,fn,gr,method="L-BFGS-B",lower=rep(-10,K+1),upper=rep(5,K+1),control=list(maxit=200))
          lambda<-exp(as.vector(X%*%opt$par))
        }
        if(k %in% nonrecurrent)lambda[states[,k]>0]<-0
        for(b in seq_len(L))rates[k,,b,g]<-lambda
      }
    }
  }
  if(isTRUE(spec$complete_followup))censor[]<-0
  if(!is.null(spec$initial_rates)){if(!identical(dim(spec$initial_rates),dim(rates))) .rc_error("targeting_initial","initial_rates dimensions must be marks by states by grid by fold-stratum.");rates<-spec$initial_rates}
  if(!is.null(spec$censor_rates)){if(length(spec$censor_rates)!=G) .rc_error("targeting_censor","censor_rates need one value per fold-stratum.");censor<-spec$censor_rates}
  minimum_G<-exp(-max(censor)*tau);if(minimum_G<.rc_target_spec_value(spec,"min_censor_survival",.001)) .rc_error("targeting_positivity","Estimated censor survival is below the declared minimum; no clipping is applied.")
  if(isTRUE(spec$complete_followup)&&any(censor!=0)) .rc_error("complete_followup","A declared uncensored target has censor hazard zero; supplied censor_rates contradict it.")
  initial_rates<-rates
  out<-.rc_core_call(rc_cpp_continuous_target(as.numeric(rates),censor,nxt,as.integer(group-1L),records,as.integer(L),tau,
      as.integer(.rc_target_spec_value(spec,"maxit",60L)),.rc_target_spec_value(spec,"tolerance",1e-7/sqrt(n))))
  targets<-c(paste0("burden_",schema$event_types),"risk","rmst");names(out$estimate)<-names(out$se)<-targets;colnames(out$EIF)<-targets;dimnames(out$covariance)<-list(targets,targets)
  out$rates<-array(out$rates,dim(initial_rates));out$values<-array(out$values,c(J,S+1L,L+1L,G))
  out$EIF_mean<-colMeans(out$EIF);out$EIF_variance<-apply(out$EIF,2,stats::var)
  event_counts<-c(vapply(seq_len(K),function(k)sum(dd$.rc_mark==k),numeric(1)),sum(dd$.rc_terminal))
  out$observed_event_counts<-stats::setNames(event_counts,c(paste0("burden_",schema$event_types),"risk"))
  out$zero_event_boundary<-c(event_counts==0,event_counts[K+1L]==0)
  out$convergence<-list(valid=isTRUE(out$converged)&&all(out$se>0)&&!any(out$zero_event_boundary),score_solved=isTRUE(out$converged),code=if(!out$converged)"target_score_not_solved" else if(any(out$zero_event_boundary))"zero_event_boundary" else if(any(out$se<=0))"nonregular_variance" else "converged",score=out$score,iterations=out$iterations)
  out$target_names<-targets;out$target_limits<-rbind(rep(0,J),c(rep(Inf,K),1,tau))
  if(length(nonrecurrent))out$target_limits[2,nonrecurrent]<-1
  lower_observed<-upper_observed<-matrix(0,n,J)
  for(i in seq_len(n)){ix<-which(dd$.rc_unit==i);last<-tail(ix,1);censored_early<-dd[[spec$censor]][last]==1&&dd$.rc_stop[last]<tau
    lower_observed[i,]<-c(vapply(seq_len(K),function(k)sum(dd$.rc_mark[ix]==k),numeric(1)),sum(dd$.rc_terminal[ix]),dd$.rc_stop[last])
    upper_observed[i,]<-lower_observed[i,]
    if(censored_early){upper_observed[i,K+1L]<-1;upper_observed[i,J]<-tau;upper_observed[i,seq_len(K)]<-Inf;if(length(nonrecurrent))upper_observed[i,nonrecurrent]<-1}
  }
  out$observation_bounds<-list(lower=lower_observed,upper=upper_observed)

  if(isTRUE(spec$complete_followup)){out$bounded_data<-matrix(0,n,J);for(i in seq_len(n)){ix<-which(dd$.rc_unit==i);out$bounded_data[i,]<-c(vapply(seq_len(K),function(k)sum(dd$.rc_mark[ix]==k),numeric(1)),sum(dd$.rc_terminal[ix]),max(dd$.rc_stop[ix]))}}
  out$targeting<-targeting;out$n<-n;out$data<-d;out$formula<-formula;out$targeting_spec<-spec;out$horizon<-tau;out$time_unit<-time_unit
  out$folds<-folds;out$subject_stratum<-x;out$strata<-levels;out$strata_names<-strata_names;out$training_subjects<-training;out$censor_rates<-censor;out$initial_rates<-initial_rates
  out$states<-states;out$next_state<-nxt;out$grid<-L;out$schema<-schema;out$control<-control;out$group<-group;out$nfolds<-nf
  out$fold_prior<-matrix(tabulate(group,G),nx,nf);out$fold_prior<-out$fold_prior/rowSums(out$fold_prior)
  out$assumptions<-c("independent subjects; no simultaneous distinct jumps","right censoring conditionally independent given the declared baseline strata","strictly positive censor survival","retained count state sufficient for event intensity","cross-fitted initial nuisances; regular targeting and product remainder conditions for EIF inference","finite assessment vector; controlled time-grid interpolation error")
  out$audit<-list(call=call,seed=seed,assumptions=out$assumptions,numerical_error=out$interpolation_error,minimum_censor_survival=minimum_G,working_model=working)
  class(out)<-c("coherent4cast_targeted","coherent4cast")
  if(!out$convergence$valid)warning(paste("Targeting score, variance, or event-count boundary did not meet regular inference criteria:",out$convergence$code,"; point results and diagnostics are retained."),call.=FALSE)
  out
}

.rc_target_intervals <- function(object,method="logit_gaussian",level=.95,draws=1000L,seed=1L,allow_projection=TRUE){
 .rc_scalar(level,"confidence level",0,1,TRUE);.rc_scalar(draws,"multiplier draws",50,100000)
 estimate<-as.numeric(object$estimate);J<-length(estimate);n<-nrow(object$EIF)
 if(method=="projection")return(.rc_target_projection(object,level))
 if(!isTRUE(object$convergence$valid)||any(!is.finite(object$se))||any(object$se<=0)) {
  if(allow_projection&&!is.null(object$observation_bounds))return(.rc_target_projection(object,level))
  code<-if(identical(object$convergence$code,"zero_event_boundary"))"zero_event_boundary_requires_observation_model_projection" else "target_score_or_variance_not_regular"
  return(list(lower=rep(NA_real_,J),upper=rep(NA_real_,J),critical=NA_real_,status=code))
 }
 if(method=="multiplier"){center<-estimate;derivative<-rep(1,J);inverse<-function(z)z} else {
  lo<-object$target_limits[1,];hi<-object$target_limits[2,];bounded<-is.finite(lo)&is.finite(hi);positive<-!bounded&is.finite(lo)
  if(any(bounded&(estimate<=lo|estimate>=hi))||any(positive&estimate<=lo))return(list(lower=rep(NA_real_,J),upper=rep(NA_real_,J),critical=NA_real_,status="target_at_transform_boundary"))
  center<-estimate;derivative<-rep(1,J)
  center[bounded]<-stats::qlogis((estimate[bounded]-lo[bounded])/(hi[bounded]-lo[bounded]));derivative[bounded]<-(hi[bounded]-lo[bounded])/((estimate[bounded]-lo[bounded])*(hi[bounded]-estimate[bounded]))
  center[positive]<-log(estimate[positive]-lo[positive]);derivative[positive]<-1/(estimate[positive]-lo[positive])
  inverse<-function(z){z[bounded]<-lo[bounded]+(hi[bounded]-lo[bounded])*stats::plogis(z[bounded]);z[positive]<-lo[positive]+exp(z[positive]);z}
 }
 centered<-sweep(object$EIF,2,colMeans(object$EIF));se<-object$se*abs(derivative)
 perturbation<-.rc_target_rng(seed,matrix(stats::rnorm(draws*n),draws,n)%*%sweep(centered,2,derivative,"*")/sqrt(n*(n-1)))
 critical<-as.numeric(stats::quantile(apply(abs(sweep(perturbation,2,se,"/")),1,max),level,type=1))
 list(lower=pmax(object$target_limits[1,],inverse(center-critical*se)),upper=pmin(object$target_limits[2,],inverse(center+critical*se)),critical=critical,status="regular_asymptotic_approximation")
}

.rc_targeted_predict <- function(object,newdata=NULL,landmark=0,horizon=NULL,type="all",interval="auto",level=.95,draws=1000L,seed=1L){
 if(is.null(horizon))horizon<-object$horizon
 if(!interval %in% c("auto","none","multiplier","logit_gaussian","projection")) .rc_error("target_interval","Choose none, multiplier, logit_gaussian, or projection.")
 .rc_scalar(landmark,"landmark",0,Inf);if(!is.numeric(horizon)||!length(horizon)||any(!is.finite(horizon))||any(horizon<=0)) .rc_error("target_horizon","Forecast horizons must be positive and finite.")
 population<-is.null(newdata)&&landmark==0&&length(horizon)==1L&&abs(horizon-object$horizon)<1e-12
 if(population){
  method<-if(interval=="auto")"logit_gaussian" else interval
  ci<-if(method=="none")list(lower=rep(NA_real_,length(object$estimate)),upper=rep(NA_real_,length(object$estimate)),status="not_requested",critical=NA_real_) else .rc_target_intervals(object,method,level,draws,seed,allow_projection=interval%in%c("auto","projection"))
  wanted<-seq_along(object$estimate)
  if(type!="all"){wanted<-if(type=="burden")grep("^burden_",object$target_names) else which(object$target_names==type);if(type=="survival"){wanted<-which(object$target_names=="risk");if(length(wanted)){value<-1-object$estimate[wanted];lower<-1-ci$upper[wanted];upper<-1-ci$lower[wanted]}}}
  if(!length(wanted)) .rc_error("target_type","The requested target is not in this fitted target vector.")
  reported_method<-if(is.null(ci$method))method else ci$method
  statuses<-if(length(ci$status)==1L)rep(ci$status,length(wanted)) else ci$status[wanted]
  ans<-data.frame(unit="population",landmark=0,horizon=object$horizon,target=object$target_names[wanted],estimate=unname(object$estimate[wanted]),std_error=unname(object$se[wanted]),lower=ci$lower[wanted],upper=ci$upper[wanted],interval=reported_method,inference_status=statuses,time_unit=object$time_unit)
  if(type=="survival"){ans$target<-"survival";ans$estimate<-value;ans$lower<-lower;ans$upper<-upper}
  if(!isTRUE(object$convergence$valid))ans$std_error<-NA_real_
  attr(ans,"critical_value")<-ci$critical
  attr(ans,"requested_interval")<-interval
  attr(ans,"inference_scope")<-if(identical(reported_method,"bounded_coarsening_projection")) {
    if(interval=="auto")"Observed-reward coarsening outer-set fallback. The data-dependent automatic Gaussian/projection switch does not inherit an unconditional or selection-conditional finite-sample guarantee. Request interval='projection' throughout to use that construction's unconditional theorem; uncapped burdens retain [0,Inf]."
    else "Finite-sample Hoeffding-Bonferroni coarsening outer sets when this construction is used throughout, for prespecified bounded targets under independent subjects; uncapped burdens retain [0,Inf]."
  } else "Simultaneous finite-vector asymptotic EIF inference; no finite-sample sparse-coverage guarantee."
 } else {
  if(object$targeting!="continuous") .rc_error("finite_path_query","This finite-path target is a fitted population functional; dynamic conditional queries require an explicit path-query specification.")
  if(any(landmark+horizon>object$horizon+1e-12)) .rc_error("target_horizon","A targeted intensity law is estimated only through its prespecified horizon.")
  if(!interval %in% c("auto","none")) .rc_error("conditional_inference","Population EIFs do not provide individual-history confidence intervals; request interval='none' for conditional forecasts.")
  ans<-.rc_target_conditional(object,if(is.null(newdata))object$data else newdata,landmark,horizon,type)
  attr(ans,"inference_scope")<-"Conditional forecasts use the common targeted law with Bayes-updated fold weights. Population EIF intervals are not individual-history intervals."
 }
 class(ans)<-c("coherent4cast_prediction","rc_forecast","data.frame");ans
}

.rc_target_conditional <- function(object,data,landmark,horizon,type){
 if(!is.character(type)||length(type)!=1L||is.na(type)||!type%in%c("all","burden","survival",object$target_names)) .rc_error("target_type","Choose a target in the fitted target vector, burden, survival, or all.")
 sc<-object$schema;d<-as_trajectory(data,sc$unit,sc$start,sc$stop,sc$event_type,sc$terminal,sc$observed,sc$event_types,sc$time_unit,"exact")
 if(length(attr(d,"rc_schema")$delayed_entry)) .rc_error("target_history","Conditional targeted histories must begin at time zero; delayed entry requires an additional observation model.")
 if(length(attr(d,"rc_schema")$gaps)||any(!d$.rc_observed)) .rc_error("target_history","Conditional targeted forecasting needs a contiguous observed history.")
 K<-length(sc$event_types);S<-nrow(object$states);M<-K+1L;L<-object$grid;h<-object$horizon/L;nx<-length(object$strata);nf<-object$nfolds;rowsout<-list();weightsout<-list()
 powers<-(max(object$states)+1)^rev(seq.int(0,K-1L));cap<-max(object$states)
 for(ix in split(seq_len(nrow(d)),d$.rc_unit)){
  if(any(d$.rc_terminal[ix]==1 & d$.rc_stop[ix]<=landmark)||max(d$.rc_stop[ix])<landmark)next
  known<-ix[d$.rc_start[ix]<=landmark]
  if(length(object$strata_names)&&any(vapply(d[known,object$strata_names,drop=FALSE],function(z)length(unique(z))!=1L,logical(1)))) .rc_error("target_stratum","The fitted target uses baseline strata; observed changes require an expanded retained state rather than silently freezing the first value.")
  key<-if(length(object$strata_names))do.call(paste,c(d[ix[1],object$strata_names,drop=FALSE],sep="\r")) else "all";x<-match(key,object$strata)
  if(is.na(x)) .rc_error("target_stratum","A prediction stratum was not in the fitted population.")
  logweight<-log(object$fold_prior[x,]);counts<-integer(K)
  for(r in ix){if(d$.rc_start[r]>=landmark)break;q<-as.integer(sum(pmin(counts,cap)*powers))+1L;end<-min(d$.rc_stop[r],landmark)
   cuts<-sort(unique(c(d$.rc_start[r],end,seq(0,object$horizon,length.out=L+1L)[seq(0,object$horizon,length.out=L+1L)>d$.rc_start[r]&seq(0,object$horizon,length.out=L+1L)<end])))
   for(f in seq_len(nf)){g<-(f-1L)*nx+x;for(b in seq_len(length(cuts)-1L)){bin<-min(L,floor(cuts[b]/h)+1L);logweight[f]<-logweight[f]-(cuts[b+1]-cuts[b])*sum(object$rates[,q,bin,g])}
    if(d$.rc_stop[r]<=landmark&&d$.rc_mark[r]>0){bin<-min(L,floor(d$.rc_stop[r]/h)+1L);logweight[f]<-logweight[f]+log(object$rates[d$.rc_mark[r],q,bin,g])}}
   if(d$.rc_stop[r]<=landmark&&d$.rc_mark[r]>0)counts[d$.rc_mark[r]]<-counts[d$.rc_mark[r]]+1L
  }
  weights<-exp(logweight-max(logweight));weights<-weights/sum(weights);q<-as.integer(sum(pmin(counts,cap)*powers))+1L
  ident<-as.character(d[[sc$unit]][ix[1]]);weightsout[[ident]]<-weights
  for(delta in horizon){tot<-numeric(K+2L);cuts<-sort(unique(c(landmark,landmark+delta,seq(0,object$horizon,length.out=L+1L)[seq(0,object$horizon,length.out=L+1L)>landmark&seq(0,object$horizon,length.out=L+1L)<landmark+delta])))
   for(f in seq_len(nf)){g<-(f-1L)*nx+x;p<-numeric(S+1L);p[q]<-1;rewardtotal<-numeric(K+1L)
    for(b in seq_len(length(cuts)-1L)){bin<-min(L,floor(cuts[b]/h)+1L);Q<-matrix(0,S+1L,S+1L);reward<-matrix(0,S+1L,K+1L)
     for(s in seq_len(S)){for(k in seq_len(M)){dest<-object$next_state[s,k]+1L;if(dest!=s)Q[s,dest]<-Q[s,dest]+object$rates[k,s,bin,g];if(k<=K)reward[s,k]<-object$rates[k,s,bin,g]};Q[s,s]<- -sum(Q[s,-s]);reward[s,K+1L]<-1}
     z<-rc_propagate(p,Q,cuts[b+1]-cuts[b],reward);p<-z$probability;rewardtotal<-rewardtotal+z$reward
    };tot<-tot+weights[f]*c(rewardtotal[seq_len(K)],p[S+1L],rewardtotal[K+1L])
   }
   values<-stats::setNames(tot,object$target_names);values<-c(values,survival=1-values["risk"]);names(values)[length(values)]<-"survival"
   wanted<-if(type=="all")object$target_names else if(type=="burden")grep("^burden_",names(values),value=TRUE) else type
   rowsout[[length(rowsout)+1L]]<-data.frame(unit=ident,landmark=landmark,horizon=delta,target=wanted,estimate=unname(values[wanted]),std_error=NA_real_,lower=NA_real_,upper=NA_real_,interval="none",inference_status="conditional_target_not_population_EIF",time_unit=object$time_unit)
  }
 }
 if(!length(rowsout)) .rc_error("target_eligibility","No observed terminal-event-free subject is eligible at this landmark.")
 ans<-do.call(rbind,rowsout);attr(ans,"fold_posterior")<-weightsout;ans
}

#' @export
predict.coherent4cast_targeted <- function(object,newdata=NULL,landmark=0,horizon=NULL,type="all",interval="auto",level=.95,draws=1000L,seed=1L,...) {
 if(length(list(...))) .rc_error("target_argument","An unsupported targeted prediction argument was supplied.")
 .rc_targeted_predict(object,newdata,landmark,horizon,type,interval,level,draws,seed)
}
#' @export
summary.coherent4cast_targeted <- function(object,...) {
 list(targeting=object$targeting,subjects=object$n,horizon=object$horizon,estimates=data.frame(target=object$target_names,estimate=unname(object$estimate),std_error=unname(object$se),EIF_mean=object$EIF_mean,EIF_variance=object$EIF_variance),convergence=object$convergence,assumptions=object$assumptions)
}
#' @export
print.coherent4cast_targeted <- function(x,...) {
 cat("coherent4cast",x$targeting,"targeted common law\nSubjects:",x$n,"  Horizon:",x$horizon,x$time_unit,"\n")
 print(data.frame(target=x$target_names,estimate=unname(x$estimate),std_error=unname(x$se)),row.names=FALSE)
 cat("Target score:",x$convergence$score,"  Score and variance checks:",x$convergence$valid,"\n");invisible(x)
}

#' @export
tidy.coherent4cast_targeted <- function(x,interval="auto",level=.95,draws=1000L,seed=1L,...) {
 p<-.rc_targeted_predict(x,interval=interval,level=level,draws=draws,seed=seed)
 data.frame(term=p$target,estimate=p$estimate,std.error=p$std_error,conf.low=p$lower,conf.high=p$upper,inference_status=p$inference_status)
}
#' @export
glance.coherent4cast_targeted <- function(x,...) {
 data.frame(n=x$n,targeting=x$targeting,horizon=x$horizon,targets=length(x$estimate),converged=isTRUE(x$convergence$valid),EIF_mean_max=max(abs(x$EIF_mean)),interpolation_error=if(is.null(x$interpolation_error))NA_real_ else x$interpolation_error)
}
#' @export
coef.coherent4cast_targeted <- function(object,...) object$estimate

.rc_target_projection <- function(object,level) {
  limits<-object$target_limits;finite<-which(is.finite(limits[1,])&is.finite(limits[2,]))
  J<-ncol(limits);bounds<-object$observation_bounds
  if(is.null(bounds)) .rc_error("target_projection","No explicit observed-outcome bounds are available for this target.")
  n<-nrow(bounds$lower);lower<-limits[1,];upper<-limits[2,]
  status<-rep("unbounded_burden_outer_support",J)
  if(length(finite))for(j in finite){
    a<-limits[1,j];b<-limits[2,j];lo<-bounds$lower[,j];hi<-bounds$upper[,j]
    if(any(!is.finite(c(lo,hi)))||any(lo<a|hi>b|lo>hi)) .rc_error("target_projection","Observed coarsening bounds contradict the declared reward support.")
    radius<-(b-a)*sqrt(log(2*length(finite)/(1-level))/(2*n))
    lower[j]<-max(a,mean(lo)-radius);upper[j]<-min(b,mean(hi)+radius)
    status[j]<-"finite_sample_coarsening_outer_set"
  }
  list(lower=lower,upper=upper,critical=NA_real_,status=status,method="bounded_coarsening_projection")
}

#' @export
logLik.coherent4cast_targeted <- function(object,...) {
 .rc_error("targeted_likelihood","A targeted population functional has no supplied maximized likelihood for AIC; inspect its EIF and targeting score diagnostics.")
}
