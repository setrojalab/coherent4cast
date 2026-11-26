.cs_split <- function(d, breaks) {
  out <- vector("list", nrow(d)); schema <- attr(d, "cs_schema")
  for (i in seq_len(nrow(d))) {
    cuts <- c(d$.cs_start[i], breaks[breaks > d$.cs_start[i] & breaks < d$.cs_stop[i]], d$.cs_stop[i])
    z <- d[rep(i, length(cuts)-1L), , drop=FALSE]
    z$.cs_start <- head(cuts,-1L);z$.cs_stop <- tail(cuts,-1L)
    if (nrow(z)>1L) {z$.cs_mark[-nrow(z)] <- 0L;z$.cs_terminal[-nrow(z)] <- 0L}
    out[[i]] <- z
  }
  ans <- do.call(rbind,out);rownames(ans)<-NULL;attr(ans,"cs_schema")<-schema
  ans
}

.cs_design <- function(d, formula, breaks, cap, nonrecurrent, terms_object=NULL, xlevels=NULL, contrasts=NULL) {
  d <- .cs_split(d, breaks)
  mf <- stats::model.frame(if(is.null(terms_object)) formula else terms_object, d, na.action=stats::na.fail, xlev=xlevels)
  tt <- attr(mf,"terms");mm<-stats::model.matrix(tt,mf,contrasts.arg=contrasts)
  keep<-colnames(mm)!="(Intercept)";covariates<-mm[,keep,drop=FALSE]
  history<-cs_build_history(d,cap=cap)$counts;K<-ncol(history);J<-length(breaks)-1L
  base<-matrix(0,nrow(d),J);base[cbind(seq_len(nrow(d)),findInterval(d$.cs_start,breaks,all.inside=TRUE))]<-1
  colnames(base)<-paste0("baseline_",seq_len(J));block<-cbind(base,covariates,history)
  p<-ncol(block);M<-K+1L;m<-nrow(d)
  X<-matrix(0,m*M,p*M);event<-numeric(m*M);exposure<-rep(d$.cs_stop-d$.cs_start,M)
  for(k in seq_len(M)) {
    rows<-(k-1L)*m+seq_len(m);cols<-(k-1L)*p+seq_len(p);X[rows,cols]<-block
    event[rows]<-if(k<=K) as.integer(d$.cs_mark==k) else d$.cs_terminal
    if(k<=K && k %in% nonrecurrent) {
      disabled<-history[,k]>0
      if(any(event[rows][disabled]>0)) .cs_error("eligibility", "A declared non-repeatable event recurs.")
      exposure[rows[disabled]]<-0
    }
  }
  names<-as.vector(vapply(seq_len(M),function(k) paste0(if(k<=K) paste0("event",k) else "terminal", ":",colnames(block)),character(p)))
  colnames(X)<-names
  list(unit=rep(d$.cs_unit,M),exposure=exposure,design=X,event=event,data=d,
       terms=tt,xlevels=stats::.getXlevels(tt,mf),contrasts=attr(mm,"contrasts"),
       covariate_columns=colnames(covariates),block_size=p,base_size=J,K=K,names=names)
}

.cs_inverse <- function(H, tolerance=1e-10) {
  e<-eigen((H+t(H))/2,symmetric=TRUE);good<-e$values>max(abs(e$values))*tolerance
  if(!any(good)) return(matrix(0,nrow(H),ncol(H)))
  tcrossprod(sweep(e$vectors[,good,drop=FALSE],2,sqrt(e$values[good]),"/"))
}

#' Fit a common shared-gamma recurrent and terminal event law
#' @export
cs_fit <- function(formula, data, unit="id", start="start", stop="stop", event_type="event_type",
    terminal="terminal", observed=NULL, event_types=NULL, time_unit, breaks=c(0,Inf),
    observation=c("exact","coarsened"), targeting="none", control=cs_control(),
    weights=NULL, subset=NULL, na.action=stats::na.fail, labels=NULL, nonrecurrent=character(),
    coarsening=NULL, intervention=NULL, estimand=NULL, seed=1L,
    frailty=c("gamma","lognormal","discrete","none"),frailty_points=2L,frailty_diagnostics=TRUE,
    targeting_spec=NULL,coarsening_control=list()) {
  frailty_supplied<-!missing(frailty)
  call<-match.call();observation<-match.arg(observation);frailty<-match.arg(frailty)
  frailty_requested<-frailty
  if(frailty=="none") {
    if(!is.null(control$phi_fixed)&&control$phi_fixed!=0) .cs_error("frailty_none","frailty='none' conflicts with a positive fixed dispersion.")
    control$phi_fixed<-0;frailty<-"gamma"
  }
  targeting<-match.arg(targeting,c("none","finite_path","continuous"))
  if(!is.null(intervention)) .cs_error("causal_scope","Interventions are not available before version 2.0.0.")
  if(targeting!="none") {
    if(!is.null(weights)||!is.null(coarsening)||observation!="exact"||length(nonrecurrent)||
       (!missing(subset)&&!identical(substitute(subset),quote(NULL)))||frailty_supplied||!identical(breaks,c(0,Inf)))
      .cs_error("targeting_scope","Targeting uses its declared finite-path MAR or continuous-time CAR model. Specify support, folds, strata and time grid through targeting_spec; likelihood-only weights, frailty, coarsening, eligibility and baseline-break arguments cannot be silently reused.")
    targeted_fit<-.cs_targeted_fit(formula,data,unit,start,stop,event_type,terminal,observed,event_types,time_unit,targeting,control,targeting_spec,seed,call)
    targeted_fit$audit<-c(targeted_fit$audit,.cs_capture_audit_context(data,call))
    targeted_fit$audit$estimand <- estimand
    return(targeted_fit)
  }
  if(!inherits(formula,"formula") || length(formula)!=2L) .cs_error("formula","Supply a one-sided covariate formula such as ~ treatment + age.")
  if(any(all.vars(formula) %in% c(stop,event_type,terminal))) .cs_error("formula_leakage","Interval outcomes and future stop times cannot enter the covariate formula.")
  if(length(breaks)<2L || breaks[1]!=0 || tail(breaks,1)!=Inf || any(diff(breaks)<=0) || anyNA(breaks)) .cs_error("breaks","Baseline breaks must increase from 0 to Inf.")
  if(!missing(subset) && !identical(substitute(subset),quote(NULL))) {
    selected<-eval(substitute(subset),data,parent.frame())
    if(length(selected)!=nrow(data)||anyNA(selected)) .cs_error("subset","subset must select complete units.")
    ids<-data[[unit]];selected<-as.logical(selected)
    if(any(vapply(split(selected,ids),function(x)length(unique(x))>1L,logical(1)))) .cs_error("subset","Row subsets cannot break histories; select whole units.")
    data<-data[selected,,drop=FALSE]
  }
  fibre_audit<-NULL
  if(observation=="coarsened"&&is.null(coarsening)) {
    permitted<-c("lower","upper","followup","format","rounding_unit","nodes","max_perms","draws","seed","max_paths")
    if(!is.list(coarsening_control)||(length(coarsening_control)&&(is.null(names(coarsening_control))||any(!names(coarsening_control)%in%permitted)))) .cs_error("coarsening_control","Use named event-window and numerical controls; fitted schema and data cannot be replaced through coarsening_control.")
    if(!is.null(observed)) {
      if(!observed %in% names(data)||anyNA(data[[observed]])||any(!data[[observed]] %in% c(1,TRUE))) .cs_error("unrecorded_history","Automatic time fibres require known event counts; an unmonitored interval is not a zero-event interval.")
      observed<-NULL
    }
    fibre_args<-utils::modifyList(list(data=data,unit=unit,start=start,stop=stop,event_type=event_type,terminal=terminal,event_types=event_types,time_unit=time_unit,nonrecurrent=nonrecurrent,breaks=breaks,covariates=all.vars(formula),seed=seed),coarsening_control)
    fibre<-do.call(cs_coarsen_fibre,fibre_args);data<-fibre$data;coarsening<-fibre$coarsening
    fibre_audit<-list(metadata=fibre$metadata,tie_audit=fibre$tie_audit)
  }
  needed<-unique(c(all.vars(formula),unit,start,stop,event_type,terminal,observed))
  if(!all(needed %in% names(data))) .cs_error("columns","A formula or trajectory column is absent.")
  missing_rows<-!stats::complete.cases(data[,needed,drop=FALSE]);omitted<-character()
  if(any(missing_rows)) {
    if(identical(na.action,stats::na.omit) || identical(na.action,"na.omit")) {
      omitted<-unique(as.character(data[[unit]][missing_rows]));data<-data[!as.character(data[[unit]]) %in% omitted,,drop=FALSE]
    } else .cs_error("missing","Missing values detected. na.omit removes complete units, not isolated history rows.")
  }
  d<-cs_as_trajectory(data,unit,start,stop,event_type,terminal,observed,event_types,time_unit,observation)
  schema<-attr(d,"cs_schema")
  if(length(schema$gaps)||any(!d$.cs_observed)) .cs_error("unrecorded_history","Exact fitting cannot interpret unmonitored intervals as zero events. Supply compatible paths or use finite-support bounds.")
  if(length(schema$delayed_entry)) .cs_error("left_truncation","This fit does not model frailty selection before delayed entry; start the observed history at origin.")
  if(!all(nonrecurrent %in% schema$event_types)) .cs_error("nonrecurrent","A non-repeatable label is not in event_types.")
  nr<-match(nonrecurrent,schema$event_types)
  if(length(nr)&&control$history_cap<1L) .cs_error("history_cap","Non-repeatable event eligibility needs history_cap >= 1.")
  dd<-.cs_design(d,formula,breaks,control$history_cap,nr)
  n<-length(unique(d$.cs_unit));full_size<-ncol(dd$design)
  active<-which(colSums(abs(dd$design)*(dd$exposure>0|dd$event>0))>0)
  dd$active<-active;dd$full_size<-full_size;dd$full_names<-dd$names
  dd$design<-dd$design[,active,drop=FALSE];dd$names<-dd$names[active];p<-length(active)
  if(!sum(dd$event)) .cs_error("no_events","No event of any declared type was observed; likelihood parameters are not identifiable.")
  sw<-if(is.null(weights)) rep(1,n) else weights
  if((observation=="exact" && length(sw)!=n) || any(!is.finite(sw)) || any(sw<=0)) .cs_error("weights","weights must contain one positive finite case weight per unit, in first-appearance order.")
  owner<-NULL;lw<-NULL
  if(observation=="coarsened") {
    if(is.null(coarsening)||is.null(coarsening$owner)||is.null(coarsening$log_weights)) .cs_error("coarsening","Supply path-to-subject owner and log_weights for legal quadrature paths.")
    owner<-coarsening$owner;lw<-coarsening$log_weights
    if(length(owner)!=n||length(lw)!=n) .cs_error("coarsening","One owner and log weight are required per complete path.")
    nobs<-length(unique(owner));if(!is.null(weights)&&length(weights)!=nobs) .cs_error("weights","Coarsened weights must be one per observed subject.")
    sw<-if(is.null(weights)) rep(1,nobs) else weights
  }
  start_coef<-numeric(full_size)
  for(k in seq_len(dd$K+1L)) {
    rows<-(k-1L)*nrow(dd$data)+seq_len(nrow(dd$data));rate<-(sum(dd$event[rows])+.1)/(sum(dd$exposure[rows])+.1)
    cols<-(k-1L)*dd$block_size+seq_len(dd$base_size);start_coef[cols]<-log(rate)
  }
  start_coef<-start_coef[active]
  frailty_spec<-.cs_frailty_spec(frailty,control,frailty_points)
  theta0<-c(start_coef,frailty_spec$start)
  lower<-c(rep(control$coefficient_bounds[1],p),frailty_spec$lower)
  upper<-c(rep(control$coefficient_bounds[2],p),frailty_spec$upper)
  evaluate<-function(theta,info=FALSE) {
    z<-.cs_fit_likelihood(theta,dd,control,frailty,frailty_points,info)
    if(!is.null(owner)) z<-cs_coarsen(z$subject_loglik,z$scores,if(info)z$information_individual else NULL,owner,lw)
    z$objective <- -sum(sw*z$subject_loglik)/sum(sw)
    z$gradient <- -colSums(z$scores*sw)/sum(sw);z
  }
  fn<-function(th)evaluate(th)$objective;gr<-function(th)evaluate(th)$gradient
  opt<-stats::optim(pmin(pmax(theta0,lower+.01),upper-.01),fn,gr,method="L-BFGS-B",lower=lower,upper=upper,
                    control=list(maxit=control$maxit,factr=max(1,min(100,control$tolerance/.Machine$double.eps)),pgtol=control$tolerance))
  multistart_values<-opt$value
  if(frailty=="discrete")for(scale in c(.5,1.8)) {
    init<-theta0;init[p+seq_along(frailty_spec$start)]<-frailty_spec$start*scale
    alt<-stats::optim(pmin(pmax(init,lower+.01),upper-.01),fn,gr,method="L-BFGS-B",lower=lower,upper=upper,
      control=list(maxit=control$maxit,factr=max(1,min(100,control$tolerance/.Machine$double.eps)),pgtol=control$tolerance))
    multistart_values<-c(multistart_values,alt$value);if(alt$value<opt$value)opt<-alt
  }
  quadrature_audit<-data.frame();quadrature_valid<-TRUE
  if(frailty=="lognormal"&&(is.null(control$phi_fixed)||control$phi_fixed>0)) {
    maximum<-if(is.null(control$quadrature_max))1024L else control$quadrature_max
    qt<-if(is.null(control$quadrature_tolerance))1e-7 else control$quadrature_tolerance
    quadrature_valid<-FALSE
    if(control$quadrature==maximum&&maximum>=4L) {
      high<-evaluate(opt$par);control$quadrature<-as.integer(floor(maximum/2));low<-evaluate(opt$par)
      gap<-max(abs(high$subject_loglik-low$subject_loglik));sgap<-max(abs(high$gradient-low$gradient))
      quadrature_audit<-data.frame(nodes_low=control$quadrature,nodes_high=maximum,max_subject_loglik_change=gap,max_mean_score_change=sgap)
      control$quadrature<-maximum;quadrature_valid<-gap<=qt&&sgap<=qt
    }
    while(control$quadrature<maximum) {
      low<-evaluate(opt$par);nlow<-control$quadrature
      control$quadrature<-min(maximum,2L*nlow);high<-evaluate(opt$par)
      gap<-max(abs(high$subject_loglik-low$subject_loglik));sgap<-max(abs(high$gradient-low$gradient))
      nnew<-control$quadrature
      opt<-stats::optim(opt$par,fn,gr,method="L-BFGS-B",lower=lower,upper=upper,
        control=list(maxit=control$maxit,factr=max(1,min(100,control$tolerance/.Machine$double.eps)),pgtol=control$tolerance))
      # The optimizer can move after refinement. Check the comparison at the
      # final optimized parameter, not only at the warm-start parameter.
      high<-evaluate(opt$par);control$quadrature<-nlow;low<-evaluate(opt$par);control$quadrature<-nnew
      gap<-max(abs(high$subject_loglik-low$subject_loglik));sgap<-max(abs(high$gradient-low$gradient))
      quadrature_audit<-rbind(quadrature_audit,data.frame(nodes_low=nlow,nodes_high=nnew,max_subject_loglik_change=gap,max_mean_score_change=sgap))
      if(gap<=qt&&sgap<=qt){quadrature_valid<-TRUE;break}
    }
    frailty_spec$points<-control$quadrature
  }
  z<-evaluate(opt$par,TRUE)
  H<-if(length(dim(z$information_individual))==3L) apply(sweep(z$information_individual,3,sw,"*"),c(1,2),sum) else z$information
  inv<-.cs_inverse(H);scores<-z$scores*sw;V<-inv%*%crossprod(scores)%*%inv
  eig<-eigen((H+t(H))/2,symmetric=TRUE,only.values=TRUE)$values
  at_bound<-any(abs(opt$par-lower)<1e-5|abs(opt$par-upper)<1e-5)
  valid<-opt$convergence==0L && min(eig)>1e-8 && !at_bound && max(abs(z$gradient))<1e-4 && quadrature_valid
  nm<-c(dd$names,frailty_spec$names);names(opt$par)<-nm;dimnames(V)<-list(nm,nm)
  hashfile<-tempfile();on.exit(unlink(hashfile),add=TRUE);saveRDS(data,hashfile,version=2)
  linear_coefficients<-stats::setNames(numeric(full_size),dd$full_names);linear_coefficients[active]<-opt$par[seq_len(p)]
  fit<-list(coefficients=opt$par,linear_coefficients=linear_coefficients,covariance=V,information=H,scores=z$scores,
      phi=if(frailty=="gamma") {if(is.null(control$phi_fixed)) exp(tail(opt$par,1)) else control$phi_fixed} else .cs_frailty_nodes(opt$par[-seq_len(p)],frailty_spec)$variance,
      frailty=frailty,frailty_spec=frailty_spec,targeting="none",case_weights=sw,
      loglik=sum(sw*z$subject_loglik),n=length(sw),npaths=n,data=d,design=dd,formula=formula,
      breaks=breaks,control=control,schema=schema,nonrecurrent=nr,labels=labels,coarsening=coarsening,
      convergence=list(code=if(!quadrature_valid)"quadrature_unresolved" else if(valid)"converged" else if(at_bound)"parameter_boundary" else "nonregular_fit",
        optimizer=opt$convergence,message=opt$message,gradient=max(abs(z$gradient)),information_eigenvalues=eig,at_boundary=at_bound,valid=valid),
      audit=list(version=as.character(utils::packageVersion("coherentSurv")),seed=seed,call=call,data_hash=unname(tools::md5sum(hashfile)),
       data_hash_algorithm="MD5 integrity fingerprint, not a security signature",omitted_units=omitted,structurally_inactive_coefficients=dd$full_names[-active],estimand=estimand,
       automatic_fibre=fibre_audit,multistart_objectives=multistart_values,quadrature_refinement=quadrature_audit,quadrature_stable=quadrature_valid,
       assumptions=c("independent units","declared observation model","sufficient retained state","predictable covariates","one absorbing terminal event")))
  class(fit)<-"coherentSurv"
  fit$frailty_requested<-frailty_requested
  fit$audit<-c(fit$audit,.cs_capture_audit_context(d,call))
  if(isTRUE(frailty_diagnostics))fit$frailty_diagnostics<-cs_frailty_diagnose(fit)
  if(!valid) .cs_warning("nonregular_fit", "Fit is nonregular or did not meet convergence criteria; inspect $convergence. Wald inference is not certified.")
  fit
}
