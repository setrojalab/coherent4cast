#' Automatically construct legal coarsening paths
#' @export
rc_coarsen_fibre <- function(data, unit="id", lower="event_time_lower",
    upper="event_time_upper", start="start", stop="stop", event_type="event_type",
    terminal="terminal", followup="followup", event_types=NULL, time_unit,
    rounding_unit=NULL, nonrecurrent=character(), breaks=c(0,Inf), nodes=4L,
    max_perms=120L, draws=2000L, seed=20261122L, max_paths=100000L,
    format=c("auto","windows","counting"), covariates=NULL) {
  format<-match.arg(format)
  scalar<-function(x,minimum,integer=FALSE) is.numeric(x)&&length(x)==1L&&is.finite(x)&&x>=minimum&&(!integer||x==floor(x))
  if(missing(time_unit)||length(time_unit)!=1L||is.na(time_unit)||!nzchar(time_unit))
    .rc_error("time_unit","Declare one nonempty time unit; rounding units are never inferred.")
  if(!is.data.frame(data)||!nrow(data)||anyDuplicated(names(data))) .rc_error("fibre_data","Supply a nonempty data frame with unique column names.")
  if(!all(c(unit,event_type,terminal)%in%names(data))) .rc_error("fibre_columns","Unit, event type and terminal columns are required.")
  if(!is.null(rounding_unit)&&!scalar(rounding_unit,.Machine$double.xmin)) .rc_error("rounding_unit","rounding_unit must be positive in the declared model time unit.")
  if(!scalar(max_perms,1,TRUE)||!scalar(draws,2,TRUE)||!scalar(seed,0,TRUE)||seed>2^53-1||!scalar(max_paths,1,TRUE)||
    max(c(max_perms,draws,max_paths))>.Machine$integer.max)
    .rc_error("fibre_control","Invalid permutation, sampling, seed or path budget.")
  if(!is.numeric(breaks)||length(breaks)<2L||anyNA(breaks)||breaks[1]!=0||tail(breaks,1)!=Inf||any(diff(breaks)<=0))
    .rc_error("breaks","Baseline breaks must increase from zero to infinity.")
  quad<-.rc_legendre_unit(nodes)
  if(format=="auto") format<-if(all(c(lower,upper)%in%names(data)))"windows" else "counting"
  needed<-if(format=="windows")c(lower,upper,followup) else c(start,stop)
  if(!all(needed%in%names(data))) .rc_error("fibre_columns",paste("Missing",format,"observation columns."))
  if(anyNA(data[[unit]])||any(!nzchar(as.character(data[[unit]])))) .rc_error("fibre_unit","Observed subject identifiers cannot be missing or empty.")
  rawmark<-as.character(data[[event_type]])
  if(anyNA(rawmark)||anyNA(data[[terminal]])||any(!data[[terminal]]%in%c(0,1))) .rc_error("fibre_event","Event labels must be nonmissing and terminal must be binary.")
  none<-rawmark%in%c("0","none","")
  if(is.null(event_types))event_types<-sort(unique(rawmark[!none]))
  event_types<-as.character(event_types)
  if(anyNA(event_types)||anyDuplicated(event_types)||any(event_types%in%c("0","none",""))||any(!none&!rawmark%in%event_types))
    .rc_error("event_types","Declare distinct recurrent marks covering every recorded event.")
  if(any(!nonrecurrent%in%event_types)) .rc_error("nonrecurrent","Every nonrepeatable mark must appear in event_types.")
  excluded<-unique(c(unit,lower,upper,start,stop,event_type,terminal,followup))
  if(is.null(covariates))covariates<-setdiff(names(data),excluded)
  if(anyDuplicated(covariates)||any(!covariates%in%names(data))||any(covariates%in%excluded)) .rc_error("fibre_covariates","Declare existing covariate columns distinct from observation fields.")
  if(any(covariates%in%c(".rc_source_unit",".rc_path_draw"))) .rc_error("fibre_covariates","Covariate names conflict with reserved provenance columns.")
  ids<-unique(as.character(data[[unit]]));out<-list();owners<-integer();logweights<-numeric();drawids<-integer();audit<-list()
  path_id<-0L
  for(u in seq_along(ids)) {
    rows<-which(as.character(data[[unit]])==ids[u]);g<-data[rows,,drop=FALSE]
    if(length(covariates)&&any(vapply(covariates,function(v)length(unique(g[[v]]))!=1L||anyNA(g[[v]]),logical(1))))
      .rc_error("fibre_covariates","Automatic fibre generation requires subject-constant supplied covariates. Supply a separate timed covariate schedule for varying covariates.")
    if(format=="counting") {
      if(!is.numeric(g[[start]])||!is.numeric(g[[stop]])||any(!is.finite(c(g[[start]],g[[stop]])))||any(g[[start]]<0|g[[stop]]<g[[start]]))
        .rc_error("fibre_time","Counting rows require finite 0 <= start <= stop; ties need an explicit rounding unit.")
      ord<-order(g[[start]],g[[stop]]);gs<-g[ord,,drop=FALSE]
      if(gs[[start]][1L]>0)
        .rc_error("fibre_left_truncation","Automatic timing fibres require complete follow-up from time zero; delayed entry cannot be replaced by invented pre-entry exposure.")
      if(nrow(gs)>1L&&any(gs[[start]][-1L]>cummax(gs[[stop]])[-nrow(gs)]+1e-12))
        .rc_error("fibre_gap","An unrecorded gap does not identify the number of hidden events. Supply explicit observed-event windows or another observation model.")
      end<-max(g[[stop]])
      lo<-hi<-g[[stop]]
      if(!is.null(rounding_unit))lo<-pmax(0,hi-rounding_unit)
    } else {
      if(!is.numeric(g[[followup]])||any(!is.finite(g[[followup]]))||length(unique(g[[followup]]))!=1L)
        .rc_error("fibre_followup","Each observed subject requires one finite common follow-up endpoint.")
      end<-g[[followup]][1L];lo<-g[[lower]];hi<-g[[upper]]
      if(!is.numeric(lo)||!is.numeric(hi)) .rc_error("fibre_time","Window endpoints must be numeric.")
      if(!is.null(rounding_unit)){point<-!is.na(lo)&!is.na(hi)&lo==hi;lo[point]<-pmax(0,hi[point]-rounding_unit)}
    }
    if(!is.finite(end)||end<=0) .rc_error("fibre_followup","Follow-up must be positive and finite.")
    mark<-match(rawmark[rows],event_types);mark[none[rows]]<-0L
    eventlo<-eventhi<-numeric();eventmark<-integer();eventterminal<-logical()
    for(j in seq_len(nrow(g))) {
      if(mark[j]>0){eventlo<-c(eventlo,lo[j]);eventhi<-c(eventhi,hi[j]);eventmark<-c(eventmark,mark[j]);eventterminal<-c(eventterminal,FALSE)}
      if(g[[terminal]][j]==1){eventlo<-c(eventlo,lo[j]);eventhi<-c(eventhi,hi[j]);eventmark<-c(eventmark,length(event_types)+1L);eventterminal<-c(eventterminal,TRUE)}
    }
    if(any(!is.finite(c(eventlo,eventhi)))||any(eventlo<0|eventhi<eventlo|eventhi>end|eventhi<=0))
      .rc_error("fibre_time","Recorded event windows must be finite subsets of (0, followup].")
    if(sum(eventterminal)>1) .rc_error("fibre_terminal","There can be only one observed terminal event per subject.")
    if(any(vapply(match(nonrecurrent,event_types),function(k)sum(eventmark==k)>1L,logical(1))))
      .rc_error("fibre_eligibility","A declared nonrepeatable mark occurs more than once in this observation.")
    generated<-.rc_kernel_guard(cpp_rc_coarsen_fibre(eventlo,eventhi,as.integer(eventmark),eventterminal,end,
      breaks[is.finite(breaks)],quad$x,quad$w,as.integer(max_perms),as.integer(draws),
      as.double((seed+u-1)%%(2^53-1)),as.integer(max_paths)))
    n<-nrow(generated$times);m<-length(eventmark)
    for(p in seq_len(n)) {
      path_id<-path_id+1L
      order<-if(m)generated$order[p,] else integer();times<-if(m)generated$times[p,] else numeric()
      terminal_path<-if(m)eventterminal[order] else logical()
      marked<-character()
      if(m){marked<-rep("0",m);marked[!terminal_path]<-event_types[eventmark[order[!terminal_path]]]}
      if(!m||!tail(terminal_path,1L)) {
        if(!m||tail(times,1L)<end){times<-c(times,end);marked<-c(marked,"0");terminal_path<-c(terminal_path,FALSE)}
      }
      z<-data.frame(path=rep(path_id,length(times)),left=c(0,head(times,-1L)),right=times,
        mark=marked,dead=as.integer(terminal_path),stringsAsFactors=FALSE)
      names(z)<-c(unit,start,stop,event_type,terminal)
      z$.rc_source_unit<-ids[u];z$.rc_path_draw<-generated$draw[p]
      for(v in covariates)z[[v]]<-rep(g[[v]][1L],nrow(z))
      out[[path_id]]<-z
    }
    owners<-c(owners,rep(u,n));logweights<-c(logweights,generated$log_weights);drawids<-c(drawids,generated$draw)
    audit[[u]]<-data.frame(owner=u,source_id=ids[u],recorded_events=m,terminal_events=sum(eventterminal),
      coincident_upper_endpoints=sum(duplicated(eventhi)),method=generated$method,reason=generated$reason,
      canonical_orders_enumerated=generated$permutations,paths=n,attempted_draws=generated$attempted,
      rejected_draws=generated$rejected,integration_volume=generated$volume,volume_mcse=generated$volume_mcse,
      importance_mass_ess=generated$importance_ess,seed=(seed+u-1)%%(2^53-1),stringsAsFactors=FALSE)
  }
  expanded<-do.call(rbind,out);rownames(expanded)<-NULL
  as_trajectory(expanded,unit,start,stop,event_type,terminal,event_types=event_types,time_unit=time_unit,observation="coarsened")
  result<-list(data=expanded,coarsening=list(owner=owners,log_weights=logweights),
    event_types=event_types,nonrecurrent=nonrecurrent,breaks=breaks,tie_audit=do.call(rbind,audit),
    path_draw=drawids,metadata=list(time_unit=time_unit,format=format,rounding_unit=rounding_unit,
      subjects=length(ids),source_ids=ids,paths=path_id,nodes=nodes,max_perms=max_perms,draws=draws,seed=seed,
      covariates=covariates,sampled=any(vapply(audit,function(x)x$method=="importance_sampling",logical(1))),
      observation="Recorded event counts and marks are complete; event times/order may be coarsened. Censoring times are fixed endpoints.",
      duplicate_policy="Canonical earliest-deadline matching avoids duplicate same-mark histories under overlapping windows.",
      validity="Deterministic paths use Gaussian quadrature, not an analytic integral. Sampled paths use inverse-proposal integration masses divided by attempted draws; their likelihood and Louis information are Monte Carlo approximations.",
      forecast_scope="Generate a separate history-only fibre at each landmark; never condition earlier forecasts on full follow-up."))
  class(result)<-c("recurcast_fibre","list");result
}

#' Monte Carlo error of an importance-sampled coarsening likelihood
#' @export
rc_fibre_mcse <- function(fibre, path_loglik) {
  if(!inherits(fibre,"recurcast_fibre")||!is.numeric(path_loglik)||length(path_loglik)!=length(fibre$coarsening$owner)||any(!is.finite(path_loglik)))
    .rc_error("fibre_mcse","Supply a generated fibre and one finite complete-path log likelihood per path.")
  ans<-fibre$tie_audit[,c("owner","source_id","method","attempted_draws"),drop=FALSE]
  ans$log_likelihood<-ans$relative_likelihood_mcse<-ans$log_likelihood_mcse_delta<-NA_real_
  for(i in seq_len(nrow(ans))) {
    rows<-which(fibre$coarsening$owner==ans$owner[i]);z<-path_loglik[rows]+fibre$coarsening$log_weights[rows]
    shift<-max(z);w<-exp(z-shift);total<-sum(w);ans$log_likelihood[i]<-shift+log(total)
    if(ans$method[i]=="importance_sampling") {
      M<-ans$attempted_draws[i]
      se<-sqrt(max(0,M/(M-1)*(sum(w*w)-total^2/M)))
      ans$relative_likelihood_mcse[i]<-se/total
      ans$log_likelihood_mcse_delta[i]<-se/total
    }
  }
  attr(ans,"scope")<-"Conditional integration error at fixed parameters. The log-scale error uses a first-order delta approximation; it is not sampling uncertainty of the fitted estimator. Quadrature error is not estimated here."
  ans
}
