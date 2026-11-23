#' Generate the canonical colon observation fibre automatically
#' @export
rc_colon_fibre <- function(data=NULL, nodes=4L, covariate=c("node4","nodes"),
    node_median=NULL, breaks=c(0,1,3,Inf), max_perms=120L,draws=2000L,
    seed=20261122L,max_paths=100000L) {
  covariate<-match.arg(covariate)
  raw<-.rc_source_data("colon","survival",data)
  .rc_adapter_columns(raw,c("id","etype","time","status","rx","sex","age","node4","nodes"))
  if(anyNA(raw$id)||anyNA(raw$etype)||any(!raw$etype%in%c(1,2))) .rc_error("adapter_id","Colon requires paired recurrence/death rows per patient.")
  if(covariate=="nodes") {
    if(is.null(node_median))node_median<-stats::median(raw$nodes,na.rm=TRUE)
    if(length(node_median)!=1L||!is.finite(node_median)||node_median<0).rc_error("node_median","Supply a finite nonnegative training median.")
  }
  grouped<-split(raw,raw$id);records<-vector("list",length(grouped));truncated<-character()
  for(i in seq_along(grouped)) {
    g<-grouped[[i]]
    if(nrow(g)!=2L||!setequal(g$etype,c(1,2))).rc_error("adapter_pair","Exactly one recurrence and one death row is required per patient.")
    r<-g[g$etype==1,,drop=FALSE];d<-g[g$etype==2,,drop=FALSE]
    if(anyNA(g[,c("time","status","rx","sex","age","node4")])||any(!is.finite(g$time))||any(g$time<=0)||any(!g$status%in%c(0,1)))
      .rc_error("adapter_time","Colon requires nonmissing positive times and binary status.")
    if(any(vapply(c("rx","sex","age","node4","nodes"),function(v)!identical(r[[v]],d[[v]]),logical(1))))
      .rc_error("adapter_covariates","Paired endpoint rows disagree on baseline covariates.")
    end<-d$time;dead<-d$status
    if(r$status==0&&r$time<end){end<-r$time;dead<-0;truncated<-c(truncated,as.character(r$id))}
    eventday<-c(r$time,end);observed<-c(r$status,dead)
    z<-data.frame(id=r$id,event_time_lower=ifelse(observed==1,(eventday-1)/365.25,NA_real_),
      event_time_upper=ifelse(observed==1,eventday/365.25,NA_real_),
      event_type=c(if(r$status==1)"1" else "0","0"),terminal=c(0,dead),followup=end/365.25,
      levamisole=as.integer(as.character(r$rx)=="Lev"),levamisole_5fu=as.integer(as.character(r$rx)=="Lev+5FU"),
      male=r$sex,age_per_10=(r$age-60)/10,stringsAsFactors=FALSE)
    if(covariate=="node4")z$node4<-r$node4 else {
      z$nodes_per_5<-if(is.na(r$nodes))node_median/5 else r$nodes/5;z$nodes_missing<-as.integer(is.na(r$nodes))
    }
    records[[i]]<-z
  }
  records<-do.call(rbind,records);records<-records[order(records$id),];rownames(records)<-NULL
  out<-rc_coarsen_fibre(records,event_types="1",nonrecurrent="1",time_unit="years",breaks=breaks,
    nodes=nodes,max_perms=max_perms,draws=draws,seed=seed,max_paths=max_paths,format="windows")
  out$formula<-if(covariate=="node4")~levamisole+levamisole_5fu+male+age_per_10+node4 else
    ~levamisole+levamisole_5fu+male+age_per_10+nodes_per_5+nodes_missing
  out$observed_data<-records
  out$metadata$dataset<-"survival::colon"
  out$metadata$source_version<-if(requireNamespace("survival",quietly=TRUE))as.character(utils::packageVersion("survival")) else NA_character_
  out$metadata$patients<-length(grouped);out$metadata$recurrences<-sum(raw$status[raw$etype==1]);out$metadata$source_deaths<-sum(raw$status[raw$etype==2])
  out$metadata$covariate<-covariate;out$metadata$node_median<-node_median
  out$metadata$surveillance_truncated_ids<-truncated
  out$metadata$same_day_ids<-as.numeric(out$tie_audit$source_id[out$tie_audit$coincident_upper_endpoints>0])
  out
}
