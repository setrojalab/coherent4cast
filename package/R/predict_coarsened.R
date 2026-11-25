#' Forecast directly from coarsened histories observed through a landmark
#' @export
rc_predict_coarsened <- function(object, history, landmark, horizon,
    type=c("risk","burden","survival","rmst","all"),
    quadrature=object$control$quadrature,coarsening_control=list(),...) {
  type<-match.arg(type)
  if(!inherits(object,"coherent4cast")) .rc_error("fit","object must be a coherent4cast fit.")
  if(!is.data.frame(history)||!nrow(history)) .rc_error("coarsened_history","Supply the observed history available through the landmark.")
  if(length(landmark)!=1L||!is.numeric(landmark)||!is.finite(landmark)||landmark<=0)
    .rc_error("landmark","Automatic history fibres require a positive finite landmark.")
  if(!is.list(coarsening_control)||(length(coarsening_control)&&is.null(names(coarsening_control)))||any(!nzchar(names(coarsening_control)))) .rc_error("coarsening_control","Use a named coarsening-control list.")
  permitted<-c("lower","upper","followup","format","rounding_unit","nodes","max_perms","draws","seed","max_paths")
  if(any(!names(coarsening_control)%in%permitted)) .rc_error("coarsening_control","History coarsening controls cannot replace fitted column definitions, marks, time units or baseline breaks.")
  extra<-list(...)
  if(length(extra)&&(is.null(names(extra))||any(!names(extra)%in%c("covariate_evolution","newdata_future","covariate_model","covariate_history","covariate_draws","covariate_grid","seed"))))
    .rc_error("unused_argument","An unsupported coarsened prediction argument was supplied.")
  mode<-if(is.null(extra$covariate_evolution))"freeze" else match.arg(extra$covariate_evolution,c("freeze","locf","joint_submodel"))
  if(mode=="joint_submodel"&&is.null(extra$covariate_history))
    .rc_error("covariate_history","The joint covariate submodel requires an explicit measured history; latent event times are not measurement times.")
  s<-object$schema
  args<-utils::modifyList(list(data=history,unit=s$unit,start=s$start,stop=s$stop,event_type=s$event_type,
    terminal=s$terminal,event_types=s$event_types,time_unit=s$time_unit,nonrecurrent=s$event_types[object$nonrecurrent],
    breaks=object$breaks,covariates=all.vars(object$formula),
    lower="event_time_lower",upper="event_time_upper",followup="followup",format="auto",
    seed=if(is.null(extra$seed))1L else extra$seed),coarsening_control)
  if(!all(c(s$unit,s$terminal)%in%names(history))) .rc_error("coarsened_history","The history lacks fitted subject or terminal columns.")
  if(anyNA(history[[s$terminal]])||any(history[[s$terminal]]!=0)) .rc_error("landmark_terminal","Every prediction subject must be terminal-event-free through the landmark.")
  if(!is.null(s$observed)&&s$observed%in%names(history)&&any(is.na(history[[s$observed]])|!history[[s$observed]]%in%c(1,TRUE)))
    .rc_error("coarsened_history","Unrecorded intervals with unknown event counts require another observation model.")
  format<-args$format
  if(format=="auto")format<-if(all(c(args$lower,args$upper)%in%names(history)))"windows" else "counting"
  if(format=="windows") {
    if(!all(c(args$lower,args$upper,args$followup)%in%names(history))) .rc_error("coarsened_history","Supply all window and follow-up columns.")
    ends<-history[[args$followup]]
    windows<-c(history[[args$lower]],history[[args$upper]])
    if(!is.numeric(ends)||any(!is.finite(ends))||any(abs(ends-landmark)>1e-12*max(1,landmark)))
      .rc_error("path_future","Every observed follow-up endpoint must equal the prediction landmark.")
    if(!is.numeric(windows)||any(windows[!is.na(windows)]>landmark))
      .rc_error("path_future","A future or landmark-straddling event window cannot enter the past-history fibre.")
  } else if(format=="counting") {
    if(!all(c(s$start,s$stop)%in%names(history))) .rc_error("coarsened_history","Supply start/stop history columns.")
    if(!is.numeric(history[[s$stop]])||any(!is.finite(history[[s$stop]]))||any(history[[s$stop]]>landmark))
      .rc_error("path_future","Counting rows cannot contain follow-up after the landmark.")
    ends<-vapply(split(history[[s$stop]],history[[s$unit]]),max,numeric(1))
    if(any(abs(ends-landmark)>1e-12*max(1,landmark))) .rc_error("path_history","Every history must end exactly at the landmark.")
  } else .rc_error("coarsening_control","format must be auto, windows or counting.")
  fibre<-do.call(rc_coarsen_fibre,args)
  path_map<-unique(fibre$data[,c(s$unit,".rc_source_unit"),drop=FALSE])
  replicate_schedule<-function(table,measured=FALSE) {
    if(is.null(table))return(NULL)
    if(!is.data.frame(table)||!all(c(s$unit,"time")%in%names(table))||anyNA(table[[s$unit]])||
       !is.numeric(table$time)||any(!is.finite(table$time)))
      .rc_error("covariate_history","Covariate tables require the fitted unit column and finite numeric time.")
    if(any(!as.character(table[[s$unit]])%in%fibre$metadata$source_ids))
      .rc_error("covariate_history","A covariate table contains subjects absent from the prediction histories.")
    if(measured) {
      if(any(table$time<0)) .rc_error("covariate_history","Measurement times must be nonnegative.")
      table<-table[table$time<=landmark,,drop=FALSE]
      if(!all(fibre$metadata$source_ids%in%as.character(table[[s$unit]])))
        .rc_error("covariate_history","Each forecast subject needs at least one actual measurement through the landmark.")
    }
    if(anyDuplicated(table[,c(s$unit,"time"),drop=FALSE]))
      .rc_error("covariate_history","Provide one measured or scheduled value per source subject and time; event rows are not repeated measurements.")
    copies<-lapply(seq_len(nrow(path_map)),function(i) {
      z<-table[as.character(table[[s$unit]])==path_map$.rc_source_unit[i],,drop=FALSE]
      z[[s$unit]]<-rep(path_map[[s$unit]][i],nrow(z));z
    })
    ans<-do.call(rbind,copies);rownames(ans)<-NULL;ans
  }
  if(!is.null(extra$newdata_future))extra$newdata_future<-replicate_schedule(extra$newdata_future)
  if(!is.null(extra$covariate_history))extra$covariate_history<-replicate_schedule(extra$covariate_history,TRUE)
  # A candidate complete path has no observation mask. Its likelihood is
  # subsequently averaged over the declared observation fibre.
  complete_object<-object;complete_object$schema$observed<-NULL
  out<-do.call(rc_predict_paths,c(list(object=complete_object,paths=fibre$data,owner=fibre$coarsening$owner,
    log_weights=fibre$coarsening$log_weights,landmark=landmark,horizon=horizon,type=type,quadrature=quadrature),extra))
  owners<-suppressWarnings(as.integer(out$unit))
  if(anyNA(owners)||any(owners<1|owners>length(fibre$metadata$source_ids)))
    .rc_error("path_owner","A mixture forecast lost its source subject mapping.")
  out$unit<-fibre$metadata$source_ids[owners]
  weights<-attr(out,"path_weights")
  if(!is.null(weights)){weights$source_unit<-fibre$metadata$source_ids[weights$owner];attr(out,"path_weights")<-weights}
  attr(out,"automatic_fibre")<-list(metadata=fibre$metadata,tie_audit=fibre$tie_audit)
  attr(out,"assumptions")<-c(attr(out,"assumptions"),
    "automatic fibre used only records ending at the declared landmark",
    "the count and marks of recorded events are complete; event timing/order may be coarsened",
    "original observed subject identifiers retained in unit",
    "future schedules and actual measurement histories mapped from observed subjects to candidate paths; latent event times are not measurements")
  out
}
