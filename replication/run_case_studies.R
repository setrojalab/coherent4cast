# Execute from the project root with an installed recurcast package.
# Optional environment variable RECURCAST_LIBRARY selects a local R library.
# Public patient records are loaded from their source packages; outputs are
# aggregate results. Local fit objects remain under build/case_fits/.
Sys.unsetenv("LC_ALL")
library_path <- Sys.getenv("RECURCAST_LIBRARY", "")
if (nzchar(library_path)) .libPaths(c(library_path, .libPaths()))
library(recurcast)
output <- Sys.getenv("RECURCAST_CASE_OUTPUT", "replication/results/cases")
cache_directory <- Sys.getenv("RECURCAST_CASE_CACHE", "build/case_fits")
dir.create(output, recursive=TRUE, showWarnings=FALSE)
dir.create(cache_directory, recursive=TRUE, showWarnings=FALSE)
overview <- metrics <- curves <- list()
previous_summary <- if (file.exists(file.path(output,"case_summary.csv"))) read.csv(file.path(output,"case_summary.csv")) else NULL
fit_model <- function(label, adapter, control) {
  cat("Fitting", label, "\n"); flush.console()
  warnings <- character(); began <- proc.time()[["elapsed"]]
  cached <- file.path(cache_directory,paste0(label,".rds"))
  reuse <- identical(Sys.getenv("RECURCAST_REUSE_CASE_FITS"),"1") && file.exists(cached)
  if (reuse) {
    fit <- readRDS(cached)
    checkfile <- tempfile();saveRDS(adapter$data,checkfile,version=2)
    reuse <- identical(unname(tools::md5sum(checkfile)),fit$audit$data_hash)
    unlink(checkfile)
  }
  if (!reuse) fit <- withCallingHandlers(recurcast(adapter$formula, adapter$data,
    time_unit=adapter$metadata$time_unit, event_types=adapter$event_types,
    nonrecurrent=adapter$nonrecurrent, breaks=if(is.null(adapter$breaks))c(0,Inf) else adapter$breaks,
    observation=if(is.null(adapter$coarsening))"exact" else "coarsened",
    coarsening=adapter$coarsening, control=control), warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
    })
  elapsed <- proc.time()[["elapsed"]] - began
  if (reuse && !is.null(previous_summary) && label %in% previous_summary$case)
    elapsed <- previous_summary$elapsed_seconds[match(label,previous_summary$case)]
  saveRDS(fit,file.path(cache_directory,paste0(label,".rds")))
  stats <- data.frame(case=label,n=fit$n,paths=fit$npaths,
    recurrent_events=adapter$metadata$recurrences,terminal_events=adapter$metadata$terminal_events,
    source_kind=adapter$metadata$source_kind,observation=fit$schema$observation,
    time_unit=fit$schema$time_unit,frailty_variance=fit$phi,log_likelihood=fit$loglik,
    estimated_parameters=length(fit$coefficients),convergence=fit$convergence$code,
    gradient_max=fit$convergence$gradient,min_information_eigenvalue=min(fit$convergence$information_eigenvalues),
    elapsed_seconds=elapsed,warnings=paste(unique(warnings),collapse=" | "),stringsAsFactors=FALSE)
  overview[[length(overview)+1L]] <<- stats
  write.csv(do.call(rbind,overview),file.path(output,"case_summary.csv"),row.names=FALSE)
  coef <- tidy(fit);coef$case <- label
  write.csv(coef,file.path(output,paste0(label,"_coefficients.csv")),row.names=FALSE)
  print(stats[,c("case","n","convergence","gradient_max","elapsed_seconds")])
  fit
}
collect_forecasts <- function(name, fit, pred, endpoint, landmark, horizon=1) {
  at <- pred[pred$target=="risk" & pred$horizon==horizon,,drop=FALSE]
  row <- match(as.character(at$unit),as.character(endpoint$id))
  stopifnot(!anyNA(row),all(endpoint$time[row]>landmark))
  evaluation <- validate(at$estimate,time=endpoint$time[row]-landmark,
    status=endpoint$status[row],horizon=horizon,threshold=.2)
  m <- evaluation$metrics;m$case<-name;m$landmark<-landmark;m$horizon<-horizon
  m$heldout_eligible<-nrow(at);m$cases<-evaluation$cases;m$controls<-evaluation$controls
  m$calibration_status<-evaluation$calibration_status
  metrics[[length(metrics)+1L]] <<- m
  write.csv(do.call(rbind,metrics),file.path(output,"heldout_metrics.csv"),row.names=FALSE)
  g <- split(seq_len(nrow(pred)),interaction(pred$horizon,pred$target,drop=TRUE))
  for (rows in g) curves[[length(curves)+1L]] <<- data.frame(case=name,
    landmark=landmark,horizon=pred$horizon[rows[1]],target=pred$target[rows[1]],
    mean=mean(pred$estimate[rows]),q25=unname(quantile(pred$estimate[rows],.25)),
    q75=unname(quantile(pred$estimate[rows],.75)),n=length(rows),time_unit=pred$time_unit[rows[1]],
    max_series_error_bound=max(pred$series_error_bound[rows]))
  write.csv(do.call(rbind,curves),file.path(output,"forecast_curve_summary.csv"),row.names=FALSE)
}
endpoint_exact <- function(data) {
  do.call(rbind,lapply(split(data,data$id),function(d) data.frame(id=d$id[1],time=max(d$stop),status=tail(d$terminal,1))))
}
exact_adapter <- function(data,formula,synthetic=FALSE) list(data=data,formula=formula,
  event_types=sort(unique(as.character(data$event_type[data$event_type!=0]))),nonrecurrent=character(),coarsening=NULL,
  metadata=list(time_unit=if(synthetic)attr(data,"time_unit") else "years",
    recurrences=sum(data$event_type!=0),terminal_events=sum(data$terminal),
    source_kind=if(synthetic)"synthetic demonstration" else "public clinical cohort"))
for (name in c("sim_perinatal","sim_oncology","sim_robotics")) {
  e <- new.env();data(list=name,package="recurcast",envir=e);d<-get(name,envir=e)
  a <- exact_adapter(d,~x,TRUE)
  ctl <- recurcast_control(maxit=600,tolerance=1e-10,quadrature=30,history_cap=2)
  full <- fit_model(name,a,ctl)
  ids <- sort(unique(d$id));test <- ids[seq_along(ids)%%5L==0L]
  train <- a;train$data<-d[!d$id%in%test,,drop=FALSE]
  train$metadata$recurrences<-sum(train$data$event_type!=0);train$metadata$terminal_events<-sum(train$data$terminal)
  fit <- fit_model(paste0(name,"_train"),train,ctl)
  heldout <- d[d$id%in%test,,drop=FALSE]
  pred <- predict(fit,heldout,landmark=.5,horizon=c(.25,.5,1,1.5,2),type="all")
  collect_forecasts(name,fit,pred,endpoint_exact(heldout),.5)
}
readmission <- rc_readmission()
readmission$breaks <- c(0,1,3,Inf)
readmission$metadata$source_kind <- "public clinical cohort"
readmission_full <- fit_model("readmission",readmission,recurcast_control(maxit=800,tolerance=1e-10,history_cap=3))
ids <- sort(unique(readmission$data$id));test <- ids[seq_along(ids)%%5L==0L]
train <- readmission;train$data<-readmission$data[!readmission$data$id%in%test,,drop=FALSE]
train$metadata$recurrences<-sum(train$data$event_type!="0");train$metadata$terminal_events<-sum(train$data$terminal)
readmission_train <- fit_model("readmission_train",train,recurcast_control(maxit=800,tolerance=1e-10,history_cap=3))
heldout<-readmission$data[readmission$data$id%in%test,,drop=FALSE]
pred <- predict(readmission_train,heldout,landmark=1,horizon=c(.25,.5,1,1.5,2),type="all")
collect_forecasts("readmission",readmission_train,pred,endpoint_exact(heldout),1)

colon <- rc_colon(nodes=4)
colon$metadata$terminal_events <- colon$metadata$source_deaths
colon$metadata$source_kind <- "public randomized clinical trial (prognostic analysis)"
ctl <- recurcast_control(maxit=800,tolerance=1e-10,history_cap=1)
colon_full <- fit_model("colon_node4",colon,ctl)
continuous <- rc_colon(nodes=4,covariate="nodes")
continuous$metadata$terminal_events <- continuous$metadata$source_deaths
continuous$metadata$source_kind <- colon$metadata$source_kind
colon_nodes <- fit_model("colon_nodes",continuous,ctl)
# Numerical quadrature audit at the same fitted parameters: no refitting and no
# target change. Retain per-subject maximum differences, not patient records.
subject_likelihood <- function(fit,adapter) {
  d <- as_trajectory(adapter$data,event_types="1",time_unit="years",observation="coarsened")
  dd <- getFromNamespace(".rc_design","recurcast")(d,fit$formula,fit$breaks,
    fit$control$history_cap,fit$nonrecurrent,fit$design$terms,fit$design$xlevels,fit$design$contrasts)
  X <- if(is.null(fit$design$active))dd$design else dd$design[,fit$design$active,drop=FALSE]
  z <- rc_loglik(fit$coefficients,dd$unit,dd$exposure,X,dd$event,information=FALSE)
  rc_coarsen(z,owner=adapter$coarsening$owner,
    log_weights=adapter$coarsening$log_weights)$subject_loglik
}
q8 <- rc_colon(nodes=8)
delta <- abs(subject_likelihood(colon_full,colon)-subject_likelihood(colon_full,q8))
qa <- data.frame(case="colon_node4",subjects=colon$metadata$patients,
  same_day_pairs=length(colon$metadata$same_day_ids),paths_4=colon$metadata$paths,
  paths_8=q8$metadata$paths,max_absolute_subject_loglik_difference=max(delta),
  mean_absolute_subject_loglik_difference=mean(delta),tolerance=1e-7,passed=max(delta)<1e-7)
write.csv(qa,file.path(output,"colon_quadrature_audit.csv"),row.names=FALSE)
raw <- getExportedValue("survival","colon")
ids<-sort(unique(raw$id));test<-ids[seq_along(ids)%%5L==0L]
training <- rc_colon(raw[!raw$id%in%test,],nodes=4)
training$metadata$terminal_events<-training$metadata$source_deaths
training$metadata$source_kind<-colon$metadata$source_kind
colon_train<-fit_model("colon_train",training,ctl)
# Landmark-specific observation fibre: erase every later event before creating
# quadrature paths. A whole-day landmark avoids a day interval straddling it.
day <- 365;landmark <- day/365.25
test_raw <- raw[raw$id%in%test,]
ep <- do.call(rbind,lapply(split(test_raw,test_raw$id),function(g) {
  r<-g[g$etype==1,];d<-g[g$etype==2,];end<-d$time;status<-d$status
  if(r$status==0 && r$time<end){end<-r$time;status<-0}
  data.frame(id=d$id,time=end/365.25,status=status)
}))
eligible<-ep$id[ep$time>landmark]
past<-test_raw[test_raw$id%in%eligible,]
future <- past$time>day
past$status[future]<-0;past$time[future]<-day
past$time[past$etype==2]<-day;past$status[past$etype==2]<-0
histories<-rc_colon(past,nodes=4)
pred <- rc_predict_paths(colon_train,histories$data,histories$coarsening$owner,
  histories$coarsening$log_weights,landmark,c(.25,.5,1,1.5,2),type="all")
pred$unit<-as.character(histories$metadata$source_ids[as.integer(pred$unit)])
collect_forecasts("colon_node4",colon_train,pred,ep,landmark)
saveRDS(list(readmission=readmission$metadata,colon=colon$metadata,
  colon_nodes=continuous$metadata,colon_landmark=histories$metadata),file.path(output,"adapter_metadata.rds"))

curves <- do.call(rbind,curves)
for (name in unique(curves$case)) {
  d<-curves[curves$case==name & curves$target=="risk",];d<-d[order(d$horizon),]
  draw <- function() {
    par(mar=c(4.5,4.5,2,1));plot(d$horizon,d$mean,type="n",ylim=c(0,max(.05,d$q75)*1.1),
      xlab=paste("Forecast horizon (",d$time_unit[1],")",sep=""),ylab="Predicted terminal-event probability",
      main=gsub("_"," ",name),bty="l")
    polygon(c(d$horizon,rev(d$horizon)),c(d$q25,rev(d$q75)),col="#D8E7F0",border=NA)
    lines(d$horizon,d$mean,lwd=2,col="#205072");points(d$horizon,d$mean,pch=16,col="#205072")
    legend("topleft",c("Held-out mean forecast","Between-person interquartile range"),
      col=c("#205072","#D8E7F0"),lwd=c(2,8),bty="n",cex=.8)
  }
  pdf(file.path(output,paste0(name,"_forecast.pdf")),width=6.8,height=4.6);draw();dev.off()
  png(file.path(output,paste0(name,"_forecast.png")),width=2040,height=1380,res=300);draw();dev.off()
}
writeLines(c("Five software case studies; reproducible diagnostics, not efficacy claims.",
  "Training/test split: every fifth sorted source identifier is held out. The split is deterministic, not an external validation cohort.",
  "All five full datasets are fitted; hold-out metrics use separate training fits.",
  "Three 160-unit datasets are synthetic demonstrations, not maternal, oncology or engineering benchmark observations.",
  "Colon: all 929 retained; same-day records integrated over legal recurrence-before-death triangles. Nonrepeatable recurrence submodel.",
  "Readmission: 861 interval rows represent 458 recurrences, not 861 events. Baseline Charlson is fixed at first interval.",
  "Patient-level local fit objects are under build/case_fits and are not manuscript tables or bundled package data.",
  "Forecast ribbons are between-person interquartile ranges, not confidence bands or coverage guarantees.",
  "Reported probabilities are prognostic. No intervention effect, mortality reduction or clinical deployment safety is established.",
  "IPCW validation uses marginal censoring KM; its validity requires independent censoring for this evaluation.",
  "Nonregular fit flags and undefined calibration metrics are retained."),file.path(output,"INTERPRETATION.txt"))
capture.output(sessionInfo(),file=file.path(output,"sessionInfo.txt"))
cat("All case-study computations completed.\n")
