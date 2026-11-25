# Reproduce the R-only 0.10.1 wrappers and audit export without changing kernels.
arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
root<-if(length(arg))normalizePath(file.path(dirname(sub("^--file=","",arg[1])),".."),winslash="/")else normalizePath(".",winslash="/")
local_library<-Sys.getenv("COHERENT4CAST_LIBRARY","")
if(nzchar(local_library)).libPaths(c(local_library,.libPaths()))
library(coherent4cast)
outdir<-Sys.getenv("COHERENT4CAST_STANDARD_OUTPUT",file.path(root,"replication/results/standard_regimes"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

d<-data.frame(id=1:12,time=c(.2,.3,.4,.5,.6,.8,1,1.2,1.3,1.5,1.7,2),
  status=c(1,2,0,1,0,2,1,0,2,0,1,0),x=rep(0:1,6))
write.csv(d,file.path(outdir,"toy_input.csv"),row.names=FALSE)
comp<-rc_competing(~1,d,unit="id",time_unit="year",frailty_diagnostics=FALSE)
pred<-predict(comp,newdata=data.frame(id=101),landmark=.5,horizon=c(.5,1,2),type="all")
write.csv(pred,file.path(outdir,"toy_competing_forecasts.csv"),row.names=FALSE)
rate<-as.numeric(table(factor(d$status,levels=c(1,2))))/sum(d$time)
oracle<-do.call(rbind,lapply(c(.5,1,2),function(h){r<-sum(rate);data.frame(horizon=h,
  target=c("cif_1","cif_2","survival","risk","rmst"),
  oracle=c(rate/r*(1-exp(-r*h)),exp(-r*h),1-exp(-r*h),(1-exp(-r*h))/r))}))
comparison<-merge(as.data.frame(pred),oracle,by=c("horizon","target"))
comparison$absolute_difference<-abs(comparison$estimate-comparison$oracle)
stopifnot(max(comparison$absolute_difference)<1e-7)
write.csv(comparison,file.path(outdir,"competing_closed_form_parity.csv"),row.names=FALSE)
write.csv(tidy(comp),file.path(outdir,"competing_coefficients.csv"),row.names=FALSE)
val<-validate(comp,time=c(.2,.3,2,2),status=c(1,2,0,0),horizon=1,cause="1",
  predicted=c(.8,.7,.2,.1),threshold=.5)
write.csv(val$metrics,file.path(outdir,"competing_validation_hand_example.csv"),row.names=FALSE)
write.csv(as.data.frame.matrix(val$confusion),file.path(outdir,"competing_validation_confusion.csv"))

survdata<-transform(d,status=as.integer(status>0))
survfit<-rc_surv(survival::Surv(time,status)~1,survdata,unit="id",time_unit="year",frailty_diagnostics=FALSE)
write.csv(predict(survfit,newdata=data.frame(id=101),horizon=c(.5,1,2),type="all"),
  file.path(outdir,"toy_survival_forecasts.csv"),row.names=FALSE)
audit<-rc_audit_report(comp,predictions=pred,
  source_files=file.path(root,c("package/R/wrappers.R","package/inst/include/libcoherent4cast/sparse_uniformize.hpp")),
  format="list")
audit_json<-rc_audit_report(comp,predictions=pred,source_files=file.path(root,"package/R/wrappers.R"),
  format="json",file=file.path(outdir,"toy_audit.json"))
audit_text<-rc_audit_report(comp,predictions=pred,format="text",file=file.path(outdir,"toy_audit.txt"))

# Public survival data; split before fitting. No clinical effect is inferred.
lung<-survival::lung
lung$id<-seq_len(nrow(lung));lung$death<-as.integer(lung$status==2)
train<-lung[lung$id%%5!=0,];test<-lung[lung$id%%5==0,]
model<-rc_surv(survival::Surv(time,death)~age+factor(sex),train,unit="id",time_unit="day",
  breaks=c(0,180,365,Inf),frailty_diagnostics=FALSE)
risk<-predict(model,newdata=test,landmark=0,horizon=365,type="risk")
index<-match(as.integer(risk$unit),test$id)
performance<-validate(model,time=test$time[index],status=test$death[index],horizon=365,predicted=risk$estimate,threshold=.5)
write.csv(tidy(model),file.path(outdir,"lung_coefficients.csv"),row.names=FALSE)
write.csv(glance(model),file.path(outdir,"lung_fit.csv"),row.names=FALSE)
write.csv(risk,file.path(outdir,"lung_heldout_forecasts.csv"),row.names=FALSE)
write.csv(performance$metrics,file.path(outdir,"lung_heldout_metrics.csv"),row.names=FALSE)
write.csv(data.frame(n=nrow(lung),training=nrow(train),heldout=nrow(test),eligible=nrow(risk),horizon=365,
  terminal_events_training=sum(train$death),terminal_events_test=sum(test$death)),
  file.path(outdir,"lung_design.csv"),row.names=FALSE)
writeLines(c("Toy data are constructed to check a closed-form likelihood/CIF calculation, not a sampling-performance study.",
  "Lung data use survival::lung; death=1(status==2), deterministic held-out IDs divisible by5, 365-day horizon.",
  "Held-out metrics use marginal IPCW; no external clinical, causal or regulatory validation is inferred.",
  "Competing events are absorbing alternatives; cumulative incidence is not 1 minus cause-specific survival.",
  "Audit SHA256 describes serialization-v3 data bytes and explicitly supplied source files; it is not a regulatory certification."),file.path(outdir,"INTERPRETATION.txt"))
writeLines(capture.output(sessionInfo()),file.path(outdir,"sessionInfo.txt"))
cat("Maximum competing closed-form error:",max(comparison$absolute_difference),"\n")
print(performance$metrics,row.names=FALSE)
cat("Outputs:",outdir,"\n")
