targeting_boundary_fixture <- function() {
 d<-read.csv(test_path("fixtures","continuous_target_subjects.csv"))
 list(data=d,spec=list(horizon=1,censor="censored",strata="x",folds="fold",grid=16,working="pooled"))
}
test_that("working models and supplied training nuisances generate the same explicit experiment", {
 a<-targeting_boundary_fixture();d<-a$data;spec<-a$spec
 spec$working<-"loglinear"
 z<-recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=spec)
 expect_true(all(is.finite(z$EIF)));expect_true(all(z$initial_rates>0))
 supplied<-spec;supplied$initial_rates<-z$initial_rates;supplied$censor_rates<-z$censor_rates
 zz<-recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=supplied)
 expect_equal(zz$estimate,z$estimate,tolerance=1e-12)
 expect_equal(zz$EIF,z$EIF,tolerance=1e-11)
 offset<-spec;offset$working<-"pooled";state<-z$states
 offset$rate_offset<-rbind(exp(.3*state[,2]),exp(.35*state[,1]),exp(.45*(state[,1]+state[,2])/2))
 off<-recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=offset)
 expect_lt(max(abs(off$EIF_mean)),1e-7/sqrt(off$n))
 expect_output(print(off),"targeted common law")
 expect_error(logLik(off),"no supplied maximized likelihood",class="recurcast_error")
 expect_error(AIC(off),class="recurcast_error")
 for(s in list(modifyList(spec,list(working="unknown")),modifyList(spec,list(rate_offset=matrix(1,2,2))),
   modifyList(spec,list(initial_rates=array(1,c(3,9,4,10)))),modifyList(spec,list(censor_rates=1)),
   modifyList(spec,list(censor_rates=rep(100,10))),modifyList(spec,list(grid=5.5)),
   modifyList(spec,list(history_cap=7)),modifyList(spec,list(folds=rep(NA,500)))))
   expect_error(recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=s),class="recurcast_error")
})
test_that("nonrepeatable mark eligibility is preserved by the exponential fluctuation", {
 a<-targeting_boundary_fixture();d<-a$data;spec<-a$spec
 for(i in unique(d$id)){ix<-which(d$id==i & d$event_type=="recurrent1");if(length(ix)>1)d$event_type[ix[-1]]<-"none"}
 spec$nonrecurrent<-"recurrent1"
 z<-suppressWarnings(recurcast(~x,d,time_unit="year",targeting="continuous",event_types=c("recurrent1","recurrent2"),targeting_spec=spec))
 expect_true(all(z$rates[1,z$states[,1]>0,,]==0))
 expect_equal(z$target_limits[2,1],1)
 expect_error(recurcast(~x,a$data,time_unit="year",targeting="continuous",targeting_spec=spec),"non-repeatable",class="recurcast_error")
 spec$nonrecurrent<-"absent"
 expect_error(recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=spec),class="recurcast_error")
})
test_that("observation guards prevent incompatible censoring and baseline strata", {
 a<-targeting_boundary_fixture();d<-a$data;spec<-a$spec
 bad<-d;bad$censored[1]<-2
 expect_error(recurcast(~x,bad,time_unit="year",targeting="continuous",targeting_spec=spec),class="recurcast_error")
 bad<-d;index<-which(duplicated(bad$id,fromLast=TRUE))[1];bad$censored[index]<-1
 expect_error(recurcast(~x,bad,time_unit="year",targeting="continuous",targeting_spec=spec),class="recurcast_error")
 bad<-d;bad$observed<-TRUE;bad$observed[1]<-FALSE
 expect_error(recurcast(~x,bad,observed="observed",time_unit="year",targeting="continuous",targeting_spec=spec),class="recurcast_error")
 bad<-d;bad$x<-bad$id
 expect_error(recurcast(~x,bad,time_unit="year",targeting="continuous",targeting_spec=spec),"30",class="recurcast_error")
 bad<-d;bad$x<-bad$fold
 expect_error(recurcast(~x,bad,time_unit="year",targeting="continuous",targeting_spec=spec),"training fold",class="recurcast_error")
 s<-spec;s$strata<-"absent"
 expect_error(recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=s),class="recurcast_error")
 s<-spec;s$complete_followup<-TRUE
 expect_error(recurcast(~x,d,time_unit="year",targeting="continuous",targeting_spec=s),"uncensored",class="recurcast_error")
})
test_that("finite-path schemas do not silently ignore predictors or duplicate subjects", {
 d<-data.frame(id=1:40,path=rep(1:2,20),recorded=1,x=rep(0:1,each=20));f<-matrix(c(0,1),2)
 spec<-list(path="path",recorded="recorded",features=f)
 z<-recurcast(~1,d,time_unit="window",targeting="finite_path",targeting_spec=spec)
 expect_identical(z$target_names,"target_1")
 expect_equal(sum(z$mixture),1,tolerance=1e-13)
 expect_error(recurcast(~x,d,time_unit="window",targeting="finite_path",targeting_spec=modifyList(spec,list(strata="id"))),class="recurcast_error")
 expect_error(recurcast(~1,d,time_unit="window",targeting="finite_path",targeting_spec=list(path="path")),class="recurcast_error")
 expect_error(recurcast(~1,d,time_unit="window",targeting="finite_path",targeting_spec=modifyList(spec,list(recorded="missing"))),class="recurcast_error")
 bad<-d;bad$id[2]<-bad$id[1]
 expect_error(recurcast(~1,bad,time_unit="window",targeting="finite_path",targeting_spec=spec),class="recurcast_error")
 d$path[]<-1;d$recorded[1:10]<-0
 rare<-recurcast(~1,d,time_unit="window",targeting="finite_path",targeting_spec=spec)
 expect_false(rare$convergence$valid)
 projected<-predict(rare)
 expect_identical(projected$interval,"bounded_coarsening_projection")
 expect_equal(projected$lower,0)
 expect_gt(projected$upper,.25)
 expect_error(recurcast(~1,d,time_unit="window",targeting="finite_path",targeting_spec=modifyList(spec,list(complete_followup=TRUE))),class="recurcast_error")
})
test_that("dynamic conditioning includes past recurrence marks and rejects unknown histories", {
 a<-targeting_boundary_fixture();z<-recurcast(~x,a$data,time_unit="year",targeting="continuous",targeting_spec=a$spec)
 d<-a$data;eligible<-vapply(split(d,d$id),function(q)any(q$event_type!="none"&q$stop<.5)&&max(q$stop)>=.5&&!any(q$terminal==1&q$stop<=.5),logical(1))
 id<-as.integer(names(eligible)[which(eligible)[1]]);subject<-d[d$id==id,,drop=FALSE]
 p<-predict(z,subject,landmark=.5,horizon=.2,type="burden",interval="none")
 expect_equal(nrow(p),2L);expect_true(all(p$estimate>=0))
 w<-attr(p,"fold_posterior")[[1]];expect_equal(sum(w),1,tolerance=1e-14)
 x<-z$subject_stratum[id];expect_gt(max(abs(w-z$fold_prior[x,])),1e-8)
 expect_error(predict(z,subject,landmark=.5,horizon=.2,unused=1),class="recurcast_error")
 unknown<-subject;unknown$x<-2
 expect_error(predict(z,unknown,landmark=.5,horizon=.2),"stratum",class="recurcast_error")
 dead<-d[d$id==d$id[which(d$terminal==1&d$stop<.9)[1]],,drop=FALSE]
 expect_error(predict(z,dead,landmark=.9,horizon=.1),"eligible",class="recurcast_error")
 expect_error(predict(z,type="absent"),class="recurcast_error")
 expect_error(predict(z,interval="absent"),class="recurcast_error")
 expect_error(predict(z,horizon=0),class="recurcast_error")
})
