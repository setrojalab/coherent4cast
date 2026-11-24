#' Simulate a shared-gamma recurrent and terminal event history
#' @export
rc_simulate <- function(n=100L,horizon=3,rates=c(.6,.25,.12),beta=c(-.2,-.1,.15),
    history=c(.2,.2,.25),dispersion=.45,censor_rate=.1,seed=20261122L,time_unit="years",history_cap=2L) {
  scalar<-function(x,lower=0,integer=FALSE) is.numeric(x)&&length(x)==1L&&is.finite(x)&&x>=lower&&(!integer||x==floor(x))
  if(!scalar(n,1,TRUE)||!scalar(horizon,.Machine$double.xmin)||!scalar(dispersion)||!scalar(censor_rate)||
     !scalar(history_cap,0,TRUE)||!scalar(seed,0,TRUE)||length(time_unit)!=1L||is.na(time_unit)||!nzchar(time_unit)||length(rates)<1||
     any(!is.finite(c(rates,beta,history,dispersion,censor_rate)))||any(rates<0)||any(rates==0)||dispersion<0||censor_rate<0||
     length(beta)!=length(rates)||length(history)!=length(rates)||history_cap<0) .rc_error("simulation","Invalid generating parameters.")
  previous_kind<-RNGkind();had_seed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) previous_seed<-get(".Random.seed",envir=.GlobalEnv)
  on.exit({do.call(RNGkind,as.list(previous_kind));if(had_seed)assign(".Random.seed",previous_seed,envir=.GlobalEnv)else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  RNGkind("L'Ecuyer-CMRG");set.seed(seed)
  K<-length(rates)-1L;rows<-list();index<-0L
  for(i in seq_len(n)) {
    x<-stats::rbinom(1,1,.5);z<-if(dispersion==0)1 else stats::rgamma(1,shape=1/dispersion,rate=1/dispersion)
    t<-0;counts<-numeric(K);jumps<-0L
    while(t<horizon) {
      rr<-rates*exp(beta*x+history*sum(pmin(counts,history_cap)))*z
      if(any(!is.finite(rr))||sum(rr)+censor_rate<=0) .rc_error("simulation_rate","Generating hazards are not finite and positive.")
      delta<-stats::rexp(1,sum(rr)+censor_rate);end<-min(t+delta,horizon);mark<-0L;terminal<-0L;censor<-FALSE
      if(t+delta<horizon) {
        event<-sample.int(length(rr)+1L,1,prob=c(rr,censor_rate))
        if(event<=K){mark<-event;counts[event]<-counts[event]+1}else if(event==K+1L)terminal<-1L else censor<-TRUE
      }
      index<-index+1L;rows[[index]]<-data.frame(id=i,start=t,stop=end,event_type=mark,terminal=terminal,x=x,observed=TRUE)
      t<-end;jumps<-jumps+1L
      if(terminal==1L||censor)break
      if(jumps>10000L) .rc_error("simulation_explosion","More than 10,000 intervals for one unit; the simulation is stopped without truncating its history.")
    }
  }
  d<-do.call(rbind,rows);rownames(d)<-NULL
  attr(d,"synthetic")<-TRUE;attr(d,"time_unit")<-time_unit
  attr(d,"truth")<-list(rates=rates,beta=beta,history=history,dispersion=dispersion,censor_rate=censor_rate,horizon=horizon,seed=seed,RNG="L'Ecuyer-CMRG",history_cap=history_cap)
  d
}
