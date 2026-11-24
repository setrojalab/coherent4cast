.rc_linear_coef <- function(object) {
  if(is.null(object$design$active))return(head(object$coefficients,ncol(object$design$design)))
  b<-numeric(object$design$full_size);b[object$design$active]<-object$coefficients[seq_along(object$design$active)];b
}

.rc_gamma_nodes <- function(shape, rate, n) {
  if(!is.finite(shape)||!is.finite(rate)||shape<=0||rate<=0) .rc_error("frailty","Invalid posterior gamma parameters.")
  i<-seq_len(n);J<-diag(2*i-1+(shape-1));off<-sqrt(seq_len(n-1)*(seq_len(n-1)+shape-1))
  J[cbind(seq_len(n-1),seq_len(n-1)+1)]<-off;J[cbind(seq_len(n-1)+1,seq_len(n-1))]<-off
  eig<-eigen(J,symmetric=TRUE);list(nodes=eig$values/rate,weights=eig$vectors[1,]^2)
}

.rc_states <- function(K, cap, maximum) {
  if((cap+1)^K+1>maximum) .rc_error("state_space","The retained count state exceeds max_states. Reduce the declared cap or provide a smaller scientific model.")
  if(K==0L)return(matrix(numeric(),1,0))
  as.matrix(expand.grid(rep(list(0:cap),K)))
}

.rc_generator <- function(object, states, covariates, bin, frailty) {
  K<-object$design$K;S<-nrow(states);J<-object$design$base_size;P<-object$design$block_size
  Q<-matrix(0,S+1,S+1);rewards<-matrix(0,S+1,K+1)
  base<-numeric(J);base[bin]<-1
  for(i in seq_len(S)) {
    coefficients<-.rc_linear_coef(object)
    v<-c(base,covariates,states[i,]);rates<-vapply(seq_len(K+1),function(k)exp(sum(v*coefficients[(k-1)*P+seq_len(P)]))*frailty,numeric(1))
    if(any(!is.finite(rates))) .rc_error("forecast_overflow","Forecast rates overflow; rescale covariates or shorten the model scale.")
    if(K)for(k in seq_len(K)) {
      if(k %in% object$nonrecurrent && states[i,k]>0)rates[k]<-0
      nextstate<-states[i,];nextstate[k]<-min(nextstate[k]+1,object$control$history_cap)
      target<-which(rowSums(abs(sweep(states,2,nextstate,"-")))==0)[1]
      if(target!=i)Q[i,target]<-Q[i,target]+rates[k]
      rewards[i,k]<-rates[k]
    }
    Q[i,S+1]<-rates[K+1];Q[i,i]<- -sum(Q[i,-i]);rewards[i,K+1]<-1
  }
  list(Q=Q,reward=rewards)
}

.rc_generator_sparse <- function(object,states,covariates,bin,frailty) {
  K<-object$design$K;S<-nrow(states);J<-object$design$base_size;P<-object$design$block_size
  base<-numeric(J);base[bin]<-1;coefficients<-.rc_linear_coef(object)
  block<-cbind(matrix(rep(c(base,covariates),each=S),S),states)
  rates<-sapply(seq_len(K+1L),function(k)exp(as.vector(block%*%coefficients[(k-1L)*P+seq_len(P)]))*frailty)
  if(any(!is.finite(rates))) .rc_error("forecast_overflow","Forecast rates overflow; rescale covariates or shorten the model scale.")
  if(!is.matrix(rates))rates<-matrix(rates,S,K+1L)
  from<-to<-integer();rate<-numeric();rewards<-matrix(0,S+1L,K+1L)
  if(K)for(k in seq_len(K)){
    if(k %in% object$nonrecurrent)rates[states[,k]>0,k]<-0
    move<-which(states[,k]<object$control$history_cap & rates[,k]>0)
    from<-c(from,move);to<-c(to,move+(object$control$history_cap+1)^(k-1));rate<-c(rate,rates[move,k]);rewards[seq_len(S),k]<-rates[,k]
  }
  from<-c(from,seq_len(S));to<-c(to,rep(S+1L,S));rate<-c(rate,rates[,K+1L]);rewards[seq_len(S),K+1L]<-1
  list(Q=rc_generator_csr(S+1L,from,to,rate),reward=rewards)
}

#' Forecast a fitted common-law trajectory at a landmark
#' @export
predict.recurcast <- function(object,newdata=NULL,landmark,horizon,type=c("risk","burden","survival","rmst","all"),
                              quadrature=object$control$quadrature,
                              covariate_evolution=c("freeze","locf","joint_submodel"),newdata_future=NULL,
                              covariate_model=NULL,covariate_draws=100L,covariate_grid=16L,seed=1L,covariate_history=NULL,
                              newdata_observation=c("exact","coarsened"),coarsening_control=list(),...) {
  type<-match.arg(type)
  covariate_evolution<-match.arg(covariate_evolution)
  newdata_observation<-match.arg(newdata_observation)
  if(length(list(...))) .rc_error("unused_argument","An unsupported prediction argument was supplied.")
  if(missing(landmark)||length(landmark)!=1L||!is.finite(landmark)||landmark<0) .rc_error("landmark","Declare one finite nonnegative landmark.")
  if(missing(horizon)||!length(horizon)||any(!is.finite(horizon))||any(horizon<=0)) .rc_error("horizon","Declare positive finite forecast horizons.")
  if(!is.numeric(quadrature)||length(quadrature)!=1L||!is.finite(quadrature)||quadrature<2||quadrature!=floor(quadrature)) .rc_error("quadrature","quadrature must be an integer of at least two.")
  if(newdata_observation=="coarsened")return(rc_predict_coarsened(object,newdata,landmark,horizon,type=type,quadrature=quadrature,coarsening_control=coarsening_control,
    covariate_evolution=covariate_evolution,newdata_future=newdata_future,covariate_model=covariate_model,covariate_draws=covariate_draws,covariate_grid=covariate_grid,covariate_history=covariate_history,seed=seed))
  if(object$schema$observation=="coarsened"&&is.null(newdata)) .rc_error("coarsened_forecast","Supply the history available at the landmark as newdata; complete training-data fibres contain future information and cannot be reused as landmark histories.")
  schema<-object$schema;d<-if(is.null(newdata))object$data else as_trajectory(newdata,schema$unit,schema$start,schema$stop,schema$event_type,schema$terminal,schema$observed,schema$event_types,schema$time_unit,"exact")
  if(length(attr(d,"rc_schema")$delayed_entry)) .rc_error("left_truncation","Forecast histories must begin at time zero; unobserved pre-entry history and frailty selection require an explicit observation model.")
  if(any(!d$.rc_observed)||length(attr(d,"rc_schema")$gaps)) .rc_error("unrecorded_history","A point forecast needs an observation model for unrecorded intervals.")
  known<-d[d$.rc_start<=landmark,,drop=FALSE]
  incomplete<-known$.rc_stop>landmark
  known$.rc_stop[incomplete]<-landmark;known$.rc_mark[incomplete]<-0L;known$.rc_terminal[incomplete]<-0L
  attr(known,"rc_schema")<-attr(d,"rc_schema")
  dd<-.rc_design(known,object$formula,object$breaks,object$control$history_cap,object$nonrecurrent,object$design$terms,object$design$xlevels,object$design$contrasts)
  dt<-pmax(0,pmin(dd$data$.rc_stop,landmark)-dd$data$.rc_start)
  linear_coefficients<-.rc_linear_coef(object)
  eta<-as.vector(dd$design%*%linear_coefficients)
  nrows<-nrow(dd$data);K<-dd$K;exposure<-rep(dt,K+1L)
  exposure[dd$exposure==0]<-0
  haz<-tapply(exposure*exp(eta),dd$unit,sum)
  before<-dd$data$.rc_stop<=landmark
  jumps<-tapply(as.integer(dd$data$.rc_mark>0 & before)+dd$data$.rc_terminal*before,dd$data$.rc_unit,sum)
  states<-.rc_states(K,object$control$history_cap,object$control$max_states);S<-nrow(states)
  output<-list();output_draws<-list();rejected<-character()
  for(rows in split(seq_len(nrow(d)),d$.rc_unit)) {
    ident<-as.character(d[[schema$unit]][rows[1]]);u<-d$.rc_unit[rows[1]]
    if(any(d$.rc_terminal[rows]==1 & d$.rc_stop[rows]<=landmark)) {rejected<-c(rejected,ident);next}
    if(max(d$.rc_stop[rows])<landmark) {rejected<-c(rejected,ident);next}
    observed_rows<-rows[d$.rc_start[rows]<=landmark]
    if(!length(observed_rows)) {rejected<-c(rejected,ident);next}
    current<-tail(observed_rows,1)
    mf<-stats::model.frame(object$design$terms,d[current,,drop=FALSE],na.action=stats::na.fail,xlev=object$design$xlevels)
    mm<-stats::model.matrix(object$design$terms,mf,contrasts.arg=object$design$contrasts)
    covariates<-mm[,object$design$covariate_columns,drop=TRUE]
    counts<-if(K)vapply(seq_len(K),function(k)sum(d$.rc_mark[rows]==k & d$.rc_stop[rows]<=landmark),numeric(1)) else numeric()
    state<-if(K)which(rowSums(abs(sweep(states,2,pmin(counts,object$control$history_cap),"-")))==0)[1] else 1L
    probability<-numeric(S+1L);probability[state]<-1
    quadr<-.rc_posterior_frailty(object,jumps[as.character(u)],haz[as.character(u)],as.integer(quadrature))
    for(delta in horizon) {
      scenarios<-.rc_future_scenarios(object,d,rows,current,ident,landmark,delta,covariate_evolution,newdata_future,covariate_model,covariate_draws,covariate_grid,seed,covariate_history)
      cuts<-scenarios$cuts;B<-length(scenarios$covariates)
      prob<-numeric(S+1L);reward<-numeric(K+1L);err<-0;reward_err<-numeric(K+1L);roundoff<-0
      compression<-numeric(K+3L)
      samples<-matrix(0,B,K+3L)
      for(b in seq_len(B)) {
        bp<-numeric(S+1L);br<-numeric(K+1L)
        for(g in seq_along(quadr$nodes)) {
          gen<-lapply(seq_len(length(cuts)-1L),function(j).rc_generator_sparse(object,states,scenarios$covariates[[b]][j,],findInterval(cuts[j],object$breaks,all.inside=TRUE),quadr$nodes[g]))
          if(!is.null(quadr$compression_index)&&g==quadr$compression_index) {
            burden_bound<-numeric(K)
            if(K)for(j in seq_along(gen))burden_bound<-burden_bound+diff(cuts)[j]*apply(gen[[j]]$reward[,seq_len(K),drop=FALSE],2,max)
            compression<-compression+quadr$compression_mass*c(1,1,delta,2*burden_bound)/B
          }
          ans<-rc_propagate_sparse(probability,lapply(gen,`[[`,"Q"),diff(cuts),lapply(gen,`[[`,"reward"),tolerance=object$control$propagation_tolerance)
          bp<-bp+quadr$weights[g]*ans$probability;br<-br+quadr$weights[g]*ans$reward
          err<-err+quadr$weights[g]*ans$error_bound/B;reward_err<-reward_err+quadr$weights[g]*ans$reward_error_bound/B;roundoff<-roundoff+quadr$weights[g]*ans$roundoff_bound/B
        }
        prob<-prob+bp/B;reward<-reward+br/B;samples[b,]<-c(bp[S+1L],sum(bp[seq_len(S)]),br[K+1L],br[seq_len(K)])
      }
      values<-c(risk=prob[S+1],survival=sum(prob[seq_len(S)]),rmst=reward[K+1],if(K)stats::setNames(reward[seq_len(K)],paste0("burden_",schema$event_types)))
      wanted<-switch(type,risk="risk",survival="survival",rmst="rmst",burden=names(values)[grepl("^burden_",names(values))],all=names(values))
      errors<-c(risk=err,survival=err,rmst=reward_err[K+1],if(K)stats::setNames(reward_err[seq_len(K)],paste0("burden_",schema$event_types)))
      mcse<-stats::setNames(if(B>1L)apply(samples,2,stats::sd)/sqrt(B)else rep(0,K+3L),names(values))
      names(compression)<-names(values)
      if(B>1L)output_draws[[length(output_draws)+1L]]<-samples[,match(wanted,names(values)),drop=FALSE]
      output[[length(output)+1L]]<-data.frame(unit=ident,landmark=landmark,horizon=delta,target=wanted,estimate=unname(values[wanted]),
           time_unit=schema$time_unit,series_error_bound=unname(errors[wanted]),probability_roundoff_bound=unname(roundoff),
           frailty_compression_bound=unname(compression[wanted]),
           covariate_mcse=unname(mcse[wanted]),covariate_evolution=covariate_evolution,stringsAsFactors=FALSE)
    }
  }
  if(!length(output)) .rc_error("eligibility","No unit is observed and terminal-event-free at this landmark.")
  ans<-do.call(rbind,output);rownames(ans)<-NULL;class(ans)<-c("recurcast_prediction","rc_forecast","data.frame")
  attr(ans,"excluded_units")<-rejected
  if(length(output_draws))attr(ans,"covariate_draw_forecasts")<-do.call(cbind,output_draws)
  attr(ans,"assumptions")<-c(paste("Future covariate rule:",covariate_evolution),"capped count state determines future intensity",
    "frailty quadrature and covariate-grid errors are separate from uniformization and Monte Carlo errors",if(covariate_evolution=="joint_submodel")covariate_model$assumptions)
  ans
}

#' Retain a posterior mixture over compatible complete histories
#' @export
rc_filter_paths <- function(loglik,log_weights,counts,cumulative_hazard,dispersion,owner=rep(1L,length(loglik))) {
  n<-length(loglik)
  if(!n || !is.numeric(dispersion) || length(dispersion)!=1L || !is.finite(dispersion) || dispersion<=0 ||
     any(vapply(list(log_weights,counts,cumulative_hazard,owner),length,integer(1))!=n)||
     any(!is.finite(c(loglik,log_weights,counts,cumulative_hazard,owner))) || any(counts<0) || any(counts!=floor(counts)) ||
     any(cumulative_hazard<0) || any(owner<1) || any(owner!=floor(owner)))
    .rc_error("path_filter","Supply finite aligned path quantities and positive dispersion.")
  out<-data.frame(owner=owner,shape=1/dispersion+counts,rate=1/dispersion+cumulative_hazard,weight=0)
  for(rows in split(seq_len(n),owner)) {z<-loglik[rows]+log_weights[rows];w<-exp(z-max(z));out$weight[rows]<-w/sum(w)}
  out
}
