.rc_softmax <- function(x) {z<-exp(x-max(x));z/sum(z)}
.rc_normal_cache <- new.env(parent=emptyenv())
.rc_normal_nodes <- function(n) {
  key<-as.character(n);if(exists(key,.rc_normal_cache,inherits=FALSE))return(get(key,.rc_normal_cache))
  J<-matrix(0,n,n);i<-seq_len(n-1L);J[cbind(i,i+1L)]<-sqrt(i);J[cbind(i+1L,i)]<-sqrt(i)
  e<-eigen(J,symmetric=TRUE);out<-list(nodes=e$values,weights=e$vectors[1,]^2);assign(key,out,.rc_normal_cache);out
}
.rc_frailty_spec <- function(family,control,points=2L) {
  if(family %in% c("gamma","lognormal"))return(list(family=family,points=control$quadrature,
    start=if(is.null(control$phi_fixed))log(control$start_phi)else numeric(),
    lower=if(is.null(control$phi_fixed))control$log_phi_bounds[1]else numeric(),
    upper=if(is.null(control$phi_fixed))control$log_phi_bounds[2]else numeric(),
    names=if(is.null(control$phi_fixed))"log_frailty_variance"else character(),fixed=control$phi_fixed))
  if(!is.null(control$phi_fixed)) .rc_error("frailty","phi_fixed applies to gamma and lognormal families, not estimated discrete support.")
  if(length(points)!=1L||!is.finite(points)||points<2||points>8||points!=floor(points)) .rc_error("frailty","frailty_points must be an integer between 2 and 8.")
  q<-as.integer(points);list(family=family,points=q,start=c(rep(0,q-1),seq(-1.2,-.15,length.out=q-1)),
    lower=c(rep(-8,q-1),rep(-4,q-1)),upper=c(rep(8,q-1),rep(4,q-1)),
    names=c(paste0("frailty_weight_logit_",seq_len(q-1)),paste0("frailty_location_",seq_len(q-1))),fixed=NULL)
}
.rc_frailty_nodes <- function(eta,spec,quadrature=spec$points,derivatives=TRUE) {
  if(spec$family=="lognormal") {
    phi<-if(is.null(spec$fixed))exp(eta[1])else spec$fixed
    if(phi==0)return(list(h=0,g=0,dh=matrix(0,1,0),dg=matrix(0,1,0),d2h=array(0,c(1,0,0)),d2g=array(0,c(1,0,0)),variance=0))
    z<-.rc_normal_nodes(quadrature);v<-log1p(phi);s<-sqrt(v);h<-z$nodes*s-v/2
    q<-if(is.null(spec$fixed))1L else 0L;dh<-dg<-matrix(0,quadrature,q);d2h<-d2g<-array(0,c(quadrature,q,q))
    if(q){u<-phi/(1+phi);du<-u*(1-u);dh[,1]<-z$nodes*u/(2*s)-u/2;d2h[,1,1]<-z$nodes*(du/(2*s)-u*u/(4*s^3))-du/2}
    return(list(h=h,g=log(z$weights),dh=dh,dg=dg,d2h=d2h,d2g=d2g,variance=phi))
  }
  L<-spec$points;q<-L-1L;a<-c(eta[seq_len(q)],0);v<-c(eta[q+seq_len(q)],0)
  w<-.rc_softmax(a);tt<-.rc_softmax(a+v);ls<-function(x)max(x)+log(sum(exp(x-max(x))))
  h<-v-ls(a+v)+ls(a);g<-log(w);P<-2*q
  dh<-dg<-matrix(0,L,P);d2h<-d2g<-array(0,c(L,P,P));W<-diag(w)-tcrossprod(w);T<-diag(tt)-tcrossprod(tt)
  for(r in seq_len(L)){
    ind<-as.numeric(seq_len(q)==r);dh[r,seq_len(q)]<-w[seq_len(q)]-tt[seq_len(q)]
    dh[r,q+seq_len(q)]<-ind-tt[seq_len(q)];dg[r,seq_len(q)]<-ind-w[seq_len(q)]
    d2h[r,seq_len(q),seq_len(q)]<-(W-T)[seq_len(q),seq_len(q)]
    d2h[r,seq_len(q),q+seq_len(q)]<-d2h[r,q+seq_len(q),seq_len(q)]<-d2h[r,q+seq_len(q),q+seq_len(q)]<- -T[seq_len(q),seq_len(q)]
    d2g[r,seq_len(q),seq_len(q)]<- -W[seq_len(q),seq_len(q)]
  }
  list(h=h,g=g,dh=dh,dg=dg,d2h=d2h,d2g=d2g,variance=sum(w*(exp(h)-1)^2))
}
#' Marginal likelihood for gamma, lognormal or estimated discrete frailty
#' @export
rc_frailty_loglik <- function(theta,unit,exposure,design,event,frailty=c("gamma","lognormal","discrete"),
  phi_fixed=NULL,points=2L,quadrature=30L,information=TRUE) {
  family<-match.arg(frailty)
  if(family=="gamma")return(rc_loglik(theta,unit,exposure,design,event,phi_fixed,information))
  if(!is.matrix(design)||!is.numeric(design)||!ncol(design)||any(!is.finite(design))||any(!is.finite(theta))) .rc_error("frailty","Finite numeric design and parameters required.")
  if(length(unit)!=nrow(design)||!length(unit)||any(!is.finite(unit))||any(unit<1|unit!=floor(unit))||!identical(sort(unique(as.integer(unit))),seq_len(max(unit)))) .rc_error("frailty","Contiguous positive unit indices required.")
  if(length(exposure)!=length(unit)||length(event)!=length(unit)||any(!is.finite(c(exposure,event)))||any(exposure<0)||any(!event %in% c(0,1))) .rc_error("frailty","Aligned nonnegative exposures and binary events required.")
  if(length(quadrature)!=1L||!is.finite(quadrature)||quadrature<2||quadrature>1024||quadrature!=floor(quadrature)) .rc_error("quadrature","quadrature must be an integer from 2 to 1024.")
  spec<-.rc_frailty_spec(family,coherent4cast_control(phi_fixed=phi_fixed,quadrature=quadrature),points)
  p<-ncol(design);if(length(theta)!=p+length(spec$start)) .rc_error("frailty","Parameter length does not match design and frailty family.")
  z<-.rc_frailty_nodes(theta[-seq_len(p)],spec)
  out<-.rc_kernel_guard(cpp_rc_mixture_lik(theta[seq_len(p)],as.integer(unit),exposure,design,event,z$h,z$g,z$dh,z$dg,z$d2h,z$d2g,information))
  if(!information){out$information<-NULL;out$information_individual<-NULL};out
}
.rc_fit_likelihood <- function(theta,dd,control,frailty="gamma",points=2L,information=TRUE) {
  rc_frailty_loglik(theta,dd$unit,dd$exposure,dd$design,dd$event,frailty,control$phi_fixed,points,control$quadrature,information)
}
.rc_posterior_frailty <- function(object,count,hazard,quadrature) {
  family<-if(is.null(object$frailty))"gamma"else object$frailty
  if(family=="gamma") {
    phi<-if(is.null(object$control$phi_fixed))exp(tail(object$coefficients,1))else object$control$phi_fixed
    return(if(phi==0)list(nodes=1,weights=1)else .rc_gamma_nodes(1/phi+count,1/phi+hazard,quadrature))
  }
  p<-ncol(object$design$design);z<-.rc_frailty_nodes(object$coefficients[-seq_len(p)],object$frailty_spec,quadrature)
  nodes<-exp(z$h);w<-.rc_softmax(z$g+count*z$h-exp(z$h)*hazard)
  keep<-w>0;nodes<-nodes[keep];w<-w[keep]
  small<-w<min(1e-12,object$control$propagation_tolerance/10)/length(w)
  if(any(small)) {
    mass<-sum(w[small]);center<-sum(w[small]*nodes[small])/mass
    return(list(nodes=c(nodes[!small],center),weights=c(w[!small],mass),compression_mass=mass,compression_index=sum(!small)+1L))
  }
  list(nodes=nodes,weights=w,compression_mass=0,compression_index=0L)
}
#' Compare fitted frailty families and test a no-frailty boundary
#' @export
rc_frailty_diagnose <- function(object,alternatives=FALSE) {
  if(!inherits(object,"coherent4cast")||!is.null(object$targeting)&&object$targeting!="none") .rc_error("frailty","A likelihood-based coherent4cast fit is required.")
  dd<-object$design;ctl<-object$control;ctl$phi_fixed<-0
  p<-ncol(dd$design);owner<-object$coarsening$owner;lw<-object$coarsening$log_weights
  eval0<-function(b,info=FALSE){z<-.rc_fit_likelihood(b,dd,ctl,information=info);if(!is.null(owner))z<-rc_coarsen(z,owner=owner,log_weights=lw);z}
  sw<-if(is.null(object$case_weights))rep(1,object$n)else object$case_weights
  opt<-stats::optim(object$coefficients[seq_len(p)],function(b)-sum(sw*eval0(b)$subject_loglik),function(b)-colSums(sw*eval0(b)$scores),method="L-BFGS-B",lower=rep(ctl$coefficient_bounds[1],p),upper=rep(ctl$coefficient_bounds[2],p),control=list(maxit=ctl$maxit,factr=10,pgtol=ctl$tolerance))
  variance_scores<-function(b) {
    zz<-.rc_fit_likelihood(b,dd,ctl,information=FALSE)
    score<-.5*((zz$counts-zz$cumulative_hazard)^2-zz$counts)
    if(!is.null(owner)){post<-rc_coarsen(zz,owner=owner,log_weights=lw)$posterior;score<-as.vector(rowsum(score*post,owner))}
    score
  }
  null<-eval0(opt$par,TRUE);u<-variance_scores(opt$par)*sw;S<-null$scores*sw
  H<-apply(sweep(null$information_individual,3,sw,"*"),c(1,2),sum)
  h<-1e-5*pmax(1,abs(opt$par))
  B<-vapply(seq_len(p),function(j){up<-down<-opt$par;up[j]<-up[j]+h[j];down[j]<-down[j]-h[j]
    sum(sw*(variance_scores(up)-variance_scores(down)))/(2*h[j])},numeric(1))
  weighted<-any(sw!=1)
  if(weighted) {
    # beta_hat - beta_0 = H^{-1} sum_i w_i S_i to first order.
    # The derivative of sum_i w_i u_i therefore enters with a plus sign.
    residual<-as.vector(u+S%*%(.rc_inverse(H)%*%B))
    adjustment<-"Weighted estimating-equation nuisance sensitivity: w_i u_i + B H^{-1} w_i S_i."
  } else {
    residual<-as.vector(u-S%*%(.rc_inverse(crossprod(S))%*%crossprod(S,u)))
    adjustment<-"Unweighted sample score projection under the correctly specified null."
  }
  variance<-sum(residual^2);z<-if(variance>0)sum(u)/sqrt(variance)else NA_real_
  eig<-eigen((H+t(H))/2,symmetric=TRUE,only.values=TRUE)$values
  null_gradient<-max(abs(colSums(S)))/sum(sw)
  null_at_boundary<-any(abs(opt$par-ctl$coefficient_bounds[1])<1e-5|abs(opt$par-ctl$coefficient_bounds[2])<1e-5)
  null_regular<-opt$convergence==0L&&min(eig)>1e-8&&!null_at_boundary&&null_gradient<1e-4
  tab<-data.frame(family=if(is.null(object$frailty))"gamma"else object$frailty,logLik=object$loglik,parameters=length(object$coefficients),AIC=-2*object$loglik+2*length(object$coefficients),regular=object$convergence$valid)
  fits<-list()
  if(isTRUE(alternatives))for(fam in setdiff(c("gamma","lognormal","discrete"),tab$family[1])) {
    sc<-object$schema;ct<-object$control;ct$phi_fixed<-NULL
    fit<-suppressWarnings(coherent4cast(object$formula,object$data,unit=sc$unit,start=sc$start,stop=sc$stop,event_type=sc$event_type,terminal=sc$terminal,event_types=sc$event_types,time_unit=sc$time_unit,breaks=object$breaks,observation=sc$observation,control=ct,weights=sw,nonrecurrent=sc$event_types[object$nonrecurrent],coarsening=object$coarsening,frailty=fam,frailty_diagnostics=FALSE))
    fits[[fam]]<-fit;tab<-rbind(tab,data.frame(family=fam,logLik=fit$loglik,parameters=length(fit$coefficients),AIC=-2*fit$loglik+2*length(fit$coefficients),regular=fit$convergence$valid))
  }
  list(null_variance_test=list(statistic=if(is.finite(z))max(0,z)^2 else NA_real_,z=z,p.value=if(!null_regular||!is.finite(z))NA_real_ else if(z<=0)1 else stats::pnorm(z,lower.tail=FALSE),null_converged=opt$convergence==0,
    null_regular=null_regular,null_gradient=null_gradient,null_at_boundary=null_at_boundary,null_information_eigenvalues=eig,
    null_coefficients=opt$par,nuisance_sensitivity=B,weighted_information=H,effective_score_variance=variance,nuisance_adjustment=adjustment,
    calibration=paste("One-sided boundary score: asymptotic 0.5 point mass at zero + 0.5 chi-square(1), under a correctly specified regular no-frailty null and identifiable variance direction.",
      if(weighted)"Weighted calibration additionally requires declared fixed or exogenous case weights, valid weighted estimating equations and a weighted-score central limit theorem; arbitrary outcome-dependent selection weights are not validated."else "Unweighted independent-subject calibration.",
      if(!null_regular)"The null nuisance fit is nonregular or unconverged; the p value is withheld."else "")),
    comparison=tab,alternative_fits=fits,comparison_scope=paste(c("AIC is descriptive for nonnested families; no chi-square likelihood-ratio reference is used for gamma versus lognormal or estimated mixture support.",
      if(any(sw!=1))"With nonunit case weights, -2 times the weighted log likelihood plus 2p is a descriptive weighted criterion; arbitrary survey or importance weights do not confer the standard AIC risk interpretation."),collapse=" "))
}
