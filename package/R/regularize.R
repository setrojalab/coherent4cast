.cs_penalty_score <- function(theta, design, inverse) {
  X<-design$design; ids<-design$unit; p<-ncol(X); n<-max(ids)
  total<-function(v) as.numeric(rowsum(v,ids,reorder=FALSE)[match(seq_len(n),unique(ids)),,drop=FALSE])
  v<-design$exposure*exp(as.vector(X%*%theta[seq_len(p)]));a<-exp(-theta[p+1L])
  A<-total(v);d<-total(design$event)
  A1<-vapply(seq_len(p),function(j)total(v*X[,j]),numeric(n))
  h<-a+A;w<-(a+d)/h;t<-w/h;g<- -a*(A-d)/h^2
  K<-inverse[seq_len(p),seq_len(p),drop=FALSE];k<-inverse[seq_len(p),p+1L];c0<-inverse[p+1L,p+1L]
  s<-rowSums((X%*%K)*X);C<-total(v*s);AK<-A1%*%K;B<-rowSums(AK*A1);E<-as.vector(A1%*%k)
  ga<- -a/h^2+2*a*(A-d)/h^3;ha<-a*(A-d)*(A-a)/h^3
  cf<- -t*C+2*t*B/h+2*ga*E+c0*ha
  rt<-w[ids]*s-2*t[ids]*rowSums(X*AK[ids,,drop=FALSE])+2*g[ids]*as.vector(X%*%k)
  adj<-numeric(p+1L);adj[seq_len(p)]<-.5*(as.vector(crossprod(A1,cf))+as.vector(crossprod(X,v*rt)))
  dig<-tri<-tet<-numeric(n)
  if(max(d)>0)for(j in 0:(max(d)-1L)) {yes<-d>j;dig[yes]<-dig[yes]+1/(a+j);tri[yes]<-tri[yes]-1/(a+j)^2;tet[yes]<-tet[yes]+2/(a+j)^3}
  ua<-dig+log(a)+1-log(h)-(a+d)/h
  uaa<-tri+1/a-1/h+(d-A)/h^2
  uaaa<-tet-1/a^2+1/h^2-2*(d-A)/h^3
  wa<-(A-d)/h^2;ta<-(A-a-2*d)/h^3;gga<-(A-d)*(a-A)/h^3;hha<- -ua-3*a*uaa-a*a*uaaa
  adj[p+1L]<- -.5*a*sum(wa*C-ta*B+2*gga*E+c0*hha)
  adj
}

#' Observed-information log-determinant stabilization
#' @export
cs_regularize <- function(object,maxit=120L,tolerance=1e-8) {
  if(!inherits(object,"coherentSurv")||(!is.null(object$frailty)&&object$frailty!="gamma")||object$schema$observation!="exact"||!is.null(object$control$phi_fixed)||
     !is.null(object$audit$call$weights))
    .cs_error("regularization_scope","Requires an unweighted exact-history coherentSurv fit with estimated positive frailty dispersion.")
  if(length(maxit)!=1L||!is.finite(maxit)||maxit<1||maxit!=floor(maxit)||length(tolerance)!=1L||!is.finite(tolerance)||tolerance<=0)
    .cs_error("regularization_control","Use positive finite controls and an integer maxit.")
  dd<-object$design;p<-ncol(dd$design);n<-object$n
  evaluate<-function(theta) {
    z<-cs_loglik(theta,dd$unit,dd$exposure,dd$design,dd$event)
    I<-z$information;ee<-eigen(I,symmetric=TRUE,only.values=TRUE)$values
    if(min(ee)<=1e-10)return(list(value=1e20,gradient=rep(0,length(theta)),valid=FALSE))
    inv<-solve(I);adjust<-.cs_penalty_score(theta,dd,inv)
    list(value=z$objective-.5*sum(log(ee))/n,gradient=z$gradient-adjust/n,valid=TRUE,raw=z)
  }
  if(!evaluate(object$coefficients)$valid) .cs_error("positive_information","Starting fit is outside the positive-information domain; no information ridge is silently added.")
  lo<-c(rep(object$control$coefficient_bounds[1],p),object$control$log_phi_bounds[1]);hi<-c(rep(object$control$coefficient_bounds[2],p),object$control$log_phi_bounds[2])
  opt<-stats::optim(object$coefficients,function(t)evaluate(t)$value,function(t)evaluate(t)$gradient,method="L-BFGS-B",lower=lo,upper=hi,
      control=list(maxit=maxit,factr=max(1,min(100,tolerance/.Machine$double.eps)),pgtol=tolerance))
  z<-evaluate(opt$par)
  if(!z$valid) .cs_error("positive_information","Penalized optimizer left the positive-information domain.")
  step<-.Machine$double.eps^(1/3)*pmax(1,abs(opt$par));H<-matrix(0,p+1L,p+1L)
  for(j in seq_len(p+1L)) {up<-dn<-opt$par;up[j]<-up[j]+step[j];dn[j]<-dn[j]-step[j];u<-evaluate(up);v<-evaluate(dn)
    if(!u$valid||!v$valid) .cs_error("penalized_curvature","Finite-difference curvature crosses the information-domain boundary.")
    H[,j]<-(u$gradient-v$gradient)/(2*step[j])*n}
  H<-(H+t(H))/2;ev<-eigen(H,symmetric=TRUE,only.values=TRUE)$values
  if(min(ev)<=1e-8) .cs_error("penalized_curvature","Penalized curvature is not positive definite.")
  inv<-solve(H);centered<-sweep(z$raw$scores,2,colMeans(z$raw$scores));V<-inv%*%crossprod(centered)%*%inv
  at_bound<-any(abs(opt$par-lo)<1e-5|abs(opt$par-hi)<1e-5)
  valid<-opt$convergence==0L&&!at_bound&&max(abs(z$gradient))<1e-4
  nm<-names(object$coefficients);names(opt$par)<-nm;dimnames(V)<-list(nm,nm)
  object$coefficients<-opt$par;object$covariance<-V;object$information<-H;object$scores<-z$raw$scores
  if(!is.null(object$linear_coefficients))object$linear_coefficients[dd$active]<-opt$par[seq_len(p)]
  object$phi<-exp(opt$par[p+1L]);object$loglik<-sum(z$raw$subject_loglik)
  object$convergence<-list(code=if(valid)"converged" else "nonregular_fit",optimizer=opt$convergence,message=opt$message,
    gradient=max(abs(z$gradient)),information_eigenvalues=ev,at_boundary=at_bound,valid=valid)
  object$audit$regularization<-list(method="Observed-information half log-determinant",penalized_objective=z$value,
      derivative="Analytic log-determinant trace score; central-difference penalized curvature",scope="Local positive-information-domain regularization; no general Firth bias cancellation or Bayesian posterior claim")
  if(!valid)warning("Penalized fit did not meet regular convergence criteria; inspect $convergence.",call.=FALSE)
  object
}
