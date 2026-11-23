library(recurcast)
out<-Sys.getenv("RECURCAST_BENCHMARK_OUT",unset="frailty_boundary_results")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(202611230)
B<-5000L;result<-list()
for(n in c(250L,1000L))for(phi in c(0,.3)) {
  reject<-logical(B);undefined<-0L
  for(b in seq_len(B)) {
    z<-if(phi==0)rep(1,n)else rgamma(n,shape=1/phi,rate=1/phi)
    d<-rpois(n,.6*z);A<-mean(d);S<-d-A;u<-.5*((d-A)^2-d)
    residual<-if(sum(S^2)>0)u-S*sum(S*u)/sum(S^2)else u
    variance<-sum(residual^2)
    if(variance<=0||!is.finite(variance)){undefined<-undefined+1L;reject[b]<-FALSE}else reject[b]<-sum(u)/sqrt(variance)>qnorm(.95)
  }
  p<-mean(reject);result[[length(result)+1L]]<-data.frame(n=n,frailty_variance=phi,B=B,rejection=p,MCSE=sqrt(p*(1-p)/B),undefined=undefined)
}
write.csv(do.call(rbind,result),file.path(out,"frailty_boundary_score.csv"),row.names=FALSE)
cat("Separately registered boundary-score experiment completed.\n")
