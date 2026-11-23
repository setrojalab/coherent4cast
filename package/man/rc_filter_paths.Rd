\name{rc_filter_paths}
\alias{rc_filter_paths}
\title{Posterior mixture over compatible histories}
\description{Normalize legal-path likelihood weights and retain each gamma posterior component.}
\usage{
rc_filter_paths(loglik, log_weights, counts, cumulative_hazard, 
    dispersion, owner = rep(1L, length(loglik))) 
}
\arguments{
\item{loglik}{One complete-path log likelihood per compatible path.}
\item{log_weights}{Log integration weights of the compatible paths.}
\item{counts}{Complete-path eligible event counts.}
\item{cumulative_hazard}{Integrated complete-path exposure excluding frailty.}
\item{dispersion}{Nonnegative shared-gamma variance (strictly positive for posterior-mixture filtering).}
\item{owner}{Observed-subject index of each compatible path.}
}
\details{Integration weights are multiplied by path likelihoods. The resulting mixture is not replaced by a single gamma distribution. All components for a subject must be compatible with the recorded data.}
\value{A data frame of owner, posterior gamma shape, rate and mixture weight.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
\examples{
rc_filter_paths(log(c(.2,.3)),log(c(.5,.5)),c(1,2),c(.4,.6),.5)
}
