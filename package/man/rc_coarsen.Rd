\name{rc_coarsen}
\alias{rc_coarsen}
\title{Integrate Likelihood Contributions over Legal Complete Paths}
\usage{
rc_coarsen(loglik, scores = NULL, information = NULL, owner,
  log_weights)
}
\arguments{
\item{loglik}{One finite log likelihood per enumerated legal path, or a
  result from \code{rc_loglik} in which each unit indexes one path.}
\item{scores}{Numeric matrix with one path per row and one parameter per column.
  Extracted from \code{loglik} when a likelihood result is supplied.}
\item{information}{Optional parameter by parameter by path array of observed
  information matrices. Extracted from a supplied likelihood result.}
\item{owner}{Contiguous positive integer observed-subject indices, one per path.}
\item{log_weights}{Log integration weights, one per legal path. These are
  integration measures, not separately normalized probabilities within subjects.}
}
\description{Aggregates complete-path likelihoods by observed subject using a
stable log-sum-exp calculation. Posterior path weights produce observed scores
and, when supplied, Louis observed information.}
\details{The caller supplies all admissible paths and any numerical quadrature
weights. The routine does not infer missing event orderings from arbitrary
input rows. Observed information equals posterior mean complete-path information
minus posterior covariance of complete-path scores. Validity depends on the
specified observation/coarsening model and the supplied path integration rule.}
\value{A list with the aggregated log likelihood, subject contributions,
scores, posterior path probabilities and optional observed information.}
\seealso{\code{rc_loglik}, \code{rc_filter_paths}}
\examples{
x <- rc_loglik(log(0.3), 1:2, c(1, 1), matrix(1, 2, 1),
               c(0, 1), phi_fixed = 0)
rc_coarsen(x, owner = c(1, 1), log_weights = log(c(0.5, 0.5)))
}
