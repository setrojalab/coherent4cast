\name{rc_propagate_sparse}
\alias{rc_propagate_sparse}
\title{Sparse ordered uniformization with integrated rewards}
\usage{
rc_propagate_sparse(probability, generator, duration, reward = NULL,
  tolerance = 1e-12, max_terms = 10000L, max_chunks = 100000L)
}
\arguments{
\item{probability}{Finite nonnegative initial row probability vector, with positive mass at most one.}
\item{generator}{An \code{rc_generator_csr} object, or an ordered list of such objects with identical state dimensions.}
\item{duration}{Nonnegative duration for each generator.}
\item{reward}{Optional state-by-reward matrix, a state vector for one reward, or an ordered list with one such object per segment. Signed rewards are allowed.}
\item{tolerance}{Positive total probability-series error budget, less than one.}
\item{max_terms}{Maximum Poisson-series terms per chunk.}
\item{max_chunks}{Maximum chunks per segment. Exceeding a budget returns an error instead of silently approximating.}
}
\description{Computes a row distribution propagated through ordered piecewise-constant continuous-time generators and expected integrated rewards, using only CSR matrix-vector products. No dense state-by-state transition matrix is constructed.}
\details{
For a generator \eqn{Q}, set \eqn{\nu=\max_i(-Q_{ii})} and \eqn{P=I+Q/\nu}. The distribution is the Poisson mixture of \eqn{pP^k}. Reward terms use \eqn{\sum_{l=0}^{k-1}pP^lr/\nu}. Durations are split until each Poisson mean is at most 16. Terms are retained until geometric bounds on the Poisson tail meet the allocated error budget. The terminal distribution is mass-normalized; its L1 truncation certificate includes the factor of two required by this normalization. Rewards are not renormalized.

The reported probability and reward bounds concern truncation in exact arithmetic. \code{roundoff_bound} is a separate operation-count diagnostic, not a directed-rounding interval certificate. Work is \eqn{O(K(\mathrm{nnz}+n r))} for total retained terms \eqn{K}, \eqn{n} states, and \eqn{r} reward columns. Working memory is \eqn{O(\mathrm{nnz}+n r)} including inputs. Runtime also depends on rate times duration; the algorithm is not constant-cost at arbitrarily high rates.

A shared latent frailty must be held fixed over all segments of an individual propagation. Average over frailty after composing the entire ordered path. \code{rc_propagate} remains available as an independent dense small-state implementation.
}
\value{A list with \code{probability}, integrated \code{reward}, \code{error_bound}, \code{reward_error_bound}, \code{roundoff_bound}, total \code{terms}, total \code{chunks}, estimated native temporary \code{workspace_bytes}, and the method label. The workspace estimate excludes input arrays and R conversion overhead.}
\seealso{\code{\link{rc_generator_csr}}, \code{\link{rc_propagate}}}
\examples{
q1 <- rc_generator_csr(2, 1, 2, .2)
q2 <- rc_generator_csr(2, 1, 2, .5)
z <- rc_propagate_sparse(c(1, 0), list(q1, q2), c(1, 2), c(1, 0))
stopifnot(abs(z$probability[1] - exp(-1.2)) < 1e-10)
}
