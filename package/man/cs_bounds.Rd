\name{cs_bounds}
\alias{cs_bounds}
\title{Finite-support recording-sensitivity bounds and confidence projection}
\usage{
cs_bounds(p, numerator, denominator = rep(1, length(numerator)),
  kappa, recording_fraction, counts = NULL, alpha = 0.05,
  observation_map = seq_along(numerator))
}
\arguments{
\item{p}{Probability vector of recorded categories conditional on recording.}
\item{numerator}{Reward numerator on the declared complete-path support.}
\item{denominator}{Nonnegative conditioning reward on the complete support. The default is an unconditional mean.}
\item{kappa}{Marginal recording-odds sensitivity strength, at least one.}
\item{recording_fraction}{Known complete-path recording fraction strictly between zero and one.}
\item{counts}{Optional recorded-category counts with proportions equal to p. Supplies simultaneous Clopper-Pearson confidence constraints.}
\item{alpha}{Overall confidence error probability.}
\item{observation_map}{Integer category index for each complete path. Declares a deterministic partition into observed coarsening fibres.}
}
\description{Computes finite-support recording-sensitivity bounds. Supports with at most 64 atoms use the dense Charnes--Cooper simplex implementation; larger supports automatically use the exact structural reduction in \code{\link{cs_bounds_sparse}}. The statistical constraints and estimand are identical.}
\details{The constraints are u_j/r_j between (1-rho)/(rho*kappa) and (1-rho)*kappa/rho, total recorded mass rho, and total missing mass 1-rho. Recorded category masses are fixed at rho*p, or lie in simultaneous confidence intervals. The observation map sums recorded complete-path masses within each declared fibre. The conditioning denominator must have a strictly positive recorded lower bound for LP inversion. If this condition fails only in the confidence projection, a conservative support envelope is returned and explicitly labelled. Bounds are sharp for this declared finite model, not for unknown support or unknown recording fraction. The outer projection uses Bonferroni Clopper-Pearson intervals, not EIF scores.}
\value{Lower and upper population bounds, optimizer laws, primal residual and primal-dual gap. With counts, also returns an outer confidence projection and its recorded-law confidence region.}
\seealso{\code{\link{cs_bounds_sparse}} for scalable optimization and reproducible primal/dual certificates.}
\examples{
cs_bounds(c(.5, 1/3, 1/6), c(0, 0, 1), c(0, 1, 1),
          kappa = 2, recording_fraction = .7)
cs_bounds(c(.5, .5), c(0, 1, 0), kappa = 1,
          recording_fraction = .7, observation_map = c(1, 2, 2))
}
