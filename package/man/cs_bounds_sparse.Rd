\name{cs_bounds_sparse}
\alias{cs_bounds_sparse}
\title{Scalable sharp bounds for the marginal recording-odds model}
\usage{
cs_bounds_sparse(p, numerator, denominator = rep(1, length(numerator)),
  kappa, recording_fraction, counts = NULL, alpha = .05,
  observation_map = seq_along(numerator), tolerance = 1e-10,
  max_iterations = 1000L)
}
\arguments{
\item{p}{Recorded-category probabilities, summing to one.}
\item{numerator, denominator}{One finite numerator and nonnegative denominator value per complete-history atom. A positive recorded lower bound on the conditioning denominator is required for the fractional optimization.}
\item{kappa}{Marginal recording-odds sensitivity bound, at least one.}
\item{recording_fraction}{Known overall recording fraction strictly between zero and one.}
\item{counts}{Optional nonnegative integer recorded-category counts whose proportions equal \code{p}.}
\item{alpha}{Error probability for the simultaneous recorded-law confidence region.}
\item{observation_map}{1-based category index for each complete-history atom. Every category must have at least one atom.}
\item{tolerance}{Stopping tolerance for the fractional endpoint certificate, scaled by one plus the absolute endpoint value.}
\item{max_iterations}{Maximum fractional optimization iterations.}
}
\description{Solves the same finite-support marginal recording-odds and coarsening-fibre problem as \code{cs_bounds}, through an exact structural reduction rather than a dense simplex tableau. Supports large complete-history supports without quadratic storage.}
\details{
Let \eqn{a_j} and \eqn{b_j} be recorded and unrecorded joint masses. With recording fraction \eqn{\rho}, constrain their totals to \eqn{\rho} and \eqn{1-\rho}, and \eqn{l a_j\le b_j\le u a_j}, where \eqn{l=(1-\rho)/(\rho\kappa)} and \eqn{u=(1-\rho)\kappa/\rho}. Each observation fibre constrains its total recorded mass. The target is the numerator divided by denominator under \eqn{q_j=a_j+b_j}.

The substitution \eqn{a=x+y}, \eqn{b=lx+uy} fixes the total mass of \eqn{y}. For each linearized fractional objective, both components choose the same maximizing atom within each fibre. Sorting fibre maxima then solves the bounded allocation exactly. Fractional iterations terminate with feasible complete laws and dual endpoint brackets. This is equivalent to the Charnes--Cooper LP for these constraints; it is not a solver for arbitrary sparse linear programs. Per-iteration work is \eqn{O(J+M\log M)} and memory is \eqn{O(J+M)} for \eqn{J} atoms and \eqn{M} fibres.

With \code{counts}, Bonferroni-adjusted Clopper--Pearson intervals form a simultaneous recorded-category confidence box intersected with the fixed total mass. The same solver projects the full box. If the outer conditioning-denominator lower bound is zero, a bounded support envelope is returned and explicitly labelled. The recording fraction, complete support, and sensitivity bound remain assumptions. This projection does not identify unrecorded outcomes without those restrictions.
}
\value{An \code{cs_bounds} object. The \code{plugin} and, where available, \code{outer} lists contain outward endpoint certificates \code{lower} and \code{upper}, feasible endpoint values and laws, maximum fractional \code{gap}, primal \code{residual}, and lower/upper certificate lists. Each certificate includes recorded/unrecorded masses, linearization point, objective sign, linear primal and dual values, mass multipliers, fibre dual terms, iteration count, conditioning denominator, and a fractional endpoint bracket. Floating-point certificates are checked to the reported tolerance; they are not symbolic or directed-rounding proofs.}
\references{
Charnes, A. and Cooper, W. W. (1962). Programming with linear fractional functionals. \emph{Naval Research Logistics Quarterly}, 9, 181--186. doi:10.1002/nav.3800090303.
Dinkelbach, W. (1967). On nonlinear fractional programming. \emph{Management Science}, 13, 492--498. doi:10.1287/mnsc.13.7.492.
}
\seealso{\code{\link{cs_bounds}}}
\examples{
# Exact two-atom trace: conditional bounds are 3/11 and 5/13.
z <- cs_bounds_sparse(c(.5, .5), c(0, 1), c(1, 2),
                      kappa = 2, recording_fraction = .5)
stopifnot(abs(z$lower - 3/11) < 1e-9, abs(z$upper - 5/13) < 1e-9)
}
