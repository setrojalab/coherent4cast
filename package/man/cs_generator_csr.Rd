\name{cs_generator_csr}
\alias{cs_generator_csr}
\title{Construct a compressed sparse row generator}
\usage{
cs_generator_csr(n, from = integer(), to = integer(), rate = numeric())
}
\arguments{
\item{n}{Number of states, including absorbing states.}
\item{from, to}{Equal-length integer vectors of 1-based source and destination state indices. Only off-diagonal transitions are supplied.}
\item{rate}{Finite, nonnegative transition rates. Duplicate transitions add their rates. Zero rates are discarded.}
}
\description{Creates a generator without constructing a dense state-by-state matrix. The diagonal is defined implicitly as the negative sum of the off-diagonal rates in each row. A state with no outgoing rates is absorbing.}
\value{An object of class \code{cs_generator_csr} with \code{n}, zero-based CSR \code{offsets} and \code{columns}, and off-diagonal \code{rates}. Storage is linear in the number of states plus nonzero transitions.}
\details{The constructor sorts edges by source and destination. It does not infer transitions from clinical associations. A recurrence that leaves the retained state unchanged belongs in a reward rate, not on the generator diagonal.}
\seealso{\code{\link{cs_propagate_sparse}}}
\examples{
q <- cs_generator_csr(3, c(1, 1, 2), c(2, 3, 3), c(.3, .2, .4))
cs_propagate_sparse(c(1, 0, 0), q, 2, reward = c(1, 1, 0))
}
