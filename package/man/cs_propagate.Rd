\name{cs_propagate}
\alias{cs_propagate}
\title{Ordered Uniformization of a Finite State Law}
\usage{cs_propagate(probability, generator, duration, reward = NULL,
  tolerance = 1e-12, max_terms = 10000L)}
\arguments{
\item{probability}{Initial row probability vector.}
\item{generator}{A generator matrix, or list in temporal order.}
\item{duration}{Nonnegative duration per generator.}
\item{reward}{Optional state by reward matrix or matching list; a vector specifies one reward.}
\item{tolerance}{Positive series truncation tolerance less than one.}
\item{max_terms}{Maximum Poisson terms in each scaled calculation.}}
\description{Propagates probability and integrated state rewards through ordered piecewise-constant generators using scaling, Poisson uniformization and squaring.}
\details{This implements the conditional propagation in Theorem 3 and equation (9). The core holds the latent frailty fixed: average over common frailty components only after all ordered segments. Capped state counts may encode rates but do not cap recurrent reward accumulation. Transition rows are normalized after truncation and squaring to stabilize probability conservation; the truncation bound includes the normalization contribution. Returned truncation bounds cover omitted Poisson terms; \code{roundoff_bound} is reported separately. Dense matrix operations have cubic state-space cost. At most 512 states are supported by this implementation. No exponential-to-polynomial state compression is claimed.}
\value{A list containing \code{probability}, integrated \code{reward}, \code{error_bound}, \code{reward_error_bound}, separate \code{roundoff_bound}, \code{terms} and \code{squarings}.}
\examples{
Q <- matrix(c(-0.3, 0.3, 0, 0), 2, 2, byrow = TRUE)
cs_propagate(c(1, 0), Q, 2, reward = c(1, 0))
cs_propagate(c(1, 0), list(Q, 2 * Q), c(1, 1))
}
