\name{cs_cornish_fisher}
\alias{cs_cornish_fisher}
\title{Cornish-Fisher correction for a standardized sample mean}
\usage{
cs_cornish_fisher(probability, n, skewness = 2,
  excess_kurtosis = 6, order = 2L)
}
\arguments{
\item{probability}{Probabilities strictly inside the unit interval.}
\item{n}{Positive sample size.}
\item{skewness}{Standardized third cumulant of one summand.}
\item{excess_kurtosis}{Standardized fourth cumulant of one summand.}
\item{order}{Zero for normal, one for skewness correction, or two including fourth-cumulant and squared-skewness terms.}
}
\description{Evaluates the normal, first-correction, or second-correction standardized-mean quantile formula used in the Version 5 exact-gamma diagnostic.}
\details{The second correction is kurtosis*(z^3-3*z)/(24*n) minus skewness^2*(2*z^3-5*z)/(36*n). These formulas require the corresponding cumulants and expansion regularity. A correction is not an exact confidence interval or a general third-order guarantee for a sandwich estimator.}
\value{A numeric vector of approximate standardized quantiles.}
\examples{
cs_cornish_fisher(c(.025, .975), n = 80)
}
