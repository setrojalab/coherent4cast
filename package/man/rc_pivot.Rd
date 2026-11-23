\name{rc_pivot}
\alias{rc_pivot}
\title{Exact pivot for the known-dispersion terminal-time submodel}
\usage{
rc_pivot(times, dispersion, horizons, alpha = 0.05)
}
\arguments{
\item{times}{Positive independent uncensored terminal-event times.}
\item{dispersion}{Known gamma-frailty variance, nonnegative. Zero is the exponential submodel.}
\item{horizons}{Positive prediction horizons in the time units of times.}
\item{alpha}{Simultaneous confidence error probability.}
}
\description{Inverts the exact gamma-distributed transformed-time pivot and maps its rate interval to terminal-risk curves.}
\details{For positive dispersion phi and baseline rate lambda, sum(log(1+phi*lambda*T))/phi has a Gamma(n,1) distribution. The inversion assumes a constant baseline and known dispersion. It is not an exact interval for estimated dispersion, censored observations, or general recurrent histories. Monotonicity transfers the common rate interval to simultaneous risk limits at all positive horizons.}
\value{Rate confidence interval, horizon-specific risk limits, and inferential scope.}
\examples{
rc_pivot(c(1, 2, 3, 4), dispersion = .3, horizons = c(.5, 1, 2))
}
