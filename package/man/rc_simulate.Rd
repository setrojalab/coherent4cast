\name{rc_simulate}
\alias{rc_simulate}
\title{Synthetic multitype gamma frailty histories}
\description{Generate complete recurrent and terminal paths from a capped-history model.}
\usage{
rc_simulate(n = 100L, horizon = 3, rates = c(0.6, 0.25, 
    0.12), beta = c(-0.2, -0.1, 0.15), history = c(0.2, 
    0.2, 0.25), dispersion = 0.45, censor_rate = 0.1, 
    seed = 20261122L, time_unit = "years", history_cap = 2L) 
}
\arguments{
\item{n}{Number of independent synthetic units.}
\item{horizon}{One or more positive forecast durations; a single duration for validation or simulation.}
\item{rates}{Positive baseline rates for each recurrent type followed by the terminal type.}
\item{beta}{Covariate coefficients matching rates.}
\item{history}{Coefficients on the sum of capped prior recurrent counts, matching rates.}
\item{dispersion}{Nonnegative shared-gamma variance (strictly positive for posterior-mixture filtering).}
\item{censor_rate}{Nonnegative independent censoring rate.}
\item{seed}{Integer seed recorded in provenance or used by a simulation.}
\item{time_unit}{Explicit nonempty unit of follow-up time, such as years or operating hours.}
\item{history_cap}{Count-feature cap used for estimation and finite-state forecasting.}
}
\details{The last rate is terminal. Simulations use a locally restored L'Ecuyer-CMRG stream and do not leave the global random state changed. Synthetic field labels do not imply clinical or engineering validation.}
\value{A data frame with id, start, stop, event_type, terminal, x and observed, plus true parameters and seed metadata.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
\examples{
d <- rc_simulate(n=10,seed=17,time_unit='years')
head(d)
}
