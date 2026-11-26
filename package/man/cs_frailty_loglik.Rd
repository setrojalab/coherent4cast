\name{cs_frailty_loglik}
\alias{cs_frailty_loglik}
\title{Alternative shared-frailty likelihoods}
\description{Evaluate gamma, mean-one lognormal or estimated mean-one discrete frailty likelihoods.}
\usage{
cs_frailty_loglik(theta, unit, exposure, design, event, frailty = c("gamma", 
    "lognormal", "discrete"), phi_fixed = NULL, points = 2L, 
    quadrature = 30L, information = TRUE) 
}
\arguments{
\item{theta}{Regression coefficients followed by the specified family's unconstrained frailty parameters.}
\item{unit}{Contiguous positive integer subject indices, one per expanded exposure row; repeated rows belong to the same independent subject.}
\item{exposure}{Eligible exposure duration for each expanded marked-event row.}
\item{design}{Numeric design matrix, one row per expanded exposure.}
\item{event}{Binary event indicator per expanded exposure row.}
\item{frailty}{Mean-one gamma (analytic integration), lognormal (Gaussian quadrature), or estimated discrete mass-point frailty.}
\item{phi_fixed}{A fixed nonnegative gamma or lognormal frailty variance. NULL estimates log dispersion; zero selects the no-frailty boundary. Estimated discrete support requires NULL.}
\item{points}{Number of discrete frailty mass points.}
\item{quadrature}{Fixed Gaussian quadrature node count for this lognormal likelihood evaluation, from 2 to 1024. Gamma likelihood integration is analytic; estimated discrete support uses points instead. This low-level call does not adaptively refit.}
\item{information}{Whether to return observed information and each independent unit's information matrix.}
}
\details{theta contains design coefficients followed by log variance for estimated lognormal frailty, or weight logits and log locations for discrete frailty (the last logit and location are fixed at zero before mean-one normalization). Lognormal Gaussian-quadrature nodes vary analytically with dispersion. Discrete weights and support locations are differentiated through normalization. Louis' identity combines complete-component scores and Hessians. Information is exact for the stated quadrature likelihood, not a certificate for quadrature integration error. Increase quadrature to check stability; mixtures with merging support points are nonregular.}
\value{Subject log likelihoods, scores, observed information, posterior component weights, counts and cumulative hazards.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
\examples{
cs_frailty_loglik(c(log(.4),log(.7)),1,2,matrix(1,1,1),1,
                  frailty='lognormal',quadrature=40)
}
