\name{cs_conformal}
\alias{cs_conformal}
\title{Weighted conformal rank inversion with candidate mass}
\usage{
cs_conformal(scores, candidate_scores, weights = rep(1, length(scores)),
  candidate_weights = rep(1, length(candidate_scores)), weight_lower = NULL,
  weight_upper = NULL, candidate_weight_upper = NULL, alpha = 0.05,
  band_error = 0, weight_status = c("known", "estimated", "certified"))
}
\arguments{
\item{scores}{Calibration nonconformity scores; larger means less conforming.}
\item{candidate_scores}{Candidate-specific nonconformity scores.}
\item{weights}{Known or estimated nonnegative calibration density ratios.}
\item{candidate_weights}{Positive candidate density ratios; these are included in every rank.}
\item{weight_lower,weight_upper}{Simultaneous calibration weight limits for certified inversion. Upper limits may be infinite.}
\item{candidate_weight_upper}{Candidate upper weight limits for certified inversion.}
\item{alpha}{Total error allocation.}
\item{band_error}{Error allocated to the simultaneous nuisance-weight band, smaller than alpha.}
\item{weight_status}{Known, estimated, or certified. Estimated weights alone do not receive an exact guarantee.}
}
\description{Computes weighted tail ranks, including score ties and candidate mass. Certified inversion maximizes the rank over the supplied rectangular weight intervals.}
\details{The known-weight rank is (candidate weight + calibration weights with score at least the candidate score) divided by all weights including the candidate. For certified weights the numerator uses upper tail and candidate weights and the complementary mass uses lower weights. Candidate inclusion uses p greater than alpha minus band_error. Finite-sample marginal coverage requires weighted exchangeability and a valid independently estimated simultaneous weight band when used. Arbitrary estimated or clipped weights do not inherit this statement. A supplied finite candidate grid is only a grid evaluation; no interpolation guarantee is asserted.}
\value{A list containing candidate scores, p-values, inclusion indicators, error allocations, and scope statements.}
\examples{
cs_conformal(c(.1, .5, .9), c(.2, .8), weights = c(1, 2, 1),
             candidate_weights = c(1, 1))
}
