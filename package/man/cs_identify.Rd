\name{cs_identify}
\alias{cs_identify}
\title{Observable-contraction diagnostic for frailty identification}
\usage{
cs_identify(iterate, survival, contrast, iterations = 20L)
}
\arguments{
\item{iterate}{Observable contraction mapping a survival value L(A) to L(A/contrast).}
\item{survival}{Initial survival value strictly between zero and one.}
\item{contrast}{Known multiplicative contrast greater than one; reverse the groups for a contrast below one.}
\item{iterations}{Positive integer number of contractions.}
}
\description{Calculates contrast^m times one minus the m-fold contracted survival, making the population identification limit explicit.}
\details{The identification argument requires a common mean-one positive frailty independent of covariates, support near time zero, and a nonunit multiplicative contrast. This numerical diagnostic is not an estimator with sampling inference. Large iterations amplify floating-point cancellation; the rounding_limit attribute identifies contraction values rounded to one.}
\value{A data frame of iteration, contracted survival and recovered baseline. Attributes record model scope and a floating-point rounding limit.}
\examples{
phi <- .7
contract <- function(s) (1 + (s^(-phi) - 1)/2)^(-1/phi)
cs_identify(contract, (1 + phi*.8)^(-1/phi), 2, 10)
}
