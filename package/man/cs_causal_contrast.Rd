\name{cs_causal_contrast}
\alias{cs_causal_contrast}
\alias{predict.coherentSurv_causal}
\alias{tidy.coherentSurv_causal}
\alias{coef.coherentSurv_causal}
\alias{vcov.coherentSurv_causal}
\alias{logLik.coherentSurv_causal}
\alias{glance.coherentSurv_causal}
\alias{summary.coherentSurv_causal}
\alias{print.coherentSurv_causal}
\title{Influence-Function Contrasts for Baseline Treatment Effects}
\usage{
cs_causal_contrast(object, contrast = c("difference", "ratio", "log_ratio"),
  level = .95, interval = c("wald", "multiplier", "none"),
  draws = 2000L, seed = 1L)
\method{predict}{coherentSurv_causal}(object, contrast = "difference", ...)
\method{tidy}{coherentSurv_causal}(x, ...)
\method{coef}{coherentSurv_causal}(object, ...)
\method{vcov}{coherentSurv_causal}(object, ...)
\method{logLik}{coherentSurv_causal}(object, ...)
\method{glance}{coherentSurv_causal}(x, ...)
\method{summary}{coherentSurv_causal}(object, ...)
\method{print}{coherentSurv_causal}(x, ...)
}
\arguments{
\item{object,x}{A fitted \code{coherentSurv_causal} object.}
\item{contrast}{Treatment 1 minus 0, treatment 1 divided by 0, or log of that ratio. Ratios require positive arm means.}
\item{level}{Confidence level strictly between zero and one.}
\item{interval}{Pointwise normal approximation, simultaneous Gaussian multiplier approximation over the reported regular coordinates, or no intervals.}
\item{draws}{Number of Gaussian multiplier draws, at least 100. A small covariance factor avoids storing a draws-by-subject matrix.}
\item{seed}{Resampling seed, with restoration of the caller's random state.}
\item{...}{Arguments passed to \code{cs_causal_contrast}, or unused method arguments.}
}
\description{Uses both arm-specific gradients on the same subjects to retain their covariance. Difference gradients are \eqn{D^1-D^0}; log-ratio gradients are \eqn{D^1/\psi(1)-D^0/\psi(0)}. Here \eqn{D^w} is the arm-\eqn{w} gradient and \eqn{\psi(w)} its mean outcome. Ratio intervals are constructed on the log scale and exponentiated.}
\details{Intervals are conditional on the causal identification and asymptotic regularity conditions of \code{cs_causal}. No exact finite-sample coverage is claimed. Nonconverged, zero-event or zero-variance coordinates receive missing interval limits and a status flag. A simultaneous interval covers the retained regular coordinates only. Prognostic recording-sensitivity bounds are not substituted for treatment-effect confidence intervals. The treatment-effect estimator is standardized over the cohort's empirical baseline distribution; its influence function includes baseline-population variation.}
\value{For contrasts, a data frame with target, contrast, both arm means, estimate, standard error, limits and inference status. Attributes contain the contrast EIF, covariance and critical value. Ratio standard errors are on the ratio scale, while their interval calculation uses the log scale. Other methods return the named model summaries. \code{logLik} raises an explicit error because cross-fitting does not define one maximized parametric likelihood.}
\seealso{\code{\link{cs_causal}}}
