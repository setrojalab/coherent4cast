\name{confint.coherentSurv}
\alias{confint.coherentSurv}
\alias{confint.coherentSurv_competing}
\alias{confint.coherentSurv_targeted}
\alias{confint.coherentSurv_causal}
\title{Confidence limits with model-specific regularity checks}
\description{Returns confidence limits from the existing likelihood, targeted,
or causal inference procedure while preserving its regularity guards and scope.}
\usage{
\method{confint}{coherentSurv}(object, parm = NULL, level = .95, ...)
\method{confint}{coherentSurv_competing}(object, parm = NULL, level = .95, ...)
\method{confint}{coherentSurv_targeted}(object, parm = NULL, level = .95,
  interval = "auto", draws = 1000L, seed = 1L, ...)
\method{confint}{coherentSurv_causal}(object, parm = NULL, level = .95,
  contrast = c("difference", "ratio", "log_ratio"),
  interval = c("wald", "multiplier", "none"), draws = 2000L, seed = 1L, ...)
}
\arguments{
\item{object}{A fitted coherentSurv likelihood, competing-risk, targeted or
baseline-treatment causal object. The survival wrapper inherits the ordinary
likelihood method.}
\item{parm}{Parameter or target names, or integer indices. Negative indices
omit coordinates; positive and negative indices cannot be mixed. \code{NULL}
selects all fitted coordinates. Selection is applied after the full interval
calculation, preserving its simultaneous family where applicable.}
\item{level}{Confidence level strictly between zero and one; default 0.95.}
\item{interval}{For targeted fits, \code{auto}, \code{none}, \code{multiplier},
\code{logit_gaussian}, or \code{projection}, as in targeted prediction.
For causal fits, pointwise \code{wald}, simultaneous \code{multiplier}, or
\code{none}, as in \code{cs_causal_contrast()}.}
\item{contrast}{For causal fits, treatment-one versus treatment-zero
\code{difference}, \code{ratio}, or \code{log_ratio} of each fitted target.}
\item{draws, seed}{Simulation controls passed to the existing inference
procedure; its random-number-state preservation applies.}
\item{...}{Reserved; unsupported arguments are rejected.}
}
\details{Ordinary likelihood fits receive pointwise sandwich Wald intervals
only when their convergence and regularity flag is true and their numerical
variance is finite and positive. Competing-risk fits additionally require every
cause-specific fit to meet its regularity checks. Unavailable limits are
\code{NA} with an explicit \code{inference_status} attribute.

Targeted fits call the same population-interval procedure as \code{predict()};
individual-history intervals are not provided. Its automatic Gaussian/projection
switch retains the documented limitation: it has no unconditional or
selection-conditional finite-sample guarantee. Request \code{interval = "projection"}
uniformly to invoke the observed-envelope theorem under its stated assumptions.
Explicit Gaussian requests at nonregular boundaries remain unavailable.

Causal fits call \code{cs_causal_contrast()}, preserving its coordinate-specific
arm regularity checks, delta method and simultaneous multiplier option.
These methods do not silently invoke R's ungated default Wald calculation.}
\value{A two-column numeric matrix with selected coordinate names and percentage
column labels. Attributes retain \code{level}, \code{interval},
\code{inference_status} and \code{inference_scope}; targeted outputs also retain
\code{requested_interval}, and causal outputs retain \code{contrast}.
The limits may be asymptotic intervals or explicitly labelled coarsening outer
sets, according to the selected existing procedure.}
\seealso{\code{summary.coherentSurv}, \code{cs_competing},
\code{predict.coherentSurv_targeted}, \code{cs_causal_contrast}}
\examples{
d <- data.frame(time = c(1,2,3,4,2,4), status = c(1,0,1,0,1,0))
fit <- cs_surv(~1, d, time_unit = "years", frailty_diagnostics = FALSE)
confint(fit)
}
