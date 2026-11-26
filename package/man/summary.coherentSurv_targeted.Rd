\name{summary.coherentSurv_targeted}
\alias{summary.coherentSurv_targeted}
\alias{print.coherentSurv_targeted}
\alias{tidy.coherentSurv_targeted}
\alias{glance.coherentSurv_targeted}
\alias{coef.coherentSurv_targeted}
\alias{logLik.coherentSurv_targeted}
\title{Inspect targeted estimates and numerical diagnostics}
\usage{
\method{summary}{coherentSurv_targeted}(object, ...)
\method{print}{coherentSurv_targeted}(x, ...)
\method{tidy}{coherentSurv_targeted}(x, interval = "auto", level = .95,
  draws = 1000L, seed = 1L, ...)
\method{glance}{coherentSurv_targeted}(x, ...)
\method{coef}{coherentSurv_targeted}(object, ...)
\method{logLik}{coherentSurv_targeted}(object, ...)
}
\arguments{
\item{object, x}{A targeted coherentSurv fit.}
\item{interval, level, draws, seed}{Population interval controls passed to the targeted prediction method.}
\item{...}{Reserved additional arguments.}
}
\description{Reports target estimates, EIF diagnostics, regularity status and the fitted observation experiment. coef returns the population target vector, not event-intensity regression coefficients.}
\value{summary returns a diagnostic list, tidy a target table, glance a one-row diagnostic table, and coef the named target vector. print returns its argument invisibly. logLik raises a structured scope error because a targeted population functional does not supply a maximized likelihood for AIC.}
\seealso{\code{targeting_spec}, \code{predict.coherentSurv_targeted}}
