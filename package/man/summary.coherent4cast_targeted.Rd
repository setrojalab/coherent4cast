\name{summary.coherent4cast_targeted}
\alias{summary.coherent4cast_targeted}
\alias{print.coherent4cast_targeted}
\alias{tidy.coherent4cast_targeted}
\alias{glance.coherent4cast_targeted}
\alias{coef.coherent4cast_targeted}
\alias{logLik.coherent4cast_targeted}
\title{Inspect targeted estimates and numerical diagnostics}
\usage{
\method{summary}{coherent4cast_targeted}(object, ...)
\method{print}{coherent4cast_targeted}(x, ...)
\method{tidy}{coherent4cast_targeted}(x, interval = "auto", level = .95,
  draws = 1000L, seed = 1L, ...)
\method{glance}{coherent4cast_targeted}(x, ...)
\method{coef}{coherent4cast_targeted}(object, ...)
\method{logLik}{coherent4cast_targeted}(object, ...)
}
\arguments{
\item{object, x}{A targeted coherent4cast fit.}
\item{interval, level, draws, seed}{Population interval controls passed to the targeted prediction method.}
\item{...}{Reserved additional arguments.}
}
\description{Reports target estimates, EIF diagnostics, regularity status and the fitted observation experiment. coef returns the population target vector, not event-intensity regression coefficients.}
\value{summary returns a diagnostic list, tidy a target table, glance a one-row diagnostic table, and coef the named target vector. print returns its argument invisibly. logLik raises a structured scope error because a targeted population functional does not supply a maximized likelihood for AIC.}
\seealso{\code{targeting_spec}, \code{predict.coherent4cast_targeted}}
