\name{summary.recurcast_targeted}
\alias{summary.recurcast_targeted}
\alias{print.recurcast_targeted}
\alias{tidy.recurcast_targeted}
\alias{glance.recurcast_targeted}
\alias{coef.recurcast_targeted}
\alias{logLik.recurcast_targeted}
\title{Inspect targeted estimates and numerical diagnostics}
\usage{
\method{summary}{recurcast_targeted}(object, ...)
\method{print}{recurcast_targeted}(x, ...)
\method{tidy}{recurcast_targeted}(x, interval = "auto", level = .95,
  draws = 1000L, seed = 1L, ...)
\method{glance}{recurcast_targeted}(x, ...)
\method{coef}{recurcast_targeted}(object, ...)
\method{logLik}{recurcast_targeted}(object, ...)
}
\arguments{
\item{object, x}{A targeted recurcast fit.}
\item{interval, level, draws, seed}{Population interval controls passed to the targeted prediction method.}
\item{...}{Reserved additional arguments.}
}
\description{Reports target estimates, EIF diagnostics, regularity status and the fitted observation experiment. coef returns the population target vector, not event-intensity regression coefficients.}
\value{summary returns a diagnostic list, tidy a target table, glance a one-row diagnostic table, and coef the named target vector. print returns its argument invisibly. logLik raises a structured scope error because a targeted population functional does not supply a maximized likelihood for AIC.}
\seealso{\code{targeting_spec}, \code{predict.recurcast_targeted}}
