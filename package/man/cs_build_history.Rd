\name{cs_build_history}
\alias{cs_build_history}
\title{Predictable history features}
\description{Calculate counts, recency and gaps using only earlier events.}
\usage{
cs_build_history(data, cap = 2L, decay = NULL) 
}
\arguments{
\item{data}{A data frame of complete counting-process intervals; one row is one observed interval.}
\item{cap}{Nonnegative integer count cap for predictable state features. Counts remain event counts; the cap applies to the intensity feature.}
\item{decay}{Optional nonnegative exponential decay rate for descriptive predictable recency features.}
}
\details{Counts and recency are evaluated at interval start, before the interval's own event. A missing gap means no prior event; it is not an imputed zero. Recency is descriptive unless the analyst supplies a corresponding sufficient propagation state.}
\value{A list of count, optional recency and gap matrices, with their feature settings.}
\seealso{\code{cs_loglik}, \code{cs_target}, \code{cs_bounds}}
