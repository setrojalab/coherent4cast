\name{rc_target}
\alias{rc_target}
\title{Cross-fitted targeting of a finite path law}
\usage{
rc_target(path, recorded, strata, features, folds = NULL,
  initial_laws = NULL, recording = NULL, control = list())
}
\arguments{
\item{path}{Integer complete-path indices from 1 to the number of feature rows. Values for unrecorded subjects may be NA and are never used.}
\item{recorded}{One binary complete-path recording indicator per independent subject.}
\item{strata}{Observed discrete covariate stratum per subject. This implementation uses stratum-specific nuisance fits.}
\item{features}{Numeric matrix with a row for each complete path and a column for each target functional. Columns must not create singular targeting information.}
\item{folds}{Optional subject-level cross-fitting labels. By default, subjects are assigned cyclically to up to five folds.}
\item{initial_laws}{Optional strictly positive initial probabilities, with fold-by-stratum rows and complete-path columns, fitted without the validation subjects.}
\item{recording}{Optional training-estimated recording probabilities for the same fold-by-stratum groups.}
\item{control}{List allowing folds, smoothing (default 0.5), recording_floor (0.05), maxit (60), and tolerance (1e-6/sqrt(n)).}
}
\description{Fits training-fold path and recording distributions, then solves a shared exponential tilting score. The same updated mixture law generates every reported target.}
\details{Finite-path MAR conditional on observed strata, recording positivity, independent subjects, and subject-separated nuisance fitting are required. Influence-function efficiency additionally requires the nuisance convergence and product remainder conditions stated in Version 5 Theorem 2. This is not arbitrary continuous-time hidden-history targeting. Recording estimates are truncated to [recording_floor, 0.99]. The recording floor is a numerical choice and is reported; it does not establish positivity or MAR. Initial laws supplied by users must be fitted without validation outcomes.}
\value{A list of estimates, standard errors, covariance, individual EIF values and diagnostics, updated laws, fold assignments, and convergence status. A failed score solution is retained as invalid rather than represented as valid inference.}
\examples{
features <- matrix(c(0, 1), ncol = 1)
path <- rep(c(1L, 2L, 2L, 1L), 10)
fit <- rc_target(path, rep(1, 40), rep(0, 40), features,
                 folds = rep(1:5, length.out = 40))
fit$estimate
fit$EIF_mean
}
