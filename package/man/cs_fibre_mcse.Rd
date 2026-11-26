\name{cs_fibre_mcse}
\alias{cs_fibre_mcse}
\title{Conditional integration error for a sampled observation fibre}
\usage{cs_fibre_mcse(fibre, path_loglik)}
\arguments{\item{fibre}{An object returned by \code{cs_coarsen_fibre}.}\item{path_loglik}{One finite complete-path log likelihood per expanded path at fixed parameter values.}}
\description{Computes a subject-specific Monte Carlo standard error from the inverse-proposal weighted integrand, including rejected samples as zeros.}
\value{A data frame with subject log likelihood, relative likelihood MCSE and its first-order log-likelihood MCSE approximation. MCSE is NA for deterministic quadrature, whose error requires a separate node-refinement check.}
\details{If M attempts give contributions a_b = likelihood_b / (M proposal_b), the integrated likelihood estimate is their sum. Its variance estimate is M/(M-1) times [sum(a_b squared) - sum(a_b) squared/M]. Stable rescaling avoids underflow. This quantifies numerical integration uncertainty conditional on parameters and observations, not the sampling variance of a fitted estimator or confidence coverage.}
\seealso{\code{\link{cs_coarsen_fibre}}}
