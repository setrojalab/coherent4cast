\name{rc_covariate_model}
\alias{rc_covariate_model}
\title{Gaussian longitudinal covariate submodels}
\description{Fit random-intercept and optional random-slope models for numeric covariates.}
\usage{
rc_covariate_model(data, unit = "id", time = "time", variables, 
    random_slope = TRUE, method = c("REML", "ML")) 
}
\arguments{
\item{data}{A data frame of observed longitudinal measurements, with one row per subject and measurement time. Event intervals are not automatically biomarker measurements.}
\item{unit}{Name of the subject identifier column in the measurement table.}
\item{time}{Name of the numeric measurement-time column, in the same time unit as the event model.}
\item{variables}{Names of numeric longitudinal covariates to model.}
\item{random_slope}{Include a random time slope as well as a random intercept.}
\item{method}{Restricted maximum likelihood (REML) or maximum likelihood (ML) for the Gaussian longitudinal submodel.}
}
\details{Uses nlme maximum or restricted maximum likelihood after centering and scaling time. Conditional Gaussian random-effect draws propagate future trajectories through the event model. Separate variables use conditionally independent submodels. The working joint law factorizes future covariate evolution and the event intensity given its predictable covariate trajectory; no feedback from future events to covariates is assumed. It is not a joint maximum-likelihood latent measurement-error survival model. Train on an appropriate estimation sample; prediction conditions only on measurements available at the landmark. Model-parameter uncertainty is separate from trajectory integration Monte Carlo error.}
\value{An rc_covariate_model containing fitted submodels, fixed and random effects, residual variances, scaling and assumptions.}
\seealso{\code{rc_loglik}, \code{rc_target}, \code{rc_bounds}}
