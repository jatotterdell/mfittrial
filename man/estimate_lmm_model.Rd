% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_trial_longitudinal.R
\name{estimate_lmm_model}
\alias{estimate_lmm_model}
\title{Estimate the model}
\usage{
estimate_lmm_model(
  dat,
  zero_sum = F,
  scale = "none",
  int_prior_sd = 10,
  b_prior_sd = 10
)
}
\arguments{
\item{dat}{The data (as a data.table)}

\item{zero_sum}{Use sum-to-zero coding}

\item{scale}{Character string indicating how to centre/scale responses, default is `none`, other options are `all` or `baseline`.}

\item{int_prior_sd}{Prior standard deviation on intercept parameter}

\item{b_prior_sd}{Prior standard deviation used for all other population parameters}
}
\value{
A list giving the mean and variance of the variational approximation to the fixed effects
}
\description{
Estimate the model
}
