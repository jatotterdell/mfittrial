% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_trial_with_control.R
\name{simulate_continuous_outcome}
\alias{simulate_continuous_outcome}
\title{Simulate continuous outcome data}
\usage{
simulate_continuous_outcome(
  nsubj = 400,
  accrual = function(n) {
     cumsum(rexp(n, 3))
 },
  means = rep(40, 4),
  sigma = 5,
  trunc = FALSE
)
}
\arguments{
\item{nsubj}{The number of participants to simulate}

\item{accrual}{A function which determines rate of accrual in weeks,
default assumes mean of 3 per week}

\item{means}{The mean response under each treatment}

\item{sigma}{Variance of response (assumed constant over treatments)}

\item{trunc}{Bound the outcome between 0 and 52 and make discrete.}
}
\description{
Simulate continuous outcome data
}
