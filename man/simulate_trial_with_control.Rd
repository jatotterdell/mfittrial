% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_trial_with_control.R
\name{simulate_trial_with_control}
\alias{simulate_trial_with_control}
\title{Simulate a trial with fixed allocation to control group}
\usage{
simulate_trial_with_control(
  n_seq = c(100, 200, 300, 400),
  dat = simulate_continuous_outcome(nsubj = max(n_seq)),
  eff_eps = 0.98,
  sup_eps = 0.98,
  fut_eps = 0.98,
  delta = 2,
  alloc = rep(1/4, 4),
  prior = c(int_prior_mean = 40, int_prior_sd = 5, b_prior_sd = 5, a0 = 3, b0 = 5),
  brar = 2,
  brar_k = 0.5,
  allow_stopping = TRUE,
  make_dt = TRUE,
  ctr = contr.treatment,
  ...
)
}
\arguments{
\item{n_seq}{The sequence of analyses by number of subjects with 12 week observation}

\item{dat}{The dataset for this particular simulation}

\item{eff_eps}{Effective decision threshold (relative to control)}

\item{sup_eps}{Superiority decision threshold (relative to active treatments)}

\item{fut_eps}{Futility decision threshold}

\item{delta}{Futility reference value}

\item{alloc}{Initial allocation ratio}

\item{brar}{Update allocation ratio using BRAR}
}
\value{
A list (or data.table) giving the trial results
}
\description{
Simulate a trial with fixed allocation to control group
}
