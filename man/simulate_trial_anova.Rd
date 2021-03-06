% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_trial.R
\name{simulate_trial_anova}
\alias{simulate_trial_anova}
\title{Simulate trial with ANOVA analysis}
\usage{
simulate_trial_anova(
  n_seq = c(107, 214, 320),
  n_delay = 27,
  m_baseline = rep(35, 16),
  m_outcome = rep(35, 16),
  v_bo = diag(8.3^2, 2),
  p_pref = rep(1/4, 4),
  p_alloc = matrix(0.25, 4, 4, dimnames = list(int = 1:4, pref = 1:4)),
  kappa_sup = 0.95,
  kappa_act = (1 - kappa_sup)/3,
  brar = FALSE,
  adjust_baseline = FALSE,
  M0 = rep(0, 16 + adjust_baseline),
  S0 = diag(c(100^2, rep(10, 6), rep(1, 9))),
  A0 = 0.01,
  B0 = 0.01,
  drop_arms = FALSE,
  fix_ctr = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{n_seq}{The sample size at which analysis occurs}

\item{n_delay}{The number with delayed outcome}

\item{m_baseline}{The baseline response on ordinal scale}

\item{m_outcome}{The outcome response on ordinal scale}

\item{v_bo}{The correlation between baseline and outcome measurement}

\item{p_pref}{The distribution of preferences in the population}

\item{p_alloc}{The allocation probability to each arm conditional on preference}

\item{kappa_sup}{Threshold for superiority}

\item{kappa_act}{Threshold for remaining active}

\item{brar}{Use RAR?}

\item{adjust_baseline}{Include the baseline measurement in the model}

\item{M0}{Prior mean}

\item{S0}{Prior variance}

\item{A0}{Prior scale for variance}

\item{B0}{Prior shape for variance}

\item{drop_arms}{Allow dropping of arms?}

\item{fix_ctr}{Fixed allocation to control}

\item{verbose}{Be chatty}
}
\description{
Simulate trial with ANOVA analysis
}
