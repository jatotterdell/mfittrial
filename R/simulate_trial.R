#' Probability each column is maximum
#'
#' @param M A matrix of Monte Carlo draws
#' @export
prob_max <- function(M) {
  as.numeric(prop.table(table(factor(max.col(M), 1:ncol(M)))))
}

#' @export
pairwise_superiority_all <- function(mat, delta = 0, ...) {
  pmat <- pairwise_diff_all(mat, ...)
  apply(pmat, 2, function(x) mean(x > delta))
}

#' @export
simulate_trial <- function(
  n_seq = c(107, 214, 320),
  n_delay = 27,
  arms = 4,
  m_baseline = rep(44, arms),
  m_outcome = rep(44, arms),
  v_bo = diag(9^2, 2),
  p_alloc = rep(1 / arms, arms),
  kappa_sup = 0.95,
  kappa_act = 1 - kappa_sup,
  brar = FALSE,
  adjust_baseline = FALSE,
  M0 = rep(0, arms + adjust_baseline),
  S0 = diag(10^2, arms + adjust_baseline),
  A0 = 1e-2,
  B0 = 1e-2
) {

  if(min(n_seq) <= n_delay) stop("min(n_seq) must exceed n_delay to be valid.")

  N <- max(n_seq)
  K <- length(n_seq)
  n_new <- diff(c(0, n_seq))
  n_ints <- K - 1
  active <- rep(1, arms - 1)
  pair_comp <- arrangements::combinations(arms, 2)
  pair_name <- apply(pair_comp[, 2:1], 1, paste, collapse = " - ")
  stopped <- 0

  p_best_end <- stats::setNames(rep(0, arms), arms)
  p_pair_end <- stats::setNames(rep(0, length(pair_name)), pair_name)
  cMat <- multcomp::contrMat(rep(1, arms), "Tukey")[,]

  if(adjust_baseline) {
    form <- formula( ~ 0 + z + scale(x))
  } else {
    form <- formula( ~ 0 + z)
  }

  D <- sim_data(n_new[1], m_baseline, m_outcome, v_bo, p_alloc)
  X <- stats::model.matrix(form, data = D)[,]

  # If interims cycle through them
  if(n_ints > 0) {

    p_best_interim <- stats::setNames(rep(0, arms), 1:arms)
    p_pair_interim <- stats::setNames(rep(0, length(pair_name)), pair_name)
    n_interim <- matrix(0, arms, n_ints, dimnames = list(arm = 1:arms, interim = 1:n_ints))
    is_sup_interim <- stats::setNames(rep(FALSE, arms), 1:arms)

    for(k in seq_len(n_ints)) {

      # Undertake analysis excluding the delayed outcomes
      res_vb <- varapproxr::vb_lm(
        X[1:(n_seq[k] - n_delay), ], D[1:(n_seq[k] - n_delay), "y"], M0, S0, A0, B0)

      mu <- res_vb$mu[1:arms]
      Sigma <- res_vb$Sigma[1:arms, 1:arms]
      p_best_interim <- mvnorm_prob_each_best(mu, Sigma)
      p_pair_interim <- drop(1 - stats::pnorm(0, cMat %*% mu, sqrt(diag(cMat %*% Sigma %*% t(cMat)))))
      n_interim[, k] <- table(factor(D[1:(n_seq[k] - n_delay), "z"], 1:arms))

      # Check stopping rules
      is_sup_interim[] <- p_best_interim > kappa_sup
      any_sup_interim <- any(is_sup_interim)
      if(any_sup_interim) {

        stopped <- 1
        break

      } else {

        # Update allocation ratios if using BRAR
        if(brar) {
          active <- active * (p_best_interim[-1] > kappa_act)
          variance <- diag(Sigma)[1:arms]
          p_alloc <- brar(p_best_interim[-1], variance[-1], n_interim[-1, k], active)
        }

        # Generate next lot of data
        D <- rbind.data.frame(D, sim_data(n_new[k + 1], m_baseline, m_outcome, v_bo, p_alloc))
        X <- stats::model.matrix(form, data = D)[,]

      }
    } # end interim loop
  } else {
    p_best_interim <- stats::setNames(rep(NA, arms), 1:arms)
    p_pair_interim <- stats::setNames(rep(NA, length(pair_name)), pair_name)
    is_sup_interim <- stats::setNames(rep(NA, arms), 1:arms)
    any_sup_interim <- NA
  }

  # Undertake end analysis assuming we stop
  res_vb <- varapproxr::vb_lm(X, D$y, M0, S0, A0, B0)
  mu <- res_vb$mu[1:arms, 1]
  Sigma <- res_vb$Sigma[1:arms, 1:arms]
  p_best_end <- stats::setNames(mvnorm_prob_each_best(mu, Sigma), 1:arms)
  p_pair_end <- stats::setNames(drop(1 - stats::pnorm(0, cMat %*% mu, sqrt(diag(cMat %*% Sigma %*% t(cMat))))), pair_name)
  n <- table(factor(D[, "z"], 1:arms))

  is_sup <- stats::setNames(p_best_end > kappa_sup, 1:arms)
  any_sup <- any(is_sup)

  return(list(
    n_seq = n_seq,
    kappa_act = kappa_act,
    kappa_sup = kappa_sup,
    stopped = stopped,
    analysis = ifelse(stopped, k, K),
    any_sup = any_sup,
    is_sup = is_sup,
    fit = res_vb,
    n = n,
    p_best_end = p_best_end,
    p_pair_end = p_pair_end,
    is_sup_interim = is_sup_interim,
    any_sup_interim = any_sup_interim,
    p_best_interim = p_best_interim,
    p_pair_interim = p_pair_interim
  ))
}


#' @export
simulate_scenario <- function(sims, ...) {
  tibble::tibble(sim = 1:sims, trial = lapply(1:sims, function(i) simulate_trial(...)))
}


#' @export
simulate_trial_ord <- function(
  n_seq = c(107, 214, 320),
  n_delay = 27,
  arms = 4,
  items = 13,
  # m_baseline = do.call(rbind, rep(list(c(0.02, 0.07, 0.2, 0.65)), arms)),
  m_outcome = do.call(rbind, rep(list(c(0.02, 0.07, 0.2, 0.65)), arms)),
  r_items = diag(1, items),
  p_alloc = rep(1 / arms, arms),
  kappa_sup = 0.95,
  kappa_act = 0.05,
  brar = FALSE,
  adjust_baseline = FALSE,
  M0 = rep(0, arms + adjust_baseline),
  S0 = diag(10^2, arms + adjust_baseline),
  A0 = 1e-2,
  B0 = 1e-2
) {

  if(min(n_seq) <= n_delay) stop("min(n_seq) must exceed n_delay to be valid.")

  N <- max(n_seq)
  K <- length(n_seq)
  n_new <- diff(c(0, n_seq))
  n_ints <- K - 1
  active <- rep(1, arms - 1)
  pair_comp <- arrangements::combinations(arms, 2)
  pair_name <- apply(pair_comp[, 2:1], 1, paste, collapse = " - ")
  stopped <- 0

  p_best_end <- stats::setNames(rep(0, arms), arms)
  p_pair_end <- stats::setNames(rep(0, length(pair_name)), pair_name)
  cMat <- multcomp::contrMat(rep(1, arms), "Tukey")[,]

  if(adjust_baseline) {
    stop("Not implemented")
  } else {
    form <- formula( ~ 0 + z)
  }

  D <- as.data.frame(sim_data_ord(n_new[1], m_outcome, r_items, p_alloc, agg = TRUE))
  X <- stats::model.matrix(form, data = D)[,]

  # If interims cycle through them
  if(n_ints > 0) {

    p_best_interim <- stats::setNames(rep(0, arms), 1:arms)
    p_pair_interim <- stats::setNames(rep(0, length(pair_name)), pair_name)
    n_interim <- matrix(0, arms, n_ints, dimnames = list(arm = 1:arms, interim = 1:n_ints))
    is_sup_interim <- stats::setNames(rep(FALSE, arms), 1:arms)

    for(k in seq_len(n_ints)) {

      # Undertake analysis excluding the delayed outcomes
      res_vb <- varapproxr::vb_lm(
        X[1:(n_seq[k] - n_delay), ], D[1:(n_seq[k] - n_delay), "y"], M0, S0, A0, B0)

      mu <- res_vb$mu[1:arms]
      Sigma <- res_vb$Sigma[1:arms, 1:arms]
      p_best_interim <- mvnorm_prob_each_best(mu, Sigma)
      p_pair_interim <- drop(1 - stats::pnorm(0, cMat %*% mu, sqrt(diag(cMat %*% Sigma %*% t(cMat)))))
      n_interim[, k] <- table(factor(D[1:(n_seq[k] - n_delay), "z"], 1:arms))

      # Check stopping rules
      is_sup_interim[] <- p_best_interim > kappa_sup
      any_sup_interim <- any(is_sup_interim)
      if(any_sup_interim) {

        stopped <- 1
        break

      } else {

        # Update allocation ratios if using BRAR
        if(brar) {
          active <- active * (p_best_interim[-1] > kappa_act)
          variance <- diag(Sigma)[1:arms]
          p_alloc <- brar(p_best_interim[-1], variance[-1], n_interim[-1, k], active)
        }

        # Generate next lot of data
        D <- rbind(D, sim_data_ord(n_new[k + 1], m_outcome, r_items, p_alloc, agg = TRUE))
        X <- stats::model.matrix(form, data = D)[,]

      }
    } # end interim loop
  } else {
    p_best_interim <- stats::setNames(rep(NA, arms), 1:arms)
    p_pair_interim <- stats::setNames(rep(NA, length(pair_name)), pair_name)
    is_sup_interim <- stats::setNames(rep(NA, arms), 1:arms)
    any_sup_interim <- NA
  }

  # Undertake end analysis assuming we stop
  res_vb <- varapproxr::vb_lm(X, D$y, M0, S0, A0, B0)
  mu <- res_vb$mu[1:arms, 1]
  Sigma <- res_vb$Sigma[1:arms, 1:arms]
  p_best_end <- stats::setNames(mvnorm_prob_each_best(mu, Sigma), 1:arms)
  p_pair_end <- stats::setNames(drop(1 - stats::pnorm(0, cMat %*% mu, sqrt(diag(cMat %*% Sigma %*% t(cMat))))), pair_name)
  n <- table(factor(D[, "z"], 1:arms))

  is_sup <- stats::setNames(p_best_end > kappa_sup, 1:arms)
  any_sup <- any(is_sup)

  return(list(
    n_seq = n_seq,
    kappa_act = kappa_act,
    kappa_sup = kappa_sup,
    stopped = stopped,
    analysis = ifelse(stopped, k, K),
    any_sup = any_sup,
    is_sup = is_sup,
    fit = res_vb,
    n = n,
    p_best_end = p_best_end,
    p_pair_end = p_pair_end,
    is_sup_interim = is_sup_interim,
    any_sup_interim = any_sup_interim,
    p_best_interim = p_best_interim,
    p_pair_interim = p_pair_interim,
    p_alloc = p_alloc,
    D = D
  ))
}


#' @export
simulate_scenario_ord <- function(sims, ...) {
  tibble::tibble(sim = 1:sims, trial = lapply(1:sims, function(i) simulate_trial_ord(...)))
}


#' Simulate trial with ANOVA analysis
#'
#' @param n_seq The sample size at which analysis occurs
#' @param n_delay The number with delayed outcome
#' @param m_baseline The baseline response on ordinal scale
#' @param m_outcome The outcome response on ordinal scale
#' @param v_bo The correlation between baseline and outcome measurement
#' @param p_pref The distribution of preferences in the population
#' @param p_alloc The allocation probability to each arm conditional on preference
#' @param kappa_sup Threshold for superiority
#' @param kappa_act Threshold for remaining active
#' @param brar Use RAR?
#' @param adjust_baseline Include the baseline measurement in the model
#' @param M0 Prior mean
#' @param S0 Prior variance
#' @param A0 Prior scale for variance
#' @param B0 Prior shape for variance
#' @param drop_arms Allow dropping of arms?
#' @param fix_ctr Fixed allocation to control
#' @param verbose Be chatty
#'
#' @export
#'
simulate_trial_anova <- function(
  n_seq = c(107, 214, 320),
  n_delay = 27,
  m_baseline = rep(35, 16),
  m_outcome = rep(35, 16),
  v_bo = diag(8.3^2, 2),
  p_pref = rep(1/ 4, 4),
  p_alloc = matrix(0.25, 4, 4, dimnames = list(int = 1:4, pref = 1:4)),
  kappa_sup = 0.95,
  kappa_act = (1 - kappa_sup) / 3,
  brar = FALSE,
  adjust_baseline = FALSE,
  M0 = rep(0, 16 + adjust_baseline),
  S0 = diag(c(100^2, rep(10, 6), rep(1, 9))),
  A0 = 1e-2,
  B0 = 1e-2,
  drop_arms = FALSE,
  fix_ctr = NULL,
  verbose = FALSE,
  ...
) {

  if(min(n_seq) <= n_delay) stop("min(n_seq) must exceed n_delay to be valid.")
  if(!is.null(fix_ctr)) {
    p_alloc[1, ] <- fix_ctr
    p_alloc[-1, ] <- sapply(1:4, function(i) (1 - p_alloc[1, i]) * p_alloc[-1, i] / sum(p_alloc[-1, i]))
  }
  if(adjust_baseline) S0 <- diag(c(diag(S0), 100^2))

  # Track some values
  N <- max(n_seq)
  K <- length(n_seq)
  P <- ncol(X_con) - adjust_baseline
  n_new <- diff(c(0, n_seq))
  n_ints <- K - 1
  pair_comp <- arrangements::permutations(16, 2, replace = TRUE)
  pair_sub_comp <- arrangements::permutations(4, 2, replace = TRUE)
  pair_name <- apply(pair_comp[, 2:1], 1, paste, collapse = "-")
  pair_sub_name <- apply(pair_sub_comp[, 1:2], 1, paste, collapse = "-")

  # Track active/superior within subgroups
  is_active_sub <- matrix(1, 4, 4, dimnames = list(int = 1:4, pref = 1:4))
  is_sup_sub_interim <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))
  stopped <- 0

  # Posterior probabilities of interest #

  # Probability each intervention is best within subgroup/overall
  p_best_sub_interim <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))
  p_best_sub_end <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))
  p_best_interim <- stats::setNames(rep(0, 16), 1:16)
  p_best_end <- stats::setNames(rep(0, 16), 1:16)

  # Proability each intervention is best on average/beats control
  p_best_int_interim <- setNames(rep(0, 4), 1:4)
  p_best_int_end <- setNames(rep(0, 4), 1:4)

  # Pairwise comparisons within each subgroup and overall
  p_pair_sub_end <- matrix(0, 16, 4, dimnames = list(pair = pair_sub_name, pref = 1:4))
  p_pair_end <- setNames(rep(0, 16*16), pair_name)
  p_pair_int_end <- setNames(rep(0, 4*4), pair_sub_name)

  # Probability average exercise effect is greater than zero
  p_ex_beat_ctrl_sub_end <- stats::setNames(rep(0, 4), 1:4)
  p_ex_beat_ctrl_end <- 0

  # Generate initial data
  D <- sim_data_anova(n_new[1], m_baseline, m_outcome, v_bo, p_pref, p_alloc)
  X <- X_con[D$com, ]
  if(adjust_baseline) {
    X <- cbind(X, x = scale(D$x))
  }
  P <- ncol(X) - adjust_baseline

  # If interims cycle through them
  if(n_ints > 0) {

    n_interim <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))

    for(k in seq_len(n_ints)) {

      # Undertake analysis excluding the delayed outcomes
      Dint <- D[1:(n_seq[k] - n_delay), ]
      Xint <- X[1:(n_seq[k] - n_delay), ]
      yint <- Dint[, "y"]

      # Full model
      res_vb <- varapproxr::vb_lm(Xint, yint, M0, S0, A0, B0)

      if(verbose) cat("Model", ifelse(res_vb$converged, "converged", "failed to converge"), "\n")
      draws <- mvnfast::rmvn(1e4, res_vb$mu[1:P, ], res_vb$Sigma[1:P, 1:P])
      mu <- draws %*% X_con_t
      beta <- draws %*% Q_t

      # Summarise posterior probabilities of interest
      p_best_sub_interim <- sapply(1:4, function(i) prob_max(mu[, c(1,5,9,13) + (i - 1)]))
      p_best_interim <- prob_max(mu)
      p_best_int_interim <- prob_max(beta[, 2:5])
      p_int_beat_ctrl_interim <-
      n_sub_interim <- table(D[1:(n_seq[k] - n_delay), ]$int, D[1:(n_seq[k] - n_delay), ]$pref)
      is_sup_sub_interim <- p_best_sub_interim > kappa_sup
      any_sup_sub_interim <- matrixStats::colAnys(is_sup_sub_interim)

      # Update allocation ratios if using BRAR
      if(brar) {
        if(drop_arms) {
          is_active_sub <- is_active_sub * (p_best_sub_interim > kappa_act)
          if(!is.null(fix_ctr)) is_active_sub[1, ] <- 1
        } else {
          is_active_sub <- p_best_sub_interim > kappa_act
          if(!is.null(fix_ctr)) is_active_sub[1, ] <- 1
        }
        variance <- matrix(diag(X_con %*% res_vb$Sigma[1:P, 1:P] %*% t(X_con)), 4, 4, byrow = T)
        p_alloc[] <- sapply(1:4, function(i) {
          brar2(p_best_sub_interim[, i], variance[, i], n_sub_interim[, i], is_active_sub[, i], fix_ctr)
        })
        # If one arm is superior in the subgroup, allocate everyone to it, even control
        if(any(any_sup_sub_interim)) {
          for(i in 1:4) {
            if(any_sup_sub_interim[i]) {
              p_alloc[, i] <- is_sup_sub_interim[, i]
              is_active_sub[, i] <- is_active_sub[, i] * is_sup_sub_interim[, i]
            }
          }
        }

      }

      # If every subgroup has a superior treatment then stop
      if(all(any_sup_sub_interim)) {

        stopped <- 1
        break

      } else {

        # Generate next lot of data
        D <- rbind.data.frame(D, sim_data_anova(n_new[k + 1], m_baseline, m_outcome, v_bo, p_pref, p_alloc))
        X <- X_con[D$com, ]
        if(adjust_baseline) {
          X <- cbind(X, x = scale(D$x))
        }

      }
    } # end interim loop
  } else {
    p_best_sub_interim <- NA
    p_best_interim <- NA
    p_best_int_interim <- NA
    n_sub_interim <- NA
    is_sup_sub_interim <- NA
    any_sup_sub_interim <- NA

  }

  # Undertake end analysis assuming we stop
  res_vb <- varapproxr::vb_lm(X, D$y, M0, S0, A0, B0)
  if(verbose) cat("Model", ifelse(res_vb$converged, "converged", "failed to converge"), "\n")
  draws <- mvnfast::rmvn(1e4, res_vb$mu[1:P, ], res_vb$Sigma[1:P, 1:P])
  mu <- draws %*% X_con_t
  beta <- draws %*% Q_t
  n_sub <- table(int = D$int, pref = D$pref)

  # Calculate posterior probabilities at final analysis

  # Best in sub group
  p_best_sub_end[] <- sapply(1:4, function(i) prob_max(mu[, c(1,5,9,13) + (i - 1)]))
  is_sup_sub_end <- p_best_sub_end > kappa_sup
  any_sup_sub_end <- matrixStats::colAnys(is_sup_sub_end)

  # Best overall
  p_best_end[] <- prob_max(mu)
  is_sup_end <- p_best_end > kappa_sup
  any_sup_end <- any(is_sup_end)

  # Best average treatment effect
  p_best_int_end[] <- prob_max(beta[, 2:5])
  is_sup_int_end <- p_best_int_end > kappa_sup
  any_sup_int_end <- any(is_sup_int_end)

  # Average treatment vs control
  p_avg_int_beat_ctrl <- mean(matrixStats::rowMeans2(beta[, 3:5] - beta[, 2]) > 0)

  p_pair_sub_end[] <- sapply(1:4, function(i) pairwise_superiority_all(mu[, c(1,5,9,13) + (i - 1)], replace = T))
  # p_pair_end[] <- pairwise_superiority_all(mu, replace = T)
  p_pair_int_end[] <- pairwise_superiority_all(beta[, 2:5], replace = T)
  p_int_beat_ctrl_end <- p_pair_int_end[grepl("[2-4]-1", names(p_pair_int_end))]

  return(list(
    n_seq = n_seq,
    kappa_act = kappa_act,
    kappa_sup = kappa_sup,
    stopped = stopped,
    analysis = ifelse(stopped, k, K),
    fit = res_vb,
    mu = matrix(drop(X_con %*% res_vb$mu[1:P, 1]), 4, 4, dimnames = list(int = 1:4, pref = 1:4), byrow = T),
    beta = drop(Q %*% res_vb$mu[1:P, 1]),
    p_best_end = p_best_end,
    is_sup_end = is_sup_end,
    any_sup_end = any_sup_end,
    n_sub = n_sub,
    p_alloc = p_alloc,
    is_active_sub = is_active_sub,
    p_best_sub_end = p_best_sub_end,
    p_best_sub_interim = p_best_sub_interim,
    is_sup_sub_end = is_sup_sub_end,
    any_sup_sub_end = any_sup_sub_end,
    p_best_int_end = p_best_int_end,
    is_sup_int_end = is_sup_int_end,
    any_sup_int_end = any_sup_int_end,
    p_pair_sub_end = p_pair_sub_end,
    p_pair_int_end = p_pair_int_end,
    p_avg_int_beat_ctrl = p_avg_int_beat_ctrl,
    p_int_beat_ctrl_end = p_int_beat_ctrl_end
  ))
}


#' @export
simulate_scenario_anova <- function(sims, ...) {
  tibble::tibble(sim = 1:sims, trial = lapply(1:sims, function(i) simulate_trial_anova(...)))
}

#' @export
simulate_scenario_anova2 <- function(sims, ...) {
  tibble::tibble(sim = 1:sims, trial = mclapply(1:sims, function(i) simulate_trial_anova(...), mc.cores = 14))
}

#' Simulate trial with ANOVA analysis keeping all interims
#'
#' @param n_seq The sample size at which analysis occurs
#' @param n_delay The number with delayed outcome
#' @param m_baseline The baseline response on ordinal scale
#' @param m_outcome The outcome response on ordinal scale
#' @param v_bo The correlation between baseline and outcome measurement
#' @param p_pref The distribution of preferences in the population
#' @param p_alloc The allocation probability to each arm conditional on preference
#' @param kappa_sup Threshold for superiority
#' @param kappa_act Threshold for remaining active
#' @param brar Use RAR?
#' @param adjust_baseline Include the baseline measurement in the model
#' @param M0 Prior mean
#' @param S0 Prior variance
#' @param A0 Prior scale for variance
#' @param B0 Prior shape for variance
#' @param drop_arms Allow dropping of arms?
#' @param fix_ctr Fixed allocation to control
#' @param verbose Be chatty
#'
#' @export
#'
simulate_example_trial_anova <- function(
  n_seq = c(107, 214, 320),
  n_delay = 27,
  m_baseline = rep(0, 16),
  m_outcome = rep(0, 16),
  v_bo = diag(8.3^2, 2),
  p_pref = rep(1/ 4, 4),
  p_alloc = matrix(0.25, 4, 4, dimnames = list(int = 1:4, pref = 1:4)),
  kappa_sup = 0.95,
  kappa_act = (1 - kappa_sup) / 3,
  brar = FALSE,
  adjust_baseline = FALSE,
  M0 = rep(0, 16 + adjust_baseline),
  S0 = diag(c(100^2, rep(100^2, 15))),
  A0 = 1e-2,
  B0 = 1e-2,
  drop_arms = FALSE,
  fix_ctr = NULL,
  verbose = FALSE,
  ...
) {

  if(min(n_seq) <= n_delay) stop("min(n_seq) must exceed n_delay to be valid.")
  if(!is.null(fix_ctr)) {
    p_alloc[1, ] <- fix_ctr
    p_alloc[-1, ] <- sapply(1:4, function(i) (1 - p_alloc[1, i]) * p_alloc[-1, i] / sum(p_alloc[-1, i]))
  }
  if(adjust_baseline) S0 <- diag(c(diag(S0), 100^2))

  # Track some values
  N <- max(n_seq)
  K <- length(n_seq)
  P <- ncol(X_con) - adjust_baseline
  n_new <- diff(c(0, n_seq))
  n_ints <- K - 1
  pair_comp <- arrangements::permutations(16, 2, replace = TRUE)
  pair_sub_comp <- arrangements::permutations(4, 2, replace = TRUE)
  pair_name <- apply(pair_comp[, 2:1], 1, paste, collapse = "-")
  pair_sub_name <- apply(pair_sub_comp[, 1:2], 1, paste, collapse = "-")

  interim_results <- vector("list", n_ints)

  # Track active/superior within subgroups
  is_active_sub <- matrix(1, 4, 4, dimnames = list(int = 1:4, pref = 1:4))
  is_sup_sub_interim <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))
  stopped <- 0

  # Posterior probabilities of interest #

  # Probability each intervention is best within subgroup/overall
  p_best_sub <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))
  p_best_int <- setNames(rep(0, 4), 1:4)
  p_best     <- stats::setNames(rep(0, 16), 1:16)

  # Pairwise comparisons within each subgroup and overall
  p_pair_sub <- matrix(0, 16, 4, dimnames = list(pair = pair_sub_name, pref = 1:4))
  p_pair     <- setNames(rep(0, 16*16), pair_name)
  p_pair_int <- setNames(rep(0, 4*4), pair_sub_name)

  # Probability average exercise effect is greater than zero
  p_ex_beat_ctrl_sub <- stats::setNames(rep(0, 4), 1:4)
  p_ex_beat_ctrl     <- 0

  # Generate initial data
  D <- sim_data_anova(n_new[1], m_baseline, m_outcome, v_bo, p_pref, p_alloc)
  X <- X_con[D$com, ]
  if(adjust_baseline) {
    X <- cbind(X, x = scale(D$x))
  }
  P <- ncol(X) - adjust_baseline

  # If interims cycle through them
  if(n_ints > 0) {

    n_interim <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))

    for(k in seq_len(n_ints)) {

      # Undertake analysis excluding the delayed outcomes
      Dint <- D[1:(n_seq[k] - n_delay), ]
      Xint <- X[1:(n_seq[k] - n_delay), ]
      yint <- Dint[, "y"]

      # Full model
      res_vb <- varapproxr::vb_lm(Xint, yint, M0, S0, A0, B0)

      if(verbose) cat("Model", ifelse(res_vb$converged, "converged", "failed to converge"), "\n")
      draws <- mvnfast::rmvn(1e4, res_vb$mu[1:P, ], res_vb$Sigma[1:P, 1:P])
      mu <- draws %*% X_con_t
      beta <- draws %*% Q_t

      # Summarise posterior probabilities of interest
      p_best_sub <- sapply(1:4, function(i) prob_max(mu[, c(1,5,9,13) + (i - 1)]))
      p_best_int <- prob_max(beta[, 2:5])
      p_best     <- prob_max(mu)
      p_avg_int_beat_ctrl <- mean(matrixStats::rowMeans2(beta[, 3:5] - beta[, 2]) > 0)
      n_sub <- table(D[1:(n_seq[k] - n_delay), ]$int, D[1:(n_seq[k] - n_delay), ]$pref)
      is_sup_sub <- p_best_sub > kappa_sup
      any_sup_sub <- matrixStats::colAnys(is_sup_sub)
      p_pair_sub[] <- sapply(1:4, function(i) pairwise_superiority_all(mu[, c(1,5,9,13) + (i - 1)], replace = T))
      p_pair_int[] <- pairwise_superiority_all(beta[, 2:5], replace = T)
      p_int_beat_ctrl <- p_pair_int[grepl("[2-4]-1", names(p_pair_int))]

      # Update allocation ratios if using BRAR
      variance <- matrix(diag(X_con %*% res_vb$Sigma[1:P, 1:P] %*% t(X_con)), 4, 4, byrow = T)
      if(brar) {
        if(drop_arms) {
          is_active_sub <- is_active_sub * (p_best_sub > kappa_act)
          if(!is.null(fix_ctr)) is_active_sub[1, ] <- 1
        } else {
          is_active_sub <- p_best_sub > kappa_act
          if(!is.null(fix_ctr)) is_active_sub[1, ] <- 1
        }
        p_alloc[] <- sapply(1:4, function(i) {
          brar2(p_best_sub[, i], variance[, i], n_sub[, i], is_active_sub[, i], fix_ctr)
        })
        # If one arm is superior in the subgroup, allocate everyone to it, even control
        if(any(any_sup_sub)) {
          for(i in 1:4) {
            if(any_sup_sub[i]) {
              p_alloc[, i] <- is_sup_sub[, i]
              is_active_sub[, i] <- is_active_sub[, i] * is_sup_sub[, i]
            }
          }
        }
      }

      # If every subgroup has a superior treatment then stop
      if(all(any_sup_sub)) {
        stopped <- 1
      } else {

        # Generate next lot of data
        D <- rbind.data.frame(D, sim_data_anova(n_new[k + 1], m_baseline, m_outcome, v_bo, p_pref, p_alloc))
        X <- X_con[D$com, ]
        if(adjust_baseline) {
          X <- cbind(X, x = scale(D$x))
        }

      }
      interim_results[[k]] <- list(
        dat  = Dint,
        fit  = res_vb,
        mu   = matrix(drop(X_con %*% res_vb$mu[1:P, 1]), 4, 4, dimnames = list(int = 1:4, pref = 1:4), byrow = T),
        sigma = variance,
        beta = drop(Q %*% res_vb$mu[1:P, 1]),
        n_sub = n_sub,
        p_alloc = p_alloc,
        is_active_sub = is_active_sub,
        p_best_sub = p_best_sub,
        p_best_int = p_best_int,
        p_best = p_best,
        p_avg_int_beat_ctrl = p_avg_int_beat_ctrl,
        p_int_beat_ctrl = p_int_beat_ctrl,
        p_pair_sub = p_pair_sub,
        p_pari_int = p_pair_int
      )
      if(stopped == 1) break
    } # end interim loop
  } else {
    p_best_sub_interim <- NA
    p_best_interim <- NA
    p_best_int_interim <- NA
    n_sub_interim <- NA
    is_sup_sub_interim <- NA
    any_sup_sub_interim <- NA

  }

  # Undertake end analysis assuming we stop
  res_vb <- varapproxr::vb_lm(X, D$y, M0, S0, A0, B0)
  if(verbose) cat("Model", ifelse(res_vb$converged, "converged", "failed to converge"), "\n")
  draws <- mvnfast::rmvn(1e4, res_vb$mu[1:P, ], res_vb$Sigma[1:P, 1:P])
  mu <- draws %*% X_con_t
  beta <- draws %*% Q_t
  n_sub <- table(int = D$int, pref = D$pref)
  variance <- matrix(diag(X_con %*% res_vb$Sigma[1:P, 1:P] %*% t(X_con)), 4, 4, byrow = T)

  # Calculate posterior probabilities at final analysis

  # Best in sub group
  p_best_sub[] <- sapply(1:4, function(i) prob_max(mu[, c(1,5,9,13) + (i - 1)]))
  is_sup_sub   <- p_best_sub > kappa_sup
  any_sup_sub  <- matrixStats::colAnys(is_sup_sub)

  # Best overall
  p_best[] <- prob_max(mu)
  is_sup   <- p_best > kappa_sup
  any_sup  <- any(is_sup)

  # Best average treatment effect
  p_best_int[] <- prob_max(beta[, 2:5])
  is_sup_int   <- p_best_int > kappa_sup
  any_sup_int  <- any(is_sup_int)

  # Average treatment vs control
  p_avg_int_beat_ctrl <- mean(matrixStats::rowMeans2(beta[, 3:5] - beta[, 2]) > 0)
  p_pair_sub[]    <- sapply(1:4, function(i) pairwise_superiority_all(mu[, c(1,5,9,13) + (i - 1)], replace = T))
  p_pair_int[]    <- pairwise_superiority_all(beta[, 2:5], replace = T)
  p_int_beat_ctrl <- p_pair_int[grepl("[2-4]-1", names(p_pair_int))]

  final_results <- list(
    n_seq = n_seq,
    kappa_act = kappa_act,
    kappa_sup = kappa_sup,
    stopped = stopped,
    analysis = ifelse(stopped, k, K),
    dat  = D,
    fit = res_vb,
    mu = matrix(drop(X_con %*% res_vb$mu[1:P, 1]), 4, 4, dimnames = list(int = 1:4, pref = 1:4), byrow = T),
    sigma = variance,
    beta = drop(Q %*% res_vb$mu[1:P, 1]),
    p_best = p_best,
    is_sup = is_sup,
    any_sup = any_sup,
    n_sub = n_sub,
    p_alloc = p_alloc,
    is_active_sub = is_active_sub,
    p_best_sub = p_best_sub,
    is_sup_sub = is_sup_sub,
    any_sup_sub = any_sup_sub,
    p_best_int = p_best_int,
    is_sup_int = is_sup_int,
    any_sup_int = any_sup_int,
    p_pair_sub = p_pair_sub,
    p_pair_int = p_pair_int,
    p_avg_int_beat_ctrl = p_avg_int_beat_ctrl,
    p_int_beat_ctrl = p_int_beat_ctrl
  )

  return(list(
    interim_results,
    final_results))
}
