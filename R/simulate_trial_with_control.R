#' Simulate continuous outcome data
#'
#' @param nsubj
#' The number of participants to simulate
#' @param accrual
#' A function which determines rate of accrual in weeks,
#' default assumes mean of 3 per week
#' @param means
#' The mean response under each treatment
#' @param sigma
#' Variance of response (assumed constant over treatments)
#' @param trunc
#' Bound the outcome between 0 and 52 and make discrete.
#' @import truncdist
#' @export
simulate_continuous_outcome <- function(nsubj = 400,
                                        accrual = function(n) {
                                          cumsum(rexp(n, 3))
                                        },
                                        means = rep(40, 4),
                                        sigma = 5,
                                        trunc = FALSE) {
  dat <- data.table(id = 1:nsubj, t0 = accrual(nsubj))
  dat[, tt := t0 + 12]
  arms <- length(means)
  if (!trunc) {
    simf <- function(x) {
      rnorm(nsubj, means[x], sigma)
    }
  } else {
    simf <- function(x) {
      floor(
        rtrunc(nsubj, spec = "norm", a = 0, b = 53, mean = means[x], sd = sigma)
      )
    }
  }
  out <- vapply(seq_len(arms), simf, FUN.VALUE = matrix(0, nsubj, 1))
  out <- as.data.table(out)[, .(
    id = V1, t = V2, trt = factor(V3, levels = 1:arms), y = value
  )]
  out <- out[dat, on = .(id)]
  setcolorder(out, c("id", "t0", "tt", "t", "trt", "y"))
  return(out)
}


#' Estimate linear regression model using variational approximation.
#'
#' @param dat A dataset
#' @param prior A vector giving prior parameters:
#' intercept prior mean, intercept prior sd, effect prior sd,
#' variance prior df, variance prior scale.
#' Note that the assumed prior on the variance is Half-t(df, scale)
#' and on effects is Normal(0, effect_sd).
#' @param ctr A contrast specification
#' @importFrom bayestestR contr.orthonorm
#' @export
estimate_lm_model <- function(dat,
                              prior = c(
                                int_prior_mean = 40,
                                int_prior_sd = 5,
                                b_prior_sd = 5,
                                df = 3,
                                scale = 5
                              ),
                              ctr = contr.treatment) {
  int_prior_mean <- prior[1]
  int_prior_sd <- prior[2]
  b_prior_sd <- prior[3]
  a0 <- prior[4]
  b0 <- prior[5]
  mX <- model.matrix(~trt, data = dat, contrasts = list(trt = ctr))
  y <- dat$y
  mu0 <- c(int_prior_mean, 0, 0, 0)
  Sigma0 <- diag(c(int_prior_sd, rep(b_prior_sd, 3)))^2
  # Note, implies Half-t(df = a0, scale = b0) prior
  fit <- varapproxr::vb_lm(mX, y, mu0, Sigma0, a0 = a0, b0 = b0, prior = 2)
  return(list(mu = fit$mu, Sigma = fit$Sigma, X = mX))
}


#' Simulate a trial with fixed allocation to control group
#'
#' @param n_seq
#' The sequence of analyses by number of subjects with 12 week observation
#' @param dat The dataset for this particular simulation
#' @param eff_eps Effective decision threshold (relative to control)
#' @param sup_eps Superiority decision threshold (relative to active treatments)
#' @param fut_eps Futility decision threshold
#' @param delta Futility reference value
#' @param alloc Initial allocation ratio
#' @param brar Update allocation ratio using BRAR
#' @return A list (or data.table) giving the trial results
#' @export
simulate_trial_with_control <-
  function(n_seq = c(100, 200, 300, 400),
           dat = simulate_continuous_outcome(nsubj = max(n_seq)),
           eff_eps = 0.98,
           sup_eps = 0.98,
           fut_eps = 0.98,
           delta = 2,
           alloc = rep(1 / 4, 4),
           prior = c(
             int_prior_mean = 40,
             int_prior_sd = 5,
             b_prior_sd = 5,
             a0 = 3,
             b0 = 5
           ),
           brar = 2,
           brar_k = 0.5,
           allow_stopping = TRUE,
           make_dt = TRUE,
           ctr = contr.treatment,
           ...) {
    # Data setup
    K <- length(alloc)
    trt <- factor(1:K)
    X <- model.matrix(~trt, contrasts = list(trt = ctr))
    C <- attr(X, "contrasts")$trt
    Z <- cbind(-1, diag(1, 3))
    accdat <- unique(dat[, .(id, t0, tt)])
    trtdat <- data.table(id = 1:max(n_seq), trt = NA_character_)
    moddat <- NULL
    obsdat <- NULL
    # Interim setup
    N <- length(n_seq)
    n_add <- sapply(1:N, function(a) {
      length(unique(accdat$id[accdat$t0 <= accdat$tt[n_seq[a]]]))
    })
    n_new <- diff(c(0, n_add))
    t_seq <- sapply(1:N, function(a) accdat$tt[n_seq[a]])
    idobs <- cbind(c(1, n_seq[-length(n_seq)] + 1), n_seq)
    idenr <- cbind(c(1, n_add[-length(n_add)] + 1), n_add)
    # Output storage
    trtlabs <- paste0("trt", 0:(K - 1))
    n_enr <- matrix(0, N, K,
      dimnames = list(analysis = 1:N, treatment = trtlabs)
    )
    n_obs <- n_enr
    parlabs <- c("intercept", trtlabs)
    trt_mean <- matrix(0, N, K,
      dimnames = list(analysis = 1:N, treatment = trtlabs)
    )
    trt_var <- trt_mean
    eff_mean <- matrix(0, N, K - 1,
      dimnames = list(analysis = 1:N, treatment = trtlabs[-1])
    )
    eff_var <- eff_mean
    b_mean <- matrix(0, N, ncol(X),
      dimnames = list(analysis = 1:N, parameter = colnames(X))
    )
    p_alloc <- trt_mean
    p_supr <- eff_mean
    i_supr <- eff_mean
    i_infr <- eff_mean
    p_eff <- eff_mean
    p_fut <- eff_mean
    i_eff <- eff_mean
    i_inf <- eff_mean
    i_fut <- eff_mean
    i_acti <- matrix(1, N + 1, K - 1,
      dimnames = list(analysis = 0:N, treatment = trtlabs[-1])
    )
    final <- matrix(0, N, 1,
      dimnames = list(analysis = 1:N, final = "final")
    )
    stopped <- FALSE
    for (i in 1:N) {
      p_alloc[i, ] <- alloc
      final[i] <- stopped | i == N
      # Finish follow-up of enrolled participants if stopped
      if (stopped) {
        trtdat <- trtdat[!is.na(trt)]
        obsdat <- dat[trtdat, on = .(id, trt)]
        n_enr[i, ] <- trtdat[, .N, keyby = trt][["N"]]
        n_obs[i, ] <- obsdat[, .N, keyby = trt][["N"]]
      } else { # Otherwise enrol new participants
        trtenr <- factor(
          mass_weighted_urn_design(alloc, n_new[i], 4)$trt,
          levels = 1:K
        )
        trtdat[between(id, idenr[i, 1], idenr[i, 2]), trt := trtenr]
        moddat <- dat[trtdat[between(id, 1, idenr[i, 2])], on = .(id, trt)]
        obsdat <- moddat[tt <= t_seq[i]]
        n_enr[i, ] <- trtdat[
          between(id, 1, idenr[i, 2]), .N,
          keyby = trt
        ][["N"]]
        n_obs[i, ] <- obsdat[, .N, keyby = trt][["N"]]
      }
      fit <- estimate_lm_model(obsdat, prior = prior, ctr = ctr)
      mu <- fit$mu
      Sigma <- fit$Sigma
      Xmu <- X %*% mu
      XSXt <- X %*% Sigma %*% t(X)

      # Update the inferences
      b_mean[i, ] <- mu
      trt_mean[i, ] <- Xmu
      trt_var[i, ] <- diag(XSXt)
      eff_mean[i, ] <- drop(Z %*% Xmu)
      eff_var[i, ] <- diag(Z %*% XSXt %*% t(Z))
      mean_draws <- mvnfast::rmvn(
        2e4, Xmu[2:4], XSXt[2:4, 2:4]
      )

      # Conclusions carry forward
      # i.e. once decided superior/inferior, '
      # also superior/inferior at all future analyses
      if (i == 1) {
        # Probability active treatment better than control
        p_eff[i, ] <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
        p_fut[i, ] <- pnorm(delta, eff_mean[i, ], sqrt(eff_var[i, ]))

        # Is active treatment better than control
        i_eff[i, ] <- p_eff[i, ] > eff_eps

        # Is active treatment worse than control
        i_inf[i, ] <- p_eff[i, ] < 1 - eff_eps
        i_fut[i, ] <- p_fut[i, ] > fut_eps

        # Is active treatment superior
        p_supr[i, ] <- prob_supr(mean_draws)
        i_supr[i, ] <- p_supr[i, ] > sup_eps
        i_infr[i, ] <- p_supr[i, ] < (1 - sup_eps) / (K - 2)

        # If something superior, drop everything else
        if (any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1
        i_acti[i + 1, ] <- 1 - as.integer(
          (i_infr[i, ] == 1) | (i_inf[i, ] == 1)
        )
      } else {
        p_eff[i, ] <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
        p_fut[i, ] <- pnorm(delta, eff_mean[i, ], sqrt(eff_var[i, ]))
        i_eff[i, ] <- p_eff[i, ] > eff_eps
        # If ineffective, or previously ineffective
        i_inf[i, ] <- (p_eff[i, ] < 1 - eff_eps) | (i_inf[i - 1, ] == 1)
        i_fut[i, ] <- p_fut[i, ] > fut_eps
        # Only include active in superiority assessment
        p_supr[i, i_acti[i, ] == 1] <- prob_supr(
          mean_draws[, i_acti[i, ] == 1, drop = FALSE]
        )
        i_supr[i, ] <- (p_supr[i, ] > sup_eps) | (i_supr[i - 1, ] == 1)

        # If something superior, drop all other active treatments
        if (any(i_supr[i, ] == 1)) {
          i_infr[i, i_supr[i, ] != 1] <- 1
        } else {
          i_infr[i, ] <- p_supr[i, ] < min(1, (1 - sup_eps) /
            (sum(i_acti[i, ]) - 1), na.rm = TRUE)
        }

        # If something inferior or harmful, drop it,
        # if something previously dropped, it stays dropped
        i_acti[i + 1, ] <- 1 -
          as.integer(
            (i_infr[i, ] == 1) |
              (i_inf[i, ] == 1) |
              (i_acti[i, ] == 0)
          )
      }

      # Update allocations
      # - fix control at 1 / number of active arms
      n_active <- sum(i_acti[i + 1, ]) + 1
      if (brar == 1) {
        alloc[1] <- 1 / n_active
        if (alloc[1] == 1) {
          alloc[-1] <- 0
        } else {
          ratio <- (p_supr[i, ] * i_acti[i + 1, ] / n_enr[i, -1])^brar_k
          alloc[-1] <- (1 - alloc[1]) * ratio / sum(ratio)
        }
      } else if (brar == 2) {
        alloc[1] <- 1 / n_active
        ratio <- (p_supr[i, ] * i_acti[i + 1, ])^brar_k
        alloc[-1] <- (1 - alloc[1]) * ratio / sum(ratio)
      } else {
        alloc[1] <- 1 / n_active
        if (alloc[1] == 1) {
          alloc[-1] <- 0
        } else {
          alloc[-1] <- (1 - alloc[1]) * (alloc[-1] * i_acti[i + 1, ]) /
            sum(alloc[-1] * i_acti[i + 1, ])
        }
      }

      # If we stopped last analysis, we are done
      if (stopped) break
      # Should the next analysis be final?
      # - only stop if one active arm is superior and is better than control
      if (allow_stopping) {
        if (any(i_supr[i, ] & i_eff[i, ]) || all(i_acti[i + 1, ] == 0)) {
          stopped <- TRUE
        }
      }
    }

    idx <- 1:i
    # Results
    out <- list(
      p_alloc = p_alloc[idx, , drop = FALSE],
      n_enr = n_enr[idx, , drop = FALSE],
      n_obs = n_obs[idx, , drop = FALSE],
      trt_mean = trt_mean[idx, , drop = FALSE],
      trt_var = trt_var[idx, , drop = FALSE],
      eff_mean = eff_mean[idx, , drop = FALSE],
      eff_var = eff_var[idx, , drop = FALSE],
      p_eff = p_eff[idx, , drop = FALSE],
      p_fut = p_fut[idx, , drop = FALSE],
      i_eff = i_eff[idx, , drop = FALSE],
      i_inf = i_inf[idx, , drop = FALSE],
      i_fut = i_fut[idx, , drop = FALSE],
      p_supr = p_supr[idx, , drop = FALSE],
      i_supr = i_supr[idx, , drop = FALSE],
      i_infr = i_infr[idx, , drop = FALSE],
      i_acti = i_acti[idx + 1, , drop = FALSE]
    )
    if (make_dt) {
      trial_as_dt(out)
    } else {
      out
    }
  }



#' Simulate a trial with fixed allocation to control group
#'
#' @param n_seq
#' The sequence of analyses by number of subjects with 12 week observation
#' @param dat
#' The dataset for this particular simulation
#' @param eff_eps
#' Effective decision threshold (relative to control)
#' as function of fraction of information
#' @param sup_eps
#' Superiority decision threshold (relative to active treatments)
#' as function of fraction of information
#' @param fut_eps
#' Futility threshold as a function of fraction of trial information.
#' @param delta Futility reference value
#' @param alloc Initial allocation ratio
#' @param brar Update allocation ratio using BRAR
#' @return A list (or data.table) giving the trial results
#' @export
simulate_trial_with_control2 <-
  function(n_seq = c(100, 200, 300, 400),
           dat = simulate_continuous_outcome(nsubj = max(n_seq)),
           eff_eps = function(FI) 0.98,
           sup_eps = function(FI) 0.98,
           fut_eps = function(FI) 0.98,
           delta = 2,
           alloc = rep(1 / 4, 4),
           prior = c(
             int_prior_mean = 40,
             int_prior_sd = 5,
             b_prior_sd = 5,
             a0 = 3,
             b0 = 5
           ),
           brar = 2,
           brar_k = 0.5,
           allow_stopping = TRUE,
           make_dt = TRUE,
           ctr = contr.treatment,
           ..) {
    # Data setup
    K <- length(alloc)
    trt <- factor(1:K)
    X <- model.matrix(~trt, contrasts = list(trt = ctr))
    C <- attr(X, "contrasts")$trt
    Z <- cbind(-1, diag(1, 3))
    accdat <- unique(dat[, .(id, t0, tt)])
    trtdat <- data.table(id = 1:max(n_seq), trt = NA_character_)
    moddat <- NULL
    obsdat <- NULL
    # Interim setup
    N <- length(n_seq)
    n_add <- sapply(1:N, function(a) length(unique(accdat$id[accdat$t0 <= accdat$tt[n_seq[a]]])))
    n_new <- diff(c(0, n_add))
    t_seq <- sapply(1:N, function(a) accdat$tt[n_seq[a]])
    idobs <- cbind(c(1, n_seq[-length(n_seq)] + 1), n_seq)
    idenr <- cbind(c(1, n_add[-length(n_add)] + 1), n_add)
    # Output storage
    trtlabs <- paste0("trt", 0:(K - 1))
    n_enr <- matrix(0, N, K, dimnames = list(analysis = 1:N, treatment = trtlabs))
    n_obs <- n_enr
    parlabs <- c("intercept", trtlabs)
    trt_mean <- matrix(0, N, K, dimnames = list(analysis = 1:N, treatment = trtlabs))
    trt_var <- trt_mean
    eff_mean <- matrix(0, N, K - 1, dimnames = list(analysis = 1:N, treatment = trtlabs[-1]))
    eff_var <- eff_mean
    b_mean <- matrix(0, N, ncol(X), dimnames = list(analysis = 1:N, parameter = colnames(X)))
    p_alloc <- trt_mean
    p_supr <- eff_mean
    i_supr <- eff_mean
    i_infr <- eff_mean
    p_eff <- eff_mean
    p_fut <- eff_mean
    i_eff <- eff_mean
    i_inf <- eff_mean
    i_fut <- eff_mean
    i_acti <- matrix(1, N + 1, K - 1, dimnames = list(analysis = 0:N, treatment = trtlabs[-1]))
    final <- matrix(0, N, 1, dimnames = list(analysis = 1:N, final = "final"))
    stopped <- FALSE
    for (i in 1:N) {
      p_alloc[i, ] <- alloc
      final[i] <- stopped | i == N
      # Finish follow-up of enrolled participants if stopped
      if (stopped) {
        trtdat <- trtdat[!is.na(trt)]
        obsdat <- dat[trtdat, on = .(id, trt)]
        n_enr[i, ] <- trtdat[, .N, keyby = trt][["N"]]
        n_obs[i, ] <- obsdat[, .N, keyby = trt][["N"]]
      } else { # Otherwise enrol new participants
        trtenr <- factor(mass_weighted_urn_design(alloc, n_new[i], 4)$trt, levels = 1:K)
        trtdat[between(id, idenr[i, 1], idenr[i, 2]), trt := trtenr]
        moddat <- dat[trtdat[between(id, 1, idenr[i, 2])], on = .(id, trt)]
        obsdat <- moddat[tt <= t_seq[i]]
        n_enr[i, ] <- trtdat[between(id, 1, idenr[i, 2]), .N, keyby = trt][["N"]]
        n_obs[i, ] <- obsdat[, .N, keyby = trt][["N"]]
      }
      fit <- estimate_lm_model(obsdat, prior = prior, ctr = ctr)
      mu <- fit$mu
      Sigma <- fit$Sigma
      Xmu <- X %*% mu
      XSXt <- X %*% Sigma %*% t(X)

      # Update the inferences
      b_mean[i, ] <- mu
      trt_mean[i, ] <- Xmu
      trt_var[i, ] <- diag(XSXt)
      eff_mean[i, ] <- drop(Z %*% Xmu)
      eff_var[i, ] <- diag(Z %*% XSXt %*% t(Z))
      mean_draws <- mvnfast::rmvn(2e4, Xmu[2:4], XSXt[2:4, 2:4])

      # Conclusions carry forward
      # i.e. once decided superior/inferior, also superior/inferior at all future analyses
      eff_eps_i <- eff_eps(sum(n_enr[i, ]) / max(n_seq))
      fut_eps_i <- fut_eps(sum(n_enr[i, ]) / max(n_seq))
      sup_eps_i <- sup_eps(sum(n_enr[i, ]) / max(n_seq))
      if (i == 1) {
        # Probability active treatment better than control
        p_eff[i, ] <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
        p_fut[i, ] <- pnorm(delta, eff_mean[i, ], sqrt(eff_var[i, ]))

        # Is active treatment better than control
        i_eff[i, ] <- p_eff[i, ] > eff_eps_i

        # Is active treatment worse than control
        i_inf[i, ] <- p_eff[i, ] < 1 - eff_eps_i
        i_fut[i, ] <- p_fut[i, ] > fut_eps_i

        # Is active treatment superior
        p_supr[i, ] <- prob_supr(mean_draws)
        i_supr[i, ] <- p_supr[i, ] > sup_eps_i
        i_infr[i, ] <- p_supr[i, ] < (1 - sup_eps_i) / (K - 2)

        # If something superior, drop everything else
        if (any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1
        i_acti[i + 1, ] <- 1 - as.integer((i_infr[i, ] == 1) | (i_inf[i, ] == 1))
      } else {
        p_eff[i, ] <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
        p_fut[i, ] <- pnorm(delta, eff_mean[i, ], sqrt(eff_var[i, ]))
        i_eff[i, ] <- p_eff[i, ] > eff_eps_i
        # If ineffective, or previously ineffective
        i_inf[i, ] <- (p_eff[i, ] < 1 - eff_eps_i) | (i_inf[i - 1, ] == 1)
        i_fut[i, ] <- p_fut[i, ] > fut_eps_i

        # Only include active in superiority assessment
        p_supr[i, i_acti[i, ] == 1] <- prob_supr(mean_draws[, i_acti[i, ] == 1, drop = FALSE])
        i_supr[i, ] <- (p_supr[i, ] > sup_eps_i) | (i_supr[i - 1, ] == 1)
        i_infr[i, ] <- p_supr[i, ] < min(1, (1 - sup_eps_i) / (sum(i_acti[i, ]) - 1))

        # If something superior, drop all other active treatments
        if (any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1

        # If something inferior or harmful, drop it,
        # if something previously dropped, it stays dropped
        i_acti[i + 1, ] <- 1 - as.integer((i_infr[i, ] == 1) | (i_inf[i, ] == 1) | (i_acti[i, ] == 0))
      }

      # Update allocations
      # - fix control at 1 / number of active arms
      if (brar == 1) {
        n_active <- sum(i_acti[i + 1, ]) + 1
        alloc[1] <- 1 / n_active
        ratio <- sqrt(p_supr[i, ] * i_acti[i + 1, ] / n_enr[i, -1])
        alloc[-1] <- (1 - alloc[1]) * ratio / sum(ratio)
      } else if (brar == 2) {
        n_active <- sum(i_acti[i + 1, ]) + 1
        alloc[1] <- 1 / n_active
        ratio <- (p_supr[i, ] * i_acti[i + 1, ])^brar_k
        alloc[-1] <- (1 - alloc[1]) * ratio / sum(ratio)
      } else {
        n_active <- sum(i_acti[i + 1, ]) + 1
        alloc[1] <- 1 / n_active
        alloc[-1] <- (1 - alloc[1]) * (alloc[-1] * i_acti[i + 1, ]) / sum(alloc[-1] * i_acti[i + 1, ])
      }

      # If we stopped last analysis, we are done
      if (stopped) break
      # Should the next analysis be final?
      # - only stop if one active arm is superior and is better than control
      if (allow_stopping) {
        if (any(i_supr[i, ] & i_eff[i, ]) | all(i_acti[i + 1, ] == 0)) stopped <- TRUE
      }
    }

    idx <- 1:i
    # Results
    out <- list(
      p_alloc = p_alloc[idx, , drop = FALSE],
      n_enr = n_enr[idx, , drop = FALSE],
      n_obs = n_obs[idx, , drop = FALSE],
      trt_mean = trt_mean[idx, , drop = FALSE],
      trt_var = trt_var[idx, , drop = FALSE],
      eff_mean = eff_mean[idx, , drop = FALSE],
      eff_var = eff_var[idx, , drop = FALSE],
      p_eff = p_eff[idx, , drop = FALSE],
      p_fut = p_fut[idx, , drop = FALSE],
      i_eff = i_eff[idx, , drop = FALSE],
      i_inf = i_inf[idx, , drop = FALSE],
      i_fut = i_fut[idx, , drop = FALSE],
      p_supr = p_supr[idx, , drop = FALSE],
      i_supr = i_supr[idx, , drop = FALSE],
      i_infr = i_infr[idx, , drop = FALSE],
      i_acti = i_acti[idx + 1, , drop = FALSE]
    )
    if (make_dt) {
      trial_as_dt(out)
    } else {
      out
    }
  }



#' Simulate a trial with fixed allocation to control group.
#'
#' In this function, there are four possible options: BRAR on or off, Drop effective arms on or off.
#' The other adaptations are: always drop a harmful arm, and always stop if all arms harmful.
#'
#' @param n_seq
#' The sequence of analyses by number of subjects with 12 week observation
#' @param dat The dataset for this particular simulation
#' @param eff_eps Effective decision threshold (relative to control)
#' @param sup_eps Superiority decision threshold (relative to active treatments)
#' @param fut_eps Futility decision threshold
#' @param delta Futility reference value
#' @param alloc Initial allocation ratio
#' @param brar Update allocation ratio using BRAR
#' @param brar_k BRAR scaling factor, default is 0.5
#' @param dropeff Drop effective arms to focus on less effective ones
#' @return A list (or data.table) giving the trial results
#' @export
simulate_trial_with_control3 <-
  function(n_seq = c(200, 300, 400),
           dat = simulate_continuous_outcome(nsubj = max(n_seq)),
           eff_eps = 0.98,
           sup_eps = 0.98,
           fut_eps = 0.98,
           delta = 2,
           alloc = rep(1 / 4, 4),
           prior = c(
             int_prior_mean = 40,
             int_prior_sd = 5,
             b_prior_sd = 5,
             a0 = 3,
             b0 = 5
           ),
           brar = 2,
           brar_k = 0.5,
           brar_min = 0,
           drop = "none",
           allow_stopping = TRUE,
           make_dt = TRUE,
           ctr = contr.treatment,
           ...) {
    # Data setup
    K <- length(alloc)
    trt <- factor(1:K)
    X <- model.matrix(~trt, contrasts = list(trt = ctr))
    C <- attr(X, "contrasts")$trt
    Z <- cbind(-1, diag(1, 3))
    accdat <- unique(dat[, .(id, t0, tt)])
    trtdat <- data.table(id = 1:max(n_seq), trt = NA_character_)
    moddat <- NULL
    obsdat <- NULL
    # Interim setup
    N <- length(n_seq)
    n_add <- sapply(1:N, function(a) {
      length(unique(accdat$id[accdat$t0 <= accdat$tt[n_seq[a]]]))
    })
    n_new <- diff(c(0, n_add))
    t_seq <- sapply(1:N, function(a) accdat$tt[n_seq[a]])
    idobs <- cbind(c(1, n_seq[-length(n_seq)] + 1), n_seq)
    idenr <- cbind(c(1, n_add[-length(n_add)] + 1), n_add)
    # Output storage
    trtlabs <- paste0("trt", 0:(K - 1))
    n_enr <- matrix(0, N, K,
      dimnames = list(analysis = 1:N, treatment = trtlabs)
    )
    n_obs <- n_enr
    y_obs <- n_obs
    parlabs <- c("intercept", trtlabs)
    trt_mean <- matrix(0, N, K,
      dimnames = list(analysis = 1:N, treatment = trtlabs)
    )
    trt_var <- trt_mean
    eff_mean <- matrix(0, N, K - 1,
      dimnames = list(analysis = 1:N, treatment = trtlabs[-1])
    )
    eff_var <- eff_mean
    b_mean <- matrix(0, N, ncol(X),
      dimnames = list(analysis = 1:N, parameter = colnames(X))
    )
    p_alloc <- trt_mean
    p_supr <- eff_mean
    p_supr_act1 <- eff_mean
    p_supr_act2 <- eff_mean
    p_2best <- eff_mean
    i_supr <- eff_mean
    i_infr <- eff_mean
    p_eff <- eff_mean
    p_fut <- eff_mean
    i_eff <- eff_mean
    i_inf <- eff_mean
    i_fut <- eff_mean
    i_suprandeff <- eff_mean
    i_supr_act1 <- eff_mean
    i_2best <- eff_mean
    i_2bestandeff <- eff_mean
    i_acti <- matrix(1, N + 1, K - 1,
      dimnames = list(analysis = 0:N, treatment = trtlabs[-1])
    )
    final <- matrix(0, N, 1,
      dimnames = list(analysis = 1:N, final = "final")
    )
    stopped <- FALSE
    for (i in 1:N) {
      p_alloc[i, ] <- alloc
      final[i] <- stopped | i == N
      # Finish follow-up of enrolled participants if stopped
      if (stopped) {
        trtdat <- trtdat[!is.na(trt)]
        obsdat <- dat[trtdat, on = .(id, trt)]
        n_enr[i, ] <- trtdat[, .N, keyby = trt][["N"]]
        n_obs[i, ] <- obsdat[, .N, keyby = trt][["N"]]
        y_obs[i, ] <- obsdat[, .(y = mean(y)), keyby = trt][["y"]]
      } else { # Otherwise enrol new participants
        trtenr <- factor(
          mass_weighted_urn_design(alloc, n_new[i], 4)$trt,
          levels = 1:K
        )
        trtdat[between(id, idenr[i, 1], idenr[i, 2]), trt := trtenr]
        moddat <- dat[trtdat[between(id, 1, idenr[i, 2])], on = .(id, trt)]
        obsdat <- moddat[tt <= t_seq[i]]
        n_enr[i, ] <- trtdat[
          between(id, 1, idenr[i, 2]), .N,
          keyby = trt
        ][["N"]]
        n_obs[i, ] <- obsdat[, .N, keyby = trt][["N"]]
        y_obs[i, ] <- obsdat[, .(y = mean(y)), keyby = trt][["y"]]
      }
      fit <- estimate_lm_model(obsdat, prior = prior, ctr = ctr)
      mu <- fit$mu
      Sigma <- fit$Sigma
      Xmu <- X %*% mu
      XSXt <- X %*% Sigma %*% t(X)

      # Update the inferences
      b_mean[i, ] <- mu
      trt_mean[i, ] <- Xmu
      trt_var[i, ] <- diag(XSXt)
      eff_mean[i, ] <- drop(Z %*% Xmu)
      eff_var[i, ] <- diag(Z %*% XSXt %*% t(Z))
      mean_draws <- mvnfast::rmvn(
        2e4, Xmu[2:4], XSXt[2:4, 2:4]
      )

      if (i == 1) {
        # Probability active treatment better than control
        p_eff[i, ] <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
        p_fut[i, ] <- pnorm(delta, eff_mean[i, ], sqrt(eff_var[i, ]))
        # Is active treatment better than control
        i_eff[i, ] <- p_eff[i, ] > eff_eps
        # Is active treatment worse than control
        i_inf[i, ] <- p_eff[i, ] < 1 - eff_eps
        i_fut[i, ] <- p_fut[i, ] > fut_eps
        # Overall superiority
        p_supr[i, ] <- prob_supr(mean_draws)
        i_supr[i, ] <- p_supr[i, ] > sup_eps
        # Active superiority
        p_supr_act1[i, ] <- p_supr[i, ]
        i_supr_act1[i, ] <- p_supr_act1[i, ] > sup_eps
        # Superior and effective
        i_suprandeff[i, ] <- (p_supr[i, ] > sup_eps) & (p_eff[i, ] > eff_eps)
        # Second best
        p_2best[i, -which.max(p_supr[i, ])] <- prob_supr(mean_draws[, -which.max(p_supr[i, ]), drop = FALSE])
        i_2best[i, ] <- p_2best[i, ] > sup_eps
        i_2bestandeff[i, ] <- ((p_2best[i, ] > sup_eps) & (p_eff[i, ] > eff_eps) & any(i_suprandeff[i, ] == 1))
      } else {
        # Probability active treatment better than control
        p_eff[i, ] <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
        p_fut[i, ] <- pnorm(delta, eff_mean[i, ], sqrt(eff_var[i, ]))
        # If ineffective, or previously ineffective
        i_eff[i, ] <- (p_eff[i, ] > eff_eps) | (i_eff[i - 1, ] == 1)
        i_inf[i, ] <- (p_eff[i, ] < 1 - eff_eps) | (i_inf[i - 1, ] == 1)
        i_fut[i, ] <- (p_fut[i, ] > fut_eps) | (i_fut[i - 1, ] == 1)
        # Overall superiority
        p_supr[i, ] <- prob_supr(mean_draws)
        i_supr[i, ] <- p_supr[i, ] > sup_eps
        # Only include active in superiority assessment
        p_supr_act1[i, i_acti[i, ] == 1] <- prob_supr(mean_draws[, i_acti[i, ] == 1, drop = FALSE])
        i_supr_act1[i, ] <- p_supr_act1[i, ] > sup_eps
        # Ever superior and effective flag
        i_suprandeff[i, ] <- (p_supr_act1[i, ] > sup_eps & p_eff[i, ] > eff_eps & !any(i_suprandeff[i - 1, ] == 1)) | (i_suprandeff[i - 1, ] == 1)
        # Second best
        if (any(i_suprandeff[i, ] == 1)) {
          p_2best[i, -which(i_suprandeff[i, ] == 1)] <- prob_supr(mean_draws[, -which(i_suprandeff[i, ] == 1), drop = FALSE])
        } else {
          p_2best[i, -which.max(p_supr[i, ])] <- prob_supr(mean_draws[, -which.max(p_supr[i, ]), drop = FALSE])
        }
        i_2best[i, ] <- (p_2best[i, ] > sup_eps)
        i_2bestandeff[i, ] <- ((p_2best[i, ] > sup_eps) & (p_eff[i, ] > eff_eps) & any(i_suprandeff[i, ] == 1)) | (i_2bestandeff[i - 1, ] == 1)
      }

      if (drop == "eff") {
        # Drop only if:
        # - ineffective, or
        # - effective
        i_acti[i + 1, ] <- 1 - as.integer((i_eff[i, ] == 1) | (i_inf[i, ] == 1) | (i_acti[i, ] == 0))
      } else if (drop == "sup") {
        # Drop only if:
        # - superior and effective, or
        # - ineffective
        # Where by "superior" we mean superior to the remaining arms
        i_acti[i + 1, ] <- 1 - as.integer((i_suprandeff[i, ] == 1) | (i_2bestandeff[i, ] == 1) | (i_inf[i, ] == 1) | (i_acti[i, ] == 0))
      } else {
        # Drop only if:
        # - ineffective
        i_acti[i + 1, ] <- 1 - as.integer((i_inf[i, ] == 1) | (i_acti[i, ] == 0))
      }
      # Amongst still active arms what is Pr(best)?
      p_supr_act2[i, i_acti[i+1, ] == 1] <- prob_supr(
        mean_draws[, i_acti[i+1, ] == 1, drop = FALSE]
      )

      # Update allocations
      # - fix control at 1 / number of active arms
      n_active <- sum(i_acti[i + 1, ]) + 1
      if (brar == 1) {
        alloc[1] <- 1 / n_active
        if (alloc[1] == 1) {
          alloc[-1] <- 0
        } else {
          ratio <- (p_supr_act2[i, ] * i_acti[i + 1, ] / n_enr[i, -1])^brar_k
          alloc[-1] <- (1 - alloc[1]) * ratio / sum(ratio)
        }
      } else if (brar == 2) {
        alloc[1] <- 1 / n_active
        ratio <- (p_supr_act2[i, ] * i_acti[i + 1, ])^brar_k
        alloc[-1] <- (1 - alloc[1]) * ratio / sum(ratio)
      } else {
        alloc[1] <- 1 / n_active
        if (alloc[1] == 1) {
          alloc[-1] <- 0
        } else {
          alloc[-1] <- (1 - alloc[1]) * (alloc[-1] * i_acti[i + 1, ]) /
            sum(alloc[-1] * i_acti[i + 1, ])
        }
      }
      if(alloc[1] != 1) {
        alloc <- apply_min_brar(alloc, brar_min, c(1, i_acti[i + 1, ]))
      }

      # If we stopped last analysis, we are done
      if (stopped) break

      # Should the next analysis be final?
      # - only stop if all active arms have been dropped
      # - or if only continuing exercise intervention is effective
      if (allow_stopping) {
        stop_condition1 <- all(i_acti[i + 1, ] == 0)
        stop_condition2 <- ifelse(sum(i_acti[i + 1, ]) == 1, (p_eff[i, which(i_acti[i + 1, ] == 1)] > eff_eps), FALSE)
        if (stop_condition1 | stop_condition2) {
          stopped <- TRUE
          i_acti[i + 1, ] <- 0
          alloc <- rep(0, 4)
        }
      }
    }

    idx <- 1:i
    # Results
    out <- list(
      p_alloc = p_alloc[idx, , drop = FALSE],
      n_enr = n_enr[idx, , drop = FALSE],
      n_obs = n_obs[idx, , drop = FALSE],
      y_obs = y_obs[idx, , drop = FALSE],
      trt_mean = trt_mean[idx, , drop = FALSE],
      trt_var = trt_var[idx, , drop = FALSE],
      eff_mean = eff_mean[idx, , drop = FALSE],
      eff_var = eff_var[idx, , drop = FALSE],
      p_eff = p_eff[idx, , drop = FALSE],
      p_fut = p_fut[idx, , drop = FALSE],
      i_eff = i_eff[idx, , drop = FALSE],
      i_inf = i_inf[idx, , drop = FALSE],
      i_fut = i_fut[idx, , drop = FALSE],
      p_supr = p_supr[idx, , drop = FALSE],
      p_supr_act1 = p_supr_act1[idx, , drop = FALSE],
      p_supr_act2 = p_supr_act2[idx, , drop = FALSE],
      p_2best = p_2best[idx, , drop = FALSE],
      i_supr = i_supr[idx, , drop = FALSE],
      i_supr_act1 = i_supr_act1[idx, , drop = FALSE],
      i_2best = i_2best[idx, , drop = FALSE],
      i_suprandeff = i_suprandeff[idx, , drop = FALSE],
      i_2bestandeff = i_2bestandeff[idx, , drop = FALSE],
      i_acti = i_acti[idx + 1, , drop = FALSE]
    )
    if (make_dt) {
      trial_as_dt(out)
    } else {
      out
    }
  }