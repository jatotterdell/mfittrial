# Two models
# - one where preference is ignored
# - one with full preference interaction weighted by preference distribution
# Generate the data and the run two trials:
# - one which ignores preference
# - one which treats preference as subgroups

#' Generate potential outcomes if a participant were to receive each treatment
#' @param n The number of participants
#' @param m A matrix giving the mean response on each treatment (cols) by preference (rows)
#' @param s The variability of response assumed constant across preferences/treatments
#' @param p Distribution of preferences in the population
#' @export
generate_preference_data <- function(
  n = 320,
  m = matrix(c(-3,-1,-1,-1,-2,2,0,0,-2,0,2,0,-1,1,1,3), 4, 4),
  s = 10,
  p = c(0.2, 0.7, 0.05, 0.05)
  ) {
  pref <- sample.int(length(p), n, replace = T, prob = p)
  y    <- matrix(rnorm(ncol(m)*n, m[pref, ], s), n, ncol(m))
  time <- seq(0, 156, length.out = n)
  list(
    id   = 1:n,
    pref = pref,
    t0   = time,
    t1   = time + 12,
    y    = y,
    p    = p
  )
}

#' Probability each column is maximum
#'
#' @param M A matrix of Monte Carlo draws
prob_supr <- function(M) {
  as.numeric(prop.table(table(factor(max.col(M), 1:ncol(M)))))
}


#' Estimate linear regression model
#'
#' @param X Design matrix
#' @param y Response vector
#' @param m0 Prior mean
#' @param s0 Prior SD
#' @param a0 Prior shape
#' @param b0 Prior shape
#' @param ... Other arguments to `vb_lm`
estimate_model <- function(
  X,
  y,
  m0 = rep(0, ncol(X)),
  s0 = rep(10, ncol(X)),
  a0 = 1e-2,
  b0 = 1e-2,
  ...
) {
  fit <- varapproxr::vb_lm(X, y, m0, diag(s0^2), a0, b0, ...)
  return(fit)
}

#' Update linear regression model
#'
#' @param fit A previous `vb_lm` fit
#' @param X New design matrix data
#' @param y New response vector data
#' @param ... Other arguments to `update_vb_lm`
update_model <- function(
  fit,
  X,
  y,
  ...
) {
  fitnew <- varapproxr::update_vb_lm(fit, X, y, ...)
  return(fitnew)
}


#' Simulate trial without any knowledge of preference
#'
#' @param dat Data generated from `generate_preference_data`
#' @param n_seq Sequence of interim analyses
#' @param alloc Initial allocation probabilities
#' @param brar Use BRAR?
#' @param ... Other arguments to `vb_lm` functions
#' @export
simulate_standard_trial <- function(
  dat = generate_preference_data(),
  n_seq = c(100, 150, 200, 250, nrow(dat$y)),
  alloc = rep(1/ncol(dat$y), ncol(dat$y)),
  sup_eps = 0.95,
  eff_eps = 0.95,
  brar = FALSE,
  ...
) {

  # Interim setup
  N     <- length(n_seq)
  n_enr <- sapply(1:N, function(a) length(dat$t0[dat$t0 <= dat$t1[n_seq[a]]]))
  n_new <- diff(c(0, n_enr))
  idobs <- cbind(c(1, n_seq[-length(n_seq)] + 1), n_seq)
  idenr <- cbind(c(1, n_enr[-length(n_enr)] + 1), n_enr)
  yobs  <- NULL
  Xobs  <- NULL
  trt   <- NULL
  trtobs<- NULL
  prf   <- dat$pref


  # Model design matrix
  K  <- ncol(dat$y)
  C <- as.matrix(Matrix::bdiag(1, eigen(diag(1, K) - 1/K)$vector[, 1:(K-1)]))
  Xd <- cbind(1, diag(1, 4))
  X <- Xd %*% C
  invC <- MASS::ginv(C)
  Z <- cbind(-1, diag(1, K - 1))

  # Output storage
  trtlabs <- paste0("trt", 1:K)
  parlabs <- c("intercept", trtlabs)
  n_enr  <- matrix(0, N, ncol(X), dimnames = list(analysis = 1:N, treatment = trtlabs))
  n_obs  <- n_enr
  y_obs  <- n_enr
  trt_mean <- n_enr
  trt_var  <- n_enr
  eff_mean <- matrix(0, N, ncol(X) - 1, dimnames = list(analysis = 1:N, treatment = trtlabs[-1]))
  eff_var  <- eff_mean
  p_supr <- trt_mean
  i_supr <- trt_mean
  i_infr <- trt_mean
  p_eff  <- eff_mean
  i_eff  <- eff_mean
  i_inf  <- eff_mean
  i_acti <- matrix(1, N+1, ncol(X), dimnames = list(analysis = 0:N, treatment = trtlabs))
  b_mean <- matrix(0, N, ncol(Xd), dimnames = list(analysis = 1:N, parameter = parlabs))
  b_var  <- b_mean
  p_pair <- matrix(NA, N, ncol(X), dimnames = list(analysis = 1:N, treatment = trtlabs))

  final <- matrix(0, N, 1, dimnames = list(analysis = 1:N, final = 'final'))
  stopped <- FALSE
  for(i in 1:N) {

    # Update the data
    if(stopped) {
      trtnew  <- trt[idobs[i,1]:idenr[i-1,2]]
      trtobs  <- c(trtobs, trtnew)
      yobsnew <- dat$y[cbind(idobs[i,1]:idenr[i-1,2], trtnew)]
      yobs    <- c(yobs, yobsnew)
      Xobsnew <- X[trtnew, ]
      Xobs    <- rbind(Xobs, Xobsnew)
    } else {
      trtenr  <- sample.int(K, n_new[i], TRUE, prob = alloc)
      trt     <- c(trt, trtenr)
      trtnew  <- trt[idobs[i,1]:idobs[i,2]]
      trtobs  <- c(trtobs, trtnew)
      yobsnew <- dat$y[cbind(idobs[i,1]:idobs[i,2], trtnew)]
      yobs    <- c(yobs, yobsnew)
      Xobsnew <- X[trtnew, ]
      Xobs    <- rbind(Xobs, Xobsnew)
    }
    final[i] <- stopped | i == N
    n_enr[i, ] <- as.numeric(table(factor(trt, 1:K)))
    n_obs[i, ] <- as.numeric(table(factor(trtobs, 1:K)))
    y_obs[i, ] <- aggregate(yobs, by = list(factor(trtobs, 1:K)), mean, drop = F)$x

    mean_y <- mean(yobs)
    sd_y   <- sd(yobs)
    y_std  <- (yobs - mean_y) / sd_y

    # Update the model
    # fit <- estimate_model(Xobs, yobs, ...)
    fit <- estimate_model(Xobs, y_std, ...)

    # Back transform parameters to original scale
    mu <- drop(fit$mu) * sd_y
    mu[1] <- mu[1] + mean_y
    Sigma <- fit$Sigma * sd_y^2

    # Update the inferences
    trt_mean[i, ] <- drop(X %*% mu)
    trt_var[i, ]  <- diag(X %*% Sigma %*% t(X))
    eff_mean[i, ] <- drop(Z %*% trt_mean[i, ])
    eff_var[i, ] <- diag((Z %*% X) %*% Sigma %*% t(Z %*% X))
    b_mean[i, ] <- drop(C %*% mu)
    b_var[i, ]  <- diag(C %*% Sigma %*% t(C))

    draws <- mvnfast::rmvn(2e4, mu, Sigma)
    means <- draws %*% t(X)

    # Conclusions carry forward
    # i.e. once decided effective, also effective at all future analyses
    if(i == 1) {
      p_eff[i, ]    <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
      i_eff[i, ]    <- p_eff[i, ] > eff_eps
      i_inf[i, ]    <- p_eff[i, ] < 1 - eff_eps
      p_supr[i, ]   <- prob_supr(means)
      i_supr[i, ]   <- p_supr[i, ] > sup_eps
      i_infr[i, ]   <- p_supr[i, ] < (1 - sup_eps) / (K - 1)
      if(any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1
      i_acti[i+1, ] <- 1 - i_infr[i, ]
    } else {
      p_eff[i, ]    <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
      i_eff[i, ]    <- (p_eff[i, ] > eff_eps) | (i_eff[i-1,] == 1)
      i_inf[i, ]    <- (p_eff[i, ] < 1 - eff_eps) | (i_inf[i-1, ] == 1)
      p_supr[i, i_acti[i, ] == 1] <- prob_supr(means[, i_acti[i, ] == 1, drop = F])
      i_supr[i, ]   <- (p_supr[i, ] > sup_eps) | (i_supr[i-1, ] == 1)
      i_infr[i, ]   <- p_supr[i, ] < min(1, (1 - sup_eps) / (sum(i_acti[i, ]) - 1))
      if(any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1
      i_acti[i+1, ] <- 1 - i_infr[i, ]
    }
    # Pairwise comparisons with current best
    best <- which.max(p_supr[i, ])
    P <- diag(1, K)
    P[, best] <- P[, best] - 1
    pair_mean <- (P %*% X %*% mu)
    pair_var  <- diag((P %*% X) %*% Sigma %*% t(P %*% X))
    p_pair[i, -best] <- pnorm(0, pair_mean, sqrt(pair_var))[-best]

    # Update allocations
    if(brar) {
      ratio <- sqrt(p_supr[i, ] * i_acti[i+1, ] / n_enr[i, ])
      alloc <- ratio / sum(ratio)
    } else {
      alloc <- (alloc * i_acti[i+1, ]) / sum(alloc * i_acti[i+1, ])
    }

    if(stopped) break

    if(any(i_supr[i, ] == 1)) stopped <- TRUE
  }

  idx <- 1:i

  # Results
  list(
    final = final[idx, , drop = F],
    n_enr = n_enr[idx, , drop = F],
    n_obs = n_obs[idx, , drop = F],
    y_obs = y_obs[idx, , drop = F],
    trt_mean = trt_mean[idx, , drop = F],
    trt_var  = trt_var[idx, , drop = F],
    eff_mean = eff_mean[idx, , drop = F],
    eff_var = eff_var[idx, , drop = F],
    p_eff = p_eff[idx, , drop = F],
    i_eff = i_eff[idx, , drop = F],
    i_inf = i_inf[idx, , drop = F],
    p_supr = p_supr[idx, , drop = F],
    i_supr = i_supr[idx, , drop = F],
    i_infr = i_infr[idx, , drop = F],
    i_acti = i_acti[idx + 1, , drop = F],
    b_mean = b_mean[idx, , drop = F],
    b_var  = b_var[idx, , drop = F],
    p_pair = p_pair[idx, , drop = F]
  )
}

#' Simulate trial accounting for preference
#'
#' @param dat Data generated from `generate_preference_data`
#' @param n_seq Sequence of interim analyses
#' @param alloc Initial allocation probabilities (matrix)
#' @param brar Use BRAR?
#' @param min_decision The minimum sample size in a subgroup to make decision/switch on BRAR
#' @param ... Other arguments to `vb_lm` functions
#' @export
simulate_preference_trial <- function(
  dat = generate_preference_data(),
  n_seq = c(100, 150, 200, 250, nrow(dat$y)),
  alloc = matrix(1/ncol(dat$y), length(dat$p), ncol(dat$y)),
  sup_eps = 0.95,
  eff_eps = 0.95,
  brar = FALSE,
  min_decision = 20,
  ...
) {

  Q <- length(dat$p)

  # Interim setup
  N     <- length(n_seq)
  n_enr <- sapply(1:N, function(a) length(dat$t0[dat$t0 < dat$t1[n_seq[a]]]))
  n_new <- diff(c(0, n_enr))
  idobs <- cbind(c(1, n_seq[-length(n_seq)] + 1), n_seq)
  idenr <- cbind(c(1, n_enr[-length(n_enr)] + 1), n_enr)
  yobs  <- NULL
  Xobs  <- NULL
  trt   <- NULL
  trtobs<- NULL
  prf   <- NULL
  prfobs<- NULL
  com   <- NULL
  comobs<- NULL


  # Model design matrix (preference design is data dependent)
  K  <- ncol(dat$y)
  Xtrt <- kronecker(diag(1, K), rep(1, Q))
  Xprf <- kronecker(rep(1, K), diag(1, Q))
  Xcom <- matrix(apply(Xtrt, 2, function(x) x * Xprf), Q*K)
  Xd <- cbind(1, Xtrt, Xprf, Xcom)
  Ctrt <- eigen(diag(1, K) - 1/K)$vector[, 1:(K-1)]
  Z <- cbind(-1, diag(1, K - 1))
  Zcom <- cbind(kronecker(rep(1,K-1),diag(-1, Q)), kronecker(diag(1, Q),diag(1, K-1)))

  map_com <- function(pref, treat) (treat - 1)*K + pref

  # Output storage
  trtlabs <- paste0("trt", 1:K)
  prflabs <- paste0("prf", 1:Q)
  comlabs <- paste(rep(trtlabs, each = Q), rep(prflabs, times = K), sep = "_")
  comefflabs <- paste(rep(trtlabs[-1], each = Q), rep(prflabs, times = K-1), sep = "_")
  parlabs <- c("intercept", trtlabs, prflabs, comlabs)
  n_enr_prf <- matrix(0, N, Q, dimnames = list(analysis = 1:N, preference = prflabs))
  n_obs_prf <- n_enr_prf
  y_obs_prf <- n_enr_prf
  prf_mean  <- n_enr_prf
  prf_var   <- n_enr_prf
  n_enr_trt <- matrix(0, N, K, dimnames = list(analysis = 1:N, treatment = trtlabs))
  n_obs_trt <- n_enr_trt
  y_obs_trt <- n_enr_trt
  trt_mean  <- n_enr_trt
  trt_var   <- n_enr_trt
  n_enr_com <- matrix(0, N, Q*K, dimnames = list(analysis = 1:N, treatment = comlabs))
  n_obs_com <- n_enr_com
  y_obs_com <- n_enr_com
  com_mean <- matrix(0, N, Q*K, dimnames = list(analysis = 1:N, combination = comlabs))
  com_var  <- com_mean
  eff_mean <- matrix(0, N, K - 1, dimnames = list(analysis = 1:N, treatment = trtlabs[-1]))
  eff_var <- eff_mean

  p_eff  <- eff_mean
  i_eff  <- eff_mean
  i_inf  <- eff_mean
  p_supr_trt <- trt_mean
  i_supr_trt <- trt_mean
  i_infr_trt <- trt_mean
  i_acti_trt <- matrix(1, N+1, K, dimnames = list(analysis = 0:N, treatment = trtlabs))
  p_supr_com <- com_mean
  i_supr_com <- com_mean
  i_infr_com <- com_mean
  p_eff_com  <- matrix(0, N, K*Q - Q, dimnames = list(analysis = 1:N, combination = comefflabs))
  i_eff_com  <- p_eff_com
  i_inf_com  <- p_eff_com
  i_acti_com <- matrix(1, N+1, Q*K, dimnames = list(analysis = 0:N, combination = comlabs))
  b_mean <- matrix(0, N, Q+K+Q*K+1, dimnames = list(analysis = 1:N, parameter = parlabs))
  b_var  <- b_mean

  final <- matrix(0, N, 1, dimnames = list(analysis = 1:N, final = 'final'))
  stopped <- FALSE
  for(i in 1:N) {

    # Update the data
    if(stopped) {
      prfnew  <- prf[idobs[i,1]:idenr[i-1,2]]
      prfobs  <- c(prfobs, prfnew)
      trtnew  <- trt[idobs[i,1]:idenr[i-1,2]]
      trtobs  <- c(trtobs, trtnew)
      comnew  <- map_com(prfnew, trtnew)
      comobs  <- c(comobs, comnew)
      yobsnew <- dat$y[cbind(idobs[i,1]:idenr[i-1,2], trtnew)]
      yobs    <- c(yobs, yobsnew)

      # Determine preference distribution estimates
      n_obs_prf[i, ] <- as.numeric(table(factor(prfobs, 1:Q)))
      p_prf      <- n_obs_prf[i, ] / sum(n_obs_prf[i, ])

      # Constrain by preference distribution
      Cprf <- eigen(diag(1, Q) - matrix(kronecker(p_prf, rep(1, K)), Q, K))[[2]][, 1:(Q-1)]
      Ccom <- kronecker(Ctrt, Cprf)
      C    <- as.matrix(Matrix::bdiag(1, Ctrt, Cprf, Ccom))
      X    <- Xd %*% C
      Xobs <- X[comobs, ]

    } else {
      prfenr  <- dat$pref[idenr[i,1]:idenr[i,2]]
      prf     <- c(prf, prfenr)
      prfnew  <- prf[idobs[i,1]:idobs[i,2]]
      prfobs  <- c(prfobs, prfnew)
      trtenr  <- sapply(1:n_new[i], function(a) sample.int(K, 1, prob = alloc[prfenr[a], ]))
      trt     <- c(trt, trtenr)
      trtnew  <- trt[idobs[i,1]:idobs[i,2]]
      trtobs  <- c(trtobs, trtnew)
      com     <- c(com, map_com(prfenr, trtenr))
      comnew  <- map_com(prfnew, trtnew)
      comobs  <- c(comobs, comnew)
      yobsnew <- dat$y[cbind(idobs[i,1]:idobs[i,2], trtnew)]
      yobs    <- c(yobs, yobsnew)

      # Determine preference distribution estimates
      n_obs_prf[i, ] <- as.numeric(table(factor(prfobs, 1:Q)))
      p_prf      <- n_obs_prf[i, ] / sum(n_obs_prf[i, ])

      # Constrain by preference distribution
      Cprf <- eigen(diag(1, Q) - matrix(kronecker(p_prf, rep(1, K)), Q, K))[[2]][, 1:(Q-1)]
      Ccom <- kronecker(Ctrt, Cprf)
      C    <- as.matrix(Matrix::bdiag(1, Ctrt, Cprf, Ccom))
      X    <- Xd %*% C
      Xobs <- X[comobs, ]
    }

    final[i] <- stopped | i == N

    # Aggregate summaries
    n_enr_prf[i, ] <- as.numeric(table(factor(prf, 1:Q)))
    y_obs_prf[i, ] <- aggregate(yobs, by = list(factor(prfobs, 1:Q)), mean, drop = F)$x
    n_enr_trt[i, ] <- as.numeric(table(factor(trt, 1:K)))
    n_obs_trt[i, ] <- as.numeric(table(factor(trtobs, 1:K)))
    y_obs_trt[i, ] <- aggregate(yobs, by = list(factor(trtobs, 1:K)), mean, drop = F)$x
    n_enr_com[i, ] <- as.numeric(table(factor(com, 1:(Q*K))))
    n_obs_com[i, ] <- as.numeric(table(factor(comobs, 1:(Q*K))))
    y_obs_com[i, ] <- aggregate(yobs, by = list(factor(comobs, 1:(Q*K))), mean, drop = F)$x

    mean_y <- mean(yobs)
    sd_y   <- sd(yobs)
    y_std  <- (yobs - mean_y) / sd_y

    # Update the model
    fit <- estimate_model(Xobs, y_std, ...)
    if(!fit$converged) stop("Failed to converge")

    # Back transform parameters to original scale
    mu <- drop(fit$mu) * sd_y
    mu[1] <- mu[1] + mean_y
    Sigma <- fit$Sigma * sd_y^2

    # Update the inferences
    b_mean[i, ] <- drop(C %*% mu)
    b_var[i, ]  <- diag(C %*% Sigma %*% t(C))
    com_mean[i, ] <- drop(X %*% mu)
    com_var[i, ]  <- diag(X %*% Sigma %*% t(X))
    trt_mean[i, ] <- drop(cbind(1, Ctrt) %*% mu[1:K])
    trt_var[i, ] <- diag(cbind(1, Ctrt) %*% Sigma[1:K, 1:K] %*% t(cbind(1, Ctrt)))
    prf_mean[i, ] <- drop(cbind(1, Cprf) %*% mu[c(1, (K+1):(K+Q-1))])
    prf_var[i, ] <- diag(
      cbind(1, Cprf) %*%
      fit$Sigma[c(1,(K+1):(K+Q-1)), c(1,(K+1):(K+Q-1))] %*%
      t(cbind(1, Cprf)))
    eff_mean[i, ] <-  drop(Z %*% cbind(1, Ctrt) %*% mu[1:K])
    eff_var[i, ] <- diag(Z %*% cbind(1, Ctrt) %*% Sigma[1:K, 1:K] %*% t(Z %*% cbind(1, Ctrt)))
    effcom_mean <- Zcom %*% X %*% mu
    effcom_var  <- diag((Zcom %*% X) %*% Sigma %*% t(Zcom %*% X))

    draws <- mvnfast::rmvn(2e4, mu, Sigma)
    betas <- draws %*% t(C)
    means <- draws %*% t(X)

    # Conclusions carry forward
    if(i == 1) {
      p_eff[i, ]    <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
      i_eff[i, ]    <- p_eff[i, ] > eff_eps
      i_inf[i, ]    <- p_eff[i, ] < 1 - eff_eps
      for(a in 1:Q) {
        id <- map_com(a, 1:K)
        ideff <- map_com(a, 2:K) - Q
        p_supr_com[i, id] <- prob_supr(means[, id])
        p_eff_com[i, ideff] <- 1 - pnorm(0, effcom_mean[ideff], sqrt(effcom_var[ideff]))
        if(n_obs_prf[i, a] >= min_decision) {
          i_supr_com[i, id] <- p_supr_com[i, id] > sup_eps
          i_infr_com[i, id] <- p_supr_com[i, id] < (1 - sup_eps) / (K - 1)
          if(any(i_supr_com[i, id])) i_infr_com[i, id][i_supr_com[i, id] != 1] <- 1
          i_acti_com[i+1, id] <- 1 - i_infr_com[i, id]
        }
        i_eff_com[i, ideff] <-  p_eff_com[i, ideff] > eff_eps
        i_inf_com[i, ideff] <-  p_eff_com[i, ideff] < 1 - eff_eps
      }
      p_supr_trt[i, ] <- prob_supr(betas[, 1] + betas[, 2:5])
      i_supr_trt[i, ] <- p_supr_trt[i, ] > sup_eps
      i_infr_trt[i, ] <- p_supr_trt[i, ] <  (1 - sup_eps) / (K - 1)
    } else {
      p_eff[i, ]    <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
      i_eff[i, ]    <- (p_eff[i, ] > eff_eps) | (i_eff[i-1, ] == 1)
      i_inf[i, ]    <- (p_eff[i, ] < 1 - eff_eps) | (i_inf[i-1,] == 1)
      for(a in 1:Q) {
        id <- map_com(a, 1:K)
        ideff <- map_com(a, 2:K) - Q
        act_id <- i_acti_com[i, id] == 1
        p_supr_com[i, id][act_id] <-
          prob_supr(means[, id][, act_id, drop = F])
        p_eff_com[i, ideff] <- 1 - pnorm(0, effcom_mean[ideff], sqrt(effcom_var[ideff]))
        if(n_obs_prf[i, a] >= min_decision) {
          i_supr_com[i, id] <- p_supr_com[i, id] >= sup_eps | (i_supr_com[i-1, id] == 1)
          i_infr_com[i, id] <- ((p_supr_com[i, id] < (1 - sup_eps) / (sum(i_acti_com[i, id]) - 1)) |
                                  (i_infr_com[i-1, id] == 1)) & (i_supr_com[i-1, id] != 1)
          if(any(i_supr_com[i, id])) {
            i_infr_com[i, id][i_supr_com[i, id] != 1] <- 1
            # i_infr_com[i, id][i_supr_com[i, id] == 1] <- 0
          }
          i_acti_com[i+1, id] <- 1 - i_infr_com[i, id]
        }
        i_eff_com[i, ideff] <-  p_eff_com[i, ideff] > eff_eps
        i_inf_com[i, ideff] <-  p_eff_com[i, ideff] < 1 - eff_eps
      }
      p_supr_trt[i, ] <- prob_supr(betas[, 1] + betas[, 2:5])
      i_supr_trt[i, ] <- p_supr_trt[i, ] > sup_eps
      i_infr_trt[i, ] <- p_supr_trt[i, ] <  (1 - sup_eps) / (K - 1)
    }

    # Update allocations
    if(brar) {
      for(a in 1:Q) {
        # Only BRAR if exceed minimum sample size
        if(n_obs_prf[i, a] >= min_decision) {
          id <- map_com(a, 1:K)
          ratio <- sqrt(p_supr_com[i, id] * i_acti_com[i+1, id] / (n_enr_com[i, id] + 1))
          alloc[a, ] <- ratio / sum(ratio)
        }
      }
    } else {
      for(a in 1:Q) {
        id <- map_com(a, 1:K)
        alloc[a, ] <- alloc[a, ] * i_acti_com[i+1, id] / sum(alloc[a, ] * i_acti_com[i+1, id])
      }
    }

    if(any(is.na(alloc))) return(list(alloc, i_supr_com, i_infr_com, i_acti_com))

    if(stopped) break

    if(all(sapply(1:Q, function(a) any(i_supr_com[i, map_com(a, 1:K)] == 1)))) stopped <- TRUE
  }


  idx <- 1:i

  # Results
  list(
    final = final[idx, , drop = F],
    n_enr_prf = n_enr_prf[idx, , drop = F],
    n_obs_prf = n_obs_prf[idx, , drop = F],
    y_obs_prf = y_obs_prf[idx, , drop = F],
    n_enr_trt = n_enr_trt[idx, , drop = F],
    n_obs_trt = n_obs_trt[idx, , drop = F],
    y_obs_trt = y_obs_trt[idx, , drop = F],
    n_enr_com = n_enr_com[idx, , drop = F],
    n_obs_com = n_obs_com[idx, , drop = F],
    y_obs_com = y_obs_com[idx, , drop = F],
    prf_mean = prf_mean[idx, , drop = F],
    prf_var  = prf_var[idx, , drop = F],
    trt_mean = trt_mean[idx, , drop = F],
    trt_var  = trt_var[idx, , drop = F],
    eff_mean = eff_mean[idx, , drop = F],
    eff_var = eff_var[idx, , drop = F],
    com_mean = com_mean[idx, , drop = F],
    com_var  = com_var[idx, , drop = F],
    p_eff = p_eff[idx, , drop = F],
    i_eff = i_eff[idx, , drop = F],
    i_inf = i_inf[idx, , drop = F],
    p_supr_com = p_supr_com[idx, , drop = F],
    i_supr_com = i_supr_com[idx, , drop = F],
    i_infr_com = i_infr_com[idx, , drop = F],
    i_acti_com = i_acti_com[idx + 1, , drop = F],
    p_eff_com = p_eff_com[idx, , drop = F],
    i_eff_com = i_eff_com[idx, , drop = F],
    i_inf_com = i_inf_com[idx, , drop = F],
    p_supr_trt = p_supr_trt[idx, , drop = F],
    i_supr_trt = i_supr_trt[idx, , drop = F],
    i_infr_trt = i_infr_trt[idx, , drop = F],
    b_mean = b_mean[idx, , drop = F],
    b_var  = b_var[idx, , drop = F]
  )
}
