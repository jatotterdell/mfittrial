library(designr)
library(varapproxr)
library(data.table)

#' Probability each column is maximum
#'
#' @param M A matrix of Monte Carlo draws
prob_supr <- function(M) {
  as.numeric(prop.table(table(factor(max.col(M), 1:ncol(M)))))
}


#' Create the design matrix
#'
#' @param zero_sum Use sum-to-zero permutation invariant coding
#' @return The design matrix
#' @export
make_design <- function(zero_sum = FALSE) {
  Xt  <- kronecker(diag(1, 3), rep(1, 4))
  Xxt <- diag(1, 12)
  # Constrain to be equal at baseline and sum to zero at each post-randomisation time point
  if(zero_sum) {
    C <- eigen(diag(1, 4) - 1/4)$vectors[, 1:3]
  } else {
    C <- contr.treatment(4)
  }
  Xx <- kronecker(rep(1, 3), C)
  # XC  <- cbind(1, rbind(0, 0, 0, 0, cbind(Xt, kronecker(diag(1, 3), C))))
  # colnames(XC) <- c("t0", "t1", "t2", "t3",
  #                   "t1x1", "t1x2", "t1x3", "t2x1", "t2x2", "t2x3", "t3x1", "t3x2", "t3x3")
  XC <- cbind(1, rbind(0, cbind(Xt, Xt * Xx[, 1], Xt * Xx[, 2], Xt * Xx[, 3])))
  colnames(XC) <- c("t0", "t1", "t2", "t3",
                    "t1x1", "t2x1", "t3x1", "t1x2", "t2x2", "t3x2", "t1x3", "t2x3", "t3x3")
  rownames(XC) <- c("t0",
                    "t1x0", "t1x1", "t1x2", "t1x3",
                    "t2x0", "t2x1", "t2x2", "t2x3",
                    "t3x0", "t3x1", "t3x2", "t3x3")
  return(XC)
}


#' Simulate accrual data
#'
#' @param nsubj The number of subjects
#' @param accrual An accrual function which returns `nsubj` ordered randomisation times in weeks
#' @param dropout A drop-out function which returns `nsubj` times-to-drop-out, one for each subject.
#' `Inf` means no drop-out (not implemented)
#' @return A data.table giving the accrual data
#' @export
#' @import data.table
simulate_accrual_data <- function(
  nsubj = 200,
  accrual = function(n) round(cumsum(rexp(n, 3))),
  dropout = function(n) sample(c(1, Inf), n, replace = T, prob = c(0.2, 0.8))
) {
  dat <- data.table(
    id = 1:nsubj,
    t0 = accrual(nsubj)
  )
  dat[, td  := t0 + dropout(nsubj)]
  return(dat)
}


#' Simulate outcome data
#'
#' @param nsims Number of sets of outcome data to simulate
#' @param nsubj Number of subjects
#' @param accrual Accrual function for use in `simulate_accrual_data`
#' @param beta Fixed effects parameters
#' @param sigma_e Variance of residuals for outcome
#' @param sigma_u Variance of subject specific intercepts
#' @return A data.table giving outcome data
#' @export
#' @import parallel
simulate_outcome_data <- function(
  nsims = 50,
  nsubj = 200,
  accrual = function(n) 1:n,
  beta = c(0, 2, 4, 6, 0, -2.5, -2.5, -2.5, 0, 2.5, 2.5, 2.5, 0, 0, 0, 0),
  sigma_e = 5,
  sigma_u = 1
) {
  accdat <- simulate_accrual_data(nsubj, accrual)
  simdat <- accdat[data.table::CJ(id = 1:nsubj, t = 0:3, x = 0:3), on = .(id)]
  simdat[, `:=`(
    trt  = x,
    x    = factor(ifelse(t == 0, 0L, x)),
    time = fcase(t == 0, 0, t == 1, 4, t == 2, 8, t == 3, 12),
    t    = factor(t)
  )][, tt := t0 + time]
  rtsimmat <- do.call(data.table:::cbind.data.table, parallel::mclapply(1:nsims, function(i) {
    designr::simLMM(
      ~ 0 + t + t:x + (1 | id),
      data = simdat,
      Fixef = beta,
      VC_sd = list(sigma_u, sigma_e),
      empirical = FALSE,
      verbose = FALSE
    )
  }, mc.cores = parallel::detectCores() - 1))
  return(data.table:::cbind.data.table(simdat, rtsimmat))
}

#' Simulate outcome data
#'
#' @param nsims Number of sets of outcome data to simulate
#' @param nsubj Number of subjects
#' @param means The mean at each time point in each arm
#' @param covar The co-variance structure between time points (assumed equal across all arms)
#' @return A data.table giving outcome data
#' @export
simulate_general_outcome_data <- function(
  nsims = 50,
  nsubj = 400,
  means = matrix(0, 4, 4, dimnames = list(t = 0:3, a = 0:3)),
  covar = {R <- matrix(0.5, 4, 4); diag(R) <- 1; 5^2*R},
  ...
) {
  accdat <- simulate_accrual_data(nsubj, ...)
  obsv <- nrow(means)
  arms <- ncol(means)
  out  <- vapply(seq_len(nsims), function(z) {
    vapply(seq_len(arms), function(j) {
      mvnfast::rmvn(nsubj, means[, j], covar)
    }, FUN.VALUE = matrix(0, nsubj, obsv))
    }, FUN.VALUE = array(0, dim = c(nsubj, obsv, arms)))
  out <- as.data.table(out)[, .(
    id = V1, t = factor(V2 - 1), trt = V3 - 1, sim = V4, y = value
  )][order(sim, id, trt, t)]
  out <- dcast(out, id + t + trt ~ forcats::fct_inorder(paste0("y", sim)), value.var = 'y')
  out <- out[accdat, on = .(id)]
  out[, `:=`(
    tt = fcase(t == 0, t0, t == 1, t0 + 4, t == 2, t0 + 8, t == 3, t0 + 12),
    x  = factor(ifelse(t == 0, 0L, trt))
  )]
  setcolorder(out, c("id", "trt", "t", "x", "t0", "tt", "td"))
  return(out)
}


#' Estimate the model
#'
#' @param dat The data (as a data.table)
#' @param zero_sum Use sum-to-zero coding
#' @param scale Character string indicating how to centre/scale responses, default is `none`, other options are `all` or `baseline`.
#' @param int_prior_sd Prior standard deviation on intercept parameter
#' @param b_prior_sd Prior standard deviation used for all other population parameters
#' @return A list giving the mean and variance of the variational approximation to the fixed effects
#' @export
#' @import bayestestR
estimate_lmm_model <- function(dat, zero_sum = F, scale = "none", int_prior_sd = 10, b_prior_sd = 10) {
  if(zero_sum) {
    mX <- model.matrix( ~ t + t:x, data = dat, contrasts.arg = list(x = 'contr.bayes'))[, -c(5, 9, 13)]
  } else {
    mX <- model.matrix( ~ t + t:x, data = dat)[, -c(5, 9, 13)]
  }
  mZ <- unname(model.matrix( ~ 0 + factor(id), data = dat)[,])
  y  <- dat[[grep("y", names(dat), value = T)]]

  if(scale == "all") {
    y_mean <- mean(y)
    y_sd   <- sd(y)
  } else if(scale == "baseline") {
    y_mean <- mean(y[dat$t == 0])
    y_sd   <- sd(y[dat$t == 0])
  } else {
    y_mean <- 0
    y_sd   <- 1
  }
  y_std  <- (y - y_mean) / y_sd
  fit_vb <- varapproxr::vb_lmm_randint(
    mX,
    mZ,
    y_std,
    mu_beta = rep(0, ncol(mX)),
    sigma_beta = diag(c(int_prior_sd^2, rep(b_prior_sd^2, ncol(mX) - 1))),
    mu = rep(0, ncol(mZ) + ncol(mX)),
    sigma = diag(1, ncol(mZ) + ncol(mX)),
    Aeps = 1e-2, Beps = 1e-2,
    Au = 1e-2, Bu = 1e-2,
    tol = 1e-4,
    maxiter = 500
  )
  if(!fit_vb$converged) warning("Failed to converge")
  # Back transform parameters to original scale
  mu <- drop(fit_vb$mu)[1:ncol(mX)] * y_sd
  mu[1] <- mu[1] + y_mean
  Sigma <- fit_vb$sigma[1:ncol(mX), 1:ncol(mX)] * y_sd^2

  return(list(mu = mu, Sigma = Sigma))
}


#' Simulate a longitudinal trial
#'
#' @param n_seq The sequence of analyses by number of subjects with 12 week observation
#' @param dat The dataset for this particular simulation
#' @param eff_eps Effective decision threshold
#' @param sup_eps Superiority decision threshold
#' @param alloc Initial allocation ratio
#' @param brar Update allocation ratio using BRAR
#' @return A list (or data.table) giving the trial results
#' @export
simulate_longitudinal_trial <- function(
  n_seq = c(100, 150, 200, 250, 400),
  dat = simulate_general_outcome_data(nsims = 1, nsubj = max(n_seq)),
  eff_eps = 0.95,
  sup_eps = 0.95,
  alloc   = rep(1/4, 4),
  brar    = TRUE,
  make_dt = TRUE,
  zero_sum = FALSE,
  scale = FALSE,
  ...
) {
  # Data setup
  K      <- length(alloc)
  X      <- make_design(zero_sum = zero_sum)
  C      <- X[-1, -(1:4)]
  Z      <- cbind(matrix(0, 3, 9), cbind(-1, diag(1, 3)))
  accdat <- unique(dat[t == 3, .(id, t0, tt)])
  trtdat <- data.table(id = 1:max(n_seq), trt = 99)
  moddat <- NULL
  obsdat <- NULL
  # Interim setup
  N     <- length(n_seq)
  n_add <- sapply(1:N, function(a) length(unique(dat$id[dat$t0 <= accdat$tt[n_seq[a]]])))
  n_new <- diff(c(0, n_add))
  t_seq <- sapply(1:N, function(a) accdat$tt[n_seq[a]])
  idobs <- cbind(c(1, n_seq[-length(n_seq)] + 1), n_seq)
  idenr <- cbind(c(1, n_add[-length(n_add)] + 1), n_add)
  # Output storage
  trtlabs <- paste0("trt", 0:(K - 1))
  n_enr   <- matrix(0, N, K, dimnames = list(analysis = 1:N, treatment = trtlabs))
  n_obs   <- n_enr
  parlabs <- c("intercept", trtlabs)
  trt_mean <- matrix(0, N, K, dimnames = list(analysis = 1:N, treatment = trtlabs))
  trt_var  <- trt_mean
  eff_mean <- matrix(0, N, K - 1, dimnames = list(analysis = 1:N, treatment = trtlabs[-1]))
  eff_var  <- eff_mean
  b_mean   <- matrix(0, N, ncol(X), dimnames = list(analysis = 1:N, parameter = colnames(X)))
  p_supr   <- trt_mean
  i_supr   <- trt_mean
  i_infr   <- trt_mean
  p_eff    <- eff_mean
  i_eff    <- eff_mean
  i_inf    <- eff_mean
  i_acti   <- matrix(1, N+1, K, dimnames = list(analysis = 0:N, treatment = trtlabs))
  final    <- matrix(0, N, 1, dimnames = list(analysis = 1:N, final = 'final'))
  stopped <- FALSE
  for(i in 1:N) {
    final[i] <- stopped | i == N
    # Finish follow-up of enrolled participants if stopped
    if(stopped) {
      trtdat <- trtdat[trt != 99]
      obsdat <- dat[trtdat, on = .(id, trt)]
      n_enr[i, ] <- trtdat[, .N, keyby = trt][["N"]]
      n_obs[i, ] <- obsdat[t == 3, .N, keyby = trt][["N"]]
    } else { # Otherwise enrol new participants
      trtenr <- sample.int(K, n_new[i], TRUE, prob = alloc) - 1
      trtdat[between(id, idenr[i, 1], idenr[i, 2]), trt := trtenr]
      moddat <- dat[trtdat[between(id, 1, idenr[i, 2])], on = .(id, trt)]
      obsdat <- moddat[tt <= t_seq[i] & td > tt]
      n_enr[i, ] <- trtdat[between(id, 1, idenr[i, 2]), .N, keyby = trt][["N"]]
      n_obs[i, ] <- obsdat[t == 3, .N, keyby = trt][["N"]]
    }

    fit <- estimate_lmm_model(obsdat, zero_sum = zero_sum, scale = scale, ...)
    mu <- fit$mu
    Sigma <- fit$Sigma
    mu_t3 <- drop(X %*% mu)[10:13]
    Sigma_t3 <- (X %*% Sigma %*% t(X))[10:13, 10:13]

    # Update the inferences
    b_mean[i, ] <- mu
    trt_mean[i, ] <- mu_t3
    trt_var[i, ]  <- diag(Sigma_t3)
    eff_mean[i, ] <- drop(Z %*% X %*% mu)
    eff_var[i, ]  <- diag((Z %*% X) %*% Sigma %*% t(Z %*% X))
    # b_mean[i, ]   <- drop(C %*% mu[5:13])[9:12]
    # b_var[i, ]    <- diag(C %*% Sigma %*% t(C))

    means <- mvnfast::rmvn(2e4, mu_t3, Sigma_t3)

    # Conclusions carry forward
    # i.e. once decided superior/inferior, also superior/inferior at all future analyses
    if(i == 1) {
      p_eff[i, ]  <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
      i_eff[i, ]  <- p_eff[i, ] > eff_eps
      i_inf[i, ]  <- p_eff[i, ] < 1 - eff_eps
      p_supr[i, ] <- prob_supr(means)
      i_supr[i, ] <- p_supr[i, ] > sup_eps
      i_infr[i, ] <- p_supr[i, ] < (1 - sup_eps) / (K - 1)
      # If something superior, drop everything else
      if(any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1
      i_acti[i+1, ] <- 1 - i_infr[i, ]
    } else {
      p_eff[i, ]    <- 1 - pnorm(0, eff_mean[i, ], sqrt(eff_var[i, ]))
      # i_eff[i, ]    <- (p_eff[i, ] > eff_eps) | (i_eff[i-1,] == 1)
      i_eff[i, ]    <- p_eff[i, ] > eff_eps
      # i_inf[i, ]    <- (p_eff[i, ] < 1 - eff_eps) | (i_inf[i-1, ] == 1)
      i_inf[i, ]    <- p_eff[i, ] < 1 - eff_eps
      # Only include active in superiority assessment
      p_supr[i, i_acti[i, ] == 1] <- prob_supr(means[, i_acti[i, ] == 1, drop = F])
      i_supr[i, ]   <- (p_supr[i, ] > sup_eps) | (i_supr[i-1, ] == 1)
      i_infr[i, ]   <- p_supr[i, ] < min(1, (1 - sup_eps) / (sum(i_acti[i, ]) - 1))
      # If something superior, drop everything else
      if(any(i_supr[i, ] == 1)) i_infr[i, i_supr[i, ] != 1] <- 1
      i_acti[i+1, ] <- 1 - i_infr[i, ]
    }

    # Update allocations
    if(brar) {
      ratio <- sqrt(p_supr[i, ] * i_acti[i+1, ] / n_enr[i, ])
      alloc <- ratio / sum(ratio)
    } else {
      alloc <- (alloc * i_acti[i+1, ]) / sum(alloc * i_acti[i+1, ])
    }

    # If we stopped last analysis, we are done
    if(stopped) break
    # Should the next analysis be final?
    if(any(i_supr[i, ] == 1)) stopped <- TRUE
  }

  idx <- 1:i
  # Results
  out <- list(
    n_enr = n_enr[idx, , drop = F],
    n_obs = n_obs[idx, , drop = F],
    trt_mean = trt_mean[idx, , drop = F],
    trt_var  = trt_var[idx, , drop = F],
    eff_mean = eff_mean[idx, , drop = F],
    eff_var = eff_var[idx, , drop = F],
    # b_mean = b_mean[idx ,, drop = F],
    p_supr = p_supr[idx, , drop = F],
    p_eff = p_eff[idx, , drop = F],
    i_eff = i_eff[idx, , drop = F],
    i_inf = i_inf[idx, , drop = F],
    i_supr = i_supr[idx, , drop = F],
    i_infr = i_infr[idx, , drop = F],
    i_acti = i_acti[idx + 1, , drop = F]
  )
  if(make_dt) {
    trial_as_dt(out)
  } else {
    out
  }
}

#' Convert trial list into trial data.table
#'
#' @param res A list which is the result of simulate_longitudinal_trial
#' @return The list converted to a data.table with one row being an analysis by treatment result
#' @export
trial_as_dt <- function(res) {
  tmp <- lapply(
    1:length(res), function(a) {
      melt(as.data.table(res[[a]], keep.rownames = "analysis"),
           id.vars = "analysis",
           value.name = names(res)[a])
  })
  tmp <- Reduce(function(x, y) merge(x, y, on = .(analysis, variable), all = TRUE), tmp)
  setkey(tmp, analysis, variable)
  return(tmp)
}

