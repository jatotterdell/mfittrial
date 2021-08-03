# Suppose that we first choose a model, and then make decisions conditional on
# the selected model

# X2 is the more complex model compared to X1
choose_model <- function(X1, X2, y, M0, S0, A0, B0) {
  model1 <- varapproxr::vb_lm(X1, y, M0, S0, A0, B0)
  model2 <- varapproxr::vb_lm(X2, y, M0[1:ncol(X2)], S0[1:ncol(X2), 1:ncol(X2)], A0, B0)
  bf <- max(model1$elbo) - max(model2$elbo)
  if(bf > 0) {
    mod <- model1
  } else {
    mod <- model2
  }
  return(list(model = (bf < 0) + 1, res = mod))
}

#' @export
run_model_choice_trial <- function(
  n_seq = c(107, 214, 320),
  n_delay = 27,
  arms = 4,
  m_baseline = rep(0, 16),
  m_outcome = rep(0, 16),
  v_bo = diag(9^2, 2),
  p_pref = rep(1/ 4, 4),
  p_alloc = matrix(0.25, 4, 4, dimnames = list(int = 1:4, pref = 1:4)),
  kappa_sup = 0.95,
  brar = FALSE,
  M0 = rep(0, 16 + adjust_baseline),
  S0 = diag(c(100^2, rep(100^2, 15))),
  A0 = 1e-2,
  B0 = 1e-2,
  perm_drop = F,
  adjust_baseline = F
) {
  if(min(n_seq) <= n_delay) stop("min(n_seq) must exceed n_delay to be valid.")

  N <- max(n_seq)
  K <- length(n_seq)
  n_new <- diff(c(0, n_seq))
  n_ints <- K - 1
  active <- rep(1, arms - 1)
  stopped <- 0

  D <- sim_data_anova(n_new[1], m_baseline, m_outcome, v_bo, p_pref, p_alloc)
  X1 <- X_con[D$com, ]
  X2 <- X_con_red[D$com, ]
  P1 <- ncol(X1)
  P2 <- ncol(X2)

  p_best_sub <- matrix(0, K, 16,
                       dimnames = list("interim" = 1:K,
                                       "arm" = apply(expand.grid(paste0("int", 1:4),
                                                                 paste0("pref", 1:4)), 1,
                                                     paste, collapse = "")))
  p_best_int <- matrix(0, K, 4, dimnames = list("interim" = 1:K, "int" = 1:4))
  chosen_mod <- integer(K)
  model_fit <- vector("list", K)
  mean_mu <- matrix(0, K, nrow(X_con), dimnames = list("interim" = 1:K, "par" = rownames(X_con)))
  mean_beta <- matrix(0, K, nrow(Q), dimnames = list("interim" = 1:K, "par" = rownames(Q)))
  is_active_sub <- matrix(0, K, nrow(X_con), dimnames = list("interim" = 1:K, "par" = rownames(X_con)))

  # If interims cycle through them
  if(n_ints > 0) {
    n_interim <- matrix(0, 4, 4, dimnames = list(trt = 1:4, pref = 1:4))
    for(k in seq_len(n_ints)) {
      # Undertake analysis excluding the delayed outcomes
      Dint <- D[1:(n_seq[k] - n_delay), ]
      Xint1 <- X1[1:(n_seq[k] - n_delay), ]
      Xint2 <- X2[1:(n_seq[k] - n_delay), ]
      yint <- Dint[, "y"]
      # Full model
      choice <- choose_model(Xint1, Xint2, yint, M0, S0, A0, B0)
      chosen_mod[k] <- choice[[1]]
      res_vb <- choice[[2]]
      if(chosen_mod[k] == 1) {
        XX <- X_con
        QQ <- Q
        P <- P1
      } else {
        XX <- X_con_red
        QQ <- Q_red
        P <- P2
      }
      model_fit[[k]] <- res_vb
      draws <- mvnfast::rmvn(1e4, res_vb$mu[1:P, ], res_vb$Sigma[1:P, 1:P])
      mu <- draws %*% t(XX)
      beta <- draws %*% t(QQ)
      variance <- matrix(diag(XX %*% res_vb$Sigma[1:P, 1:P] %*% t(XX)), 4, 4, byrow = T)

      # Summarise posterior probabilities of interest
      mean_mu[k, ] <- colMeans(mu)
      mean_beta[k, 1:nrow(QQ)] <- colMeans(beta)
      n_sub <- table(D[1:(n_seq[k] - n_delay), ]$int, D[1:(n_seq[k] - n_delay), ]$pref)
      p_best_mat <- sapply(1:4, function(i) automaticsims:::prob_max(mu[, c(1,5,9,13) + (i - 1)]))
      p_best_sub[k, ] <- c(sapply(1:4, function(i) automaticsims:::prob_max(mu[, c(1,5,9,13) + (i - 1)])))
      p_best_int[k, ] <- automaticsims:::prob_max(beta[, 2:5])
      is_sup_sub_mat <- sweep(p_best_mat, 2, kappa_sup, ">")

      # Tweak thresholds, if only one arm active, always call it superior
      act <- colSums(is_active_sub)
      kappa_act <- ifelse(act > 1, (1 - kappa_sup) / (act - 1), 0)
      kappa_sup <- ifelse(act > 1, kappa_sup, 0)
      if(brar) {
        if(perm_drop) {
          is_active_sub <- is_active_sub * sweep(p_best_mat, 2, kappa_act, ">")
        } else {
          is_active_sub <- sweep(p_best_mat, 2, kappa_act, ">")
        }
        p_alloc[] <- sapply(1:4, function(i) {
          brar2(p_best_mat[, i], 1, n_sub[, i], is_active_sub[, i])
        })
      }

      # If every subgroup has a superior treatment then stop
      if(all(apply(is_sup_sub_mat, 2, any))) {
        stopped <- 1
      } else {
        # Generate next lot of data
        D <- rbind.data.frame(D, sim_data_anova(n_new[k + 1], m_baseline, m_outcome, v_bo, p_pref, p_alloc))
        X1 <- X_con[D$com, ]
        X2 <- X_con_red[D$com, ]
      }
      if(stopped == 1) break
    }
  }

  # Undertake final analysis
  y <- D$y
  choice <- choose_model(X1, X2, y, M0, S0, A0, B0)
  chosen_mod[K] <- choice[[1]]
  res_vb <- choice[[2]]
  if(chosen_mod[K] == 1) {
    XX <- X_con
    QQ <- Q
    P <- P1
  } else {
    XX <- X_con_red
    QQ <- Q_red
    P <- P2
  }
  model_fit[[K]] <- res_vb
  draws <- mvnfast::rmvn(1e4, res_vb$mu[1:P, ], res_vb$Sigma[1:P, 1:P])
  mu <- draws %*% t(XX)
  beta <- draws %*% t(QQ)
  variance <- matrix(diag(XX %*% res_vb$Sigma[1:P, 1:P] %*% t(XX)), 4, 4, byrow = T)
  # Summarise posterior probabilities of interest
  mean_mu[k, ] <- colMeans(mu)
  mean_beta[k, 1:nrow(QQ)] <- colMeans(beta)
  n_sub <- table(D$int, D$pref)
  p_best_mat <- sapply(1:4, function(i) automaticsims:::prob_max(mu[, c(1,5,9,13) + (i - 1)]))
  p_best_sub[K, ] <- c(sapply(1:4, function(i) automaticsims:::prob_max(mu[, c(1,5,9,13) + (i - 1)])))
  p_best_int[K, ] <- automaticsims:::prob_max(beta[, 2:5])
  is_sup_sub_mat <- sweep(p_best_mat, 2, kappa_sup, ">")

  return(tibble::enframe(list(
    model = tibble::tibble(chosen_mod) %>% tibble::rownames_to_column("analysis"),
    model_fit = tibble::enframe(model_fit, "analysis", "model_fit"),
    mean_mu = tidyr::gather(tibble::as_tibble(mean_mu, rownames = "analysis"), "parameter", mean_mu, -analysis),
    mean_beta = tidyr::gather(tibble::as_tibble(mean_beta, rownames = "analysis"), "parameter", mean_beta, -analysis),
    p_best_sub = tidyr::gather(tibble::as_tibble(p_best_sub, rownames = "analysis"), "arm", p_best_sub, -analysis),
    p_best_int = tidyr::gather(tibble::as_tibble(p_best_int, rownames = "analysis"), "intervention", p_best_int, -analysis)
  )))
}