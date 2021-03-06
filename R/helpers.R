#' Creates C matrix for use in mvnorm_prob_each_best
#'
#' @param K The dimension of C
#' @return A matrix
create_C_mat <- function(K) {
  cbind(1:K, sapply(2:K, function(i) c(2:i, c(1, setdiff(2:K, 2:i)))))
}

mvnorm_prob_each_best <- function(mu, Sigma, delta = 0) {
  K <- length(mu)
  C <- create_C_mat(K)
  A <- cbind(1, diag(-1, K - 1))
  P <- as.numeric(K)
  for(i in 1:K) {
    B <- A[, C[, i]]
    P[i] <- mvtnorm::pmvnorm(delta, Inf, drop(B %*% mu), sigma = B %*% (Sigma %*% t(B)))
  }
  return(P)
}

sim_data <- function(
  n, mu_x, mu_y, V_xy, p_z
) {
  g <- length(p_z)
  z <- sample.int(g, n, replace = TRUE, prob = p_z)
  L <- t(chol(V_xy))
  mu <- rbind(mu_x[z], mu_y[z])
  xy <- t(mu + L %*% matrix(stats::rnorm(2*n), 2, n))
  data.frame(z = factor(z), x = xy[, 1], y = xy[, 2])
}

#'
#'
#' @importFrom dplyr "%>%"
sim_data_ord <- function(
  n, mu_y, R, p_z, agg = FALSE
) {
  g <- length(p_z)
  l <- nrow(R)
  N <- stats::rmultinom(1, n, p_z)[, 1]
  out <- lapply(1:g, function(a) {
    if(N[a] > 0) {
      tmp <- GenOrd::ordsample(N[a], rep(list(mu_y[a, ]), l), R, cormat = "continuous") - 1
      colnames(tmp) <- paste0("i", 1:l)
      return(tibble::as_tibble(tmp))
    } else {
      return(NULL)
    }
  })
  names(out) <- 1:g
  out <- dplyr::bind_rows(out, .id = "z")
  out <- out %>%
    dplyr::mutate(id = 1:dplyr::n(), z = factor(z, levels = 1:g)) %>%
    tidyr::gather(item, y, -id, -z)
  if(agg) {
    out <- out %>%
      dplyr::group_by(z, id) %>%
      dplyr::summarise(y = sum(y)) %>%
      dplyr::ungroup()
  }
  return(out[sample(nrow(out)), ])
}


sample_cond <- function(N, p_pref, p_int) {
  stopifnot(length(p_pref) == nrow(p_int))
  pref <- factor(sample.int(length(p_pref), N, prob = p_pref, replace=TRUE), 1:length(p_pref))
  int <- split(data.frame(pref), pref)
  dat <- unsplit(Map(function(data, r) {
    int <- factor(sample.int(ncol(p_int), nrow(data), replace = TRUE, prob = p_int[, r]), 1:ncol(p_int))
    data.frame(data, int)
  }, int, as.numeric(names(int))), pref)
  dat$com <- (as.numeric(dat$int) - 1) * ncol(p_int) + as.numeric(dat$pref)
  dat$trt <- factor(as.numeric(dat$int != 1))
  dat$gotpref <- factor(ifelse(dat$pref == 1, 1, ifelse(dat$int == dat$pref, 2, 3)), 1:3)
  # dat$gotpref <- factor(as.numeric((dat$int == dat$pref) & dat$pref != 1))
  # dat$havepref <- as.numeric(dat$pref != 1)
  return(dat)
}


sim_data_anova <- function(
  n, mu_x, mu_y, V_xy, p_pref, p_int
) {
  pref_int <- sample_cond(n, p_pref, p_int)
  L <- t(chol(V_xy))
  mu <- rbind(mu_x[pref_int$com], mu_y[pref_int$com])
  xy <- t(mu + L %*% matrix(stats::rnorm(2*n), 2, n))
  cbind.data.frame(pref_int, x = xy[, 1], y = xy[, 2], chg = xy[, 2] - xy[, 1])
}


brar <- function(pbest, variance, n, active) {
  m <- length(active) + 1
  p <- rep(1 / m, m)
  w <- sqrt(pbest * variance / (n + 1))
  w[!active] <- 0
  if(sum(w) > 0) {
    p[-1] <- (1 - p[1]) * w / sum(w)
  } else {
    p[-1] <- 0
    p[1] <- 1
  }
  return(p)
}


brar2 <- function(pbest, variance, n, active, fix_ctr = NULL) {
  m <- length(pbest)
  # w <- sqrt(pbest * variance / (n + 1))
  w <- sqrt(pbest / (n + 1))
  w[!active] <- 0
  # If we aren't fixing the control group allocation, just straight BRAR
  if(is.null(fix_ctr)) {
    # w[!active] <- 0
    p <- w / sum(w)
  } else { # If we are fixing...
    # Check that at least one non-control is active
    p <- rep(fix_ctr, m)
    if(sum(w[-1]) > 0) {
      p[-1] <- (1 - p[1]) * w[-1] / sum(w[-1])
    } else { # Otherwise, assign everyone to control
      p[-1] <- 0
      p[1] <- 1
    }
  }
  return(p)
}


odds_trans_cdf <- function(cdf, theta) {
  return(cdf / (cdf + (1 - cdf) * exp(-theta)))
}

unmatrix <- function(m) {
  nam <- dimnames(m)
  res <- c(m)
  names(res) <- c(outer(paste0(names(nam)[1], nam[[1]]),
                        paste0(names(nam)[2], nam[[2]]), paste, sep = "_"))
  return(res)
}
