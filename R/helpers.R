#' Mass-weighted urn randomisation

#' Zhao W. (2015).
#' Mass weighted urn design--A new randomization
#' algorithm for unequal allocations.
#' Contemporary clinical trials, 43, 209–216.
#' https://doi.org/10.1016/j.cct.2015.06.008
#'
#' @param target_alloc The target allocation ratios
#' @param sample_size The number of allocations to generate
#' @param alpha Parameter to control imbalance between arms
#' @return A list detailing the mass-weighted-urn process.
#' @importFrom stats aggregate contr.treatment formula
#' @importFrom stats model.matrix pnorm rexp rnorm runif sd setNames time
mass_weighted_urn_design <- function(target_alloc,
                                     sample_size,
                                     alpha = 4) {
  arms <- length(target_alloc)
  prob_alloc <- target_alloc / sum(target_alloc)
  # Masses
  x <- matrix(0, sample_size + 1, arms)
  x[1, ] <- alpha * prob_alloc
  # Sample size
  n <- matrix(0, sample_size + 1, arms)
  # Random number
  y <- runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size + 1, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation Predictability
  g <- rep(0, sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)

  imbalance_cap <- sqrt(sum(((alpha - 1) * (1 - prob_alloc) + (arms - 1))^2))

  for (i in 2:(sample_size + 1)) {
    # Update allocation probabilities
    ptmp <- alpha * prob_alloc - n[i - 1, ] + sum(n[i - 1, ]) * prob_alloc
    p[i - 1, ] <- pmax(ptmp, 0)
    p[i - 1, ] <- p[i - 1, ] / sum(p[i - 1, ])
    trt[i - 1] <- findInterval(y[i - 1], c(0, cumsum(p[i - 1, ])))
    # Update sample sizes
    n[i, ] <- n[i - 1, ]
    n[i, trt[i - 1]] <- n[i, trt[i - 1]] + 1
    # Update urn masses
    x[i, trt[i - 1]] <- x[i - 1, trt[i - 1]] - 1 + prob_alloc[trt[i - 1]]
    x[i, -trt[i - 1]] <- x[i - 1, -trt[i - 1]] + prob_alloc[-trt[i - 1]]
    # Calculate imbalance
    d[i - 1] <- sqrt(sum((n[i, ] - (i - 1) * prob_alloc)^2))
    # Calculate allocation predictability
    g[i] <- d[i - 1] / alpha
  }
  return(list(
    max_imbalance_bound = imbalance_cap,
    imbalance = d,
    alloc_predict = g,
    rand_num = y,
    trt = trt,
    mass = x,
    sample_size = n,
    selection_prob = p
  ))
}

tnmean_inner <- function(a, b) {
  (dnorm(b) - dnorm(a)) / (pnorm(b) - pnorm(a))
}

#' Calculate mean of truncated normal
#'
#' @param a Lower bound
#' @param b Upper bound
#' @param mu Mean of normal
#' @param sigma Standard deviation of normal
#' @return Mean of truncated normal
#' @export
tnmean <- function(a, b, mu, sigma) {
  if(a > b) stop("a must be less than or equal to b")
  if(mu < a || mu > b) stop("mu outside boundary, algorithm unstable")
  if(b == a) return(b)
  alpha <- (a - mu) / sigma
  beta <- (b - mu) / sigma
  return(mu - tnmean_inner(alpha, beta) * sigma)
}


#' Calculate variance of truncated normal
#'
#' @param a Lower bound
#' @param b Upper bound
#' @param mu Mean of normal
#' @param sigma Standard deviation of normal
#' @return Variance of truncated normal
#' @export
tnvar <- function(a, b, mu, sigma) {
  if(a > b) stop("a must be less than or equal to b")
  if(mu < a || mu > b) stop("mu outside boundary, algorithm unstable")
  if(b == a) return(b)
  alpha <- (a - mu) / sigma
  beta <- (b - mu) / sigma
  tmp1 <- (beta * dnorm(beta) - alpha * dnorm(alpha)) / (pnorm(beta) - pnorm(alpha))
  tmp2 <- ((dnorm(beta) - dnorm(alpha)) / (pnorm(beta) - pnorm(alpha))) ^ 2
  return(sigma ^ 2 * (1 - tmp1 - tmp2))
}


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
  for (i in 1:K) {
    B <- A[, C[, i]]
    P[i] <- mvtnorm::pmvnorm(
      delta, Inf, drop(B %*% mu),
      sigma = B %*% (Sigma %*% t(B))
    )
  }
  return(P)
}

sim_data <- function(n, mu_x, mu_y, V_xy, p_z) {
  g <- length(p_z)
  z <- sample.int(g, n, replace = TRUE, prob = p_z)
  L <- t(chol(V_xy))
  mu <- rbind(mu_x[z], mu_y[z])
  xy <- t(mu + L %*% matrix(stats::rnorm(2 * n), 2, n))
  data.frame(z = factor(z), x = xy[, 1], y = xy[, 2])
}

#'
#'
#' @importFrom dplyr "%>%"
sim_data_ord <- function(n, mu_y, R, p_z, agg = FALSE) {
  g <- length(p_z)
  l <- nrow(R)
  N <- stats::rmultinom(1, n, p_z)[, 1]
  out <- lapply(1:g, function(a) {
    if (N[a] > 0) {
      tmp <- GenOrd::ordsample(
        N[a], rep(list(mu_y[a, ]), l), R,
        cormat = "continuous"
      ) - 1
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
  if (agg) {
    out <- out %>%
      dplyr::group_by(z, id) %>%
      dplyr::summarise(y = sum(y)) %>%
      dplyr::ungroup()
  }
  return(out[sample(nrow(out)), ])
}


sample_cond <- function(N, p_pref, p_int) {
  stopifnot(length(p_pref) == nrow(p_int))
  pref <- factor(
    sample.int(
      length(p_pref), N,
      prob = p_pref, replace = TRUE
    ), seq_along(p_pref)
  )
  int <- split(data.frame(pref), pref)
  dat <- unsplit(Map(function(data, r) {
    int <- factor(
      sample.int(
        ncol(p_int), nrow(data),
        replace = TRUE, prob = p_int[, r]
      ), seq_len(ncol(p_int))
    )
    data.frame(data, int)
  }, int, as.numeric(names(int))), pref)
  dat$com <- (as.numeric(dat$int) - 1) * ncol(p_int) + as.numeric(dat$pref)
  dat$trt <- factor(as.numeric(dat$int != 1))
  dat$gotpref <- factor(
    ifelse(dat$pref == 1, 1, ifelse(dat$int == dat$pref, 2, 3)), 1:3
  )
  return(dat)
}

sim_data_anova <- function(n, mu_x, mu_y, V_xy, p_pref, p_int) {
  pref_int <- sample_cond(n, p_pref, p_int)
  L <- t(chol(V_xy))
  mu <- rbind(mu_x[pref_int$com], mu_y[pref_int$com])
  xy <- t(mu + L %*% matrix(stats::rnorm(2 * n), 2, n))
  cbind.data.frame(pref_int, x = xy[, 1], y = xy[, 2], chg = xy[, 2] - xy[, 1])
}

brar <- function(pbest, variance, n, active) {
  m <- length(active) + 1
  p <- rep(1 / m, m)
  w <- sqrt(pbest * variance / (n + 1))
  w[!active] <- 0
  if (sum(w) > 0) {
    p[-1] <- (1 - p[1]) * w / sum(w)
  } else {
    p[-1] <- 0
    p[1] <- 1
  }
  return(p)
}

brar2 <- function(pbest, variance, n, active, fix_ctr = NULL) {
  m <- length(pbest)
  w <- sqrt(pbest / (n + 1))
  w[!active] <- 0
  # If we aren't fixing the control group allocation, just straight BRAR
  if (is.null(fix_ctr)) {
    p <- w / sum(w)
  } else { # If we are fixing...
    # Check that at least one non-control is active
    p <- rep(fix_ctr, m)
    if (sum(w[-1]) > 0) {
      p[-1] <- (1 - p[1]) * w[-1] / sum(w[-1])
    } else { # Otherwise, assign everyone to control
      p[-1] <- 0
      p[1] <- 1
    }
  }
  return(p)
}

apply_min_brar <- function(p, brar_min, a) {
  # stopifnot("sum(p) != 1" = sum(p) == 1)
  stopifnot("brar_min > 1 / sum(a)" = brar_min <= 1 / sum(a))
  stopifnot("p > 0 but a == 0" = all(a[p > 0] > 0))
  stopifnot(brar_min * sum(a[-1]) <= (1 - p[1]))
  id1 <- (p < brar_min) & a
  id1[1] <- FALSE
  id2 <- p >= brar_min
  id2[1] <- FALSE
  R <- sum(pmax(0, brar_min - p[-1]) * a[-1])
  p[id1] <- brar_min
  r <- (p[-1] - brar_min) * a[-1]
  p[-1] <- p[-1] - R * r / sum(r)
  return(p)
}

odds_trans_cdf <- function(cdf, theta) {
  return(cdf / (cdf + (1 - cdf) * exp(-theta)))
}

unmatrix <- function(m) {
  nam <- dimnames(m)
  res <- c(m)
  names(res) <- c(outer(paste0(names(nam)[1], nam[[1]]),
    paste0(names(nam)[2], nam[[2]]), paste,
    sep = "_"
  ))
  return(res)
}
