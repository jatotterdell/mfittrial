# Two models
# - one where preference is ignored
# - one with full preference interaction weighted by preference distribution

#' Generate potential outcomes if a participant were to receive each treatment
#' @param n The number of participants
#' @param m The mean response on each treatment
#' @param s The variability of response on each treatment
generate_data <- function(n, m, s) {
  y <- matrix(rnorm(length(m)*n, m, s), n, length(m), byrow = T)
  data.frame(
    id = 1:n,
    y  = y
  )
}


