get_ord_dist <- function() {
    x <- 0:52
    p <- c(
        seq(0.001, 0.0095, length.out = 15),
        seq(0.01, 0.0245, length.out = 10),
        seq(0.025, 0.0295, length.out = 15),
        seq(0.03, 0.015, length.out = 13)
    )
    p <- p / sum(p)
    return(cbind(x = x, p = p))
}
