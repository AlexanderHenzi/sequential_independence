#' Sequential rank test for independence
#' 
#' @param X X-variable (numeric vector with pairwise distinct entries).
#' @param Y Y variable (numeric vector with pairwise distinct entries).
#' @param d grid sizes (cell widths are 1/d).
#' @param n0 number of prior observations to include in each cell (reasonable
#'     choices are 1 or 1/2).
#'
#' @return 
#' A list of length 2*length(d), containing the multiplicative increments of the
#' test martingales in the first length(d) components, and matrices of estimated
#' cell probabilities in the last length(d) components (the [k,l]th entry 
#' corresponds to the estimated probability for the cell 
#' [(k-1)/d, k/d]x[(l-1)/d, l/d], k,l=1,...,d-1).
srt_simple <- function(X, Y, d = c(2, 4, 8, 16), n0 = 1) {
  out <- sequential_histogram(X = X, Y = Y, d = d, n0 = n0)
  structure(out, type = "srt_simple", nd = length(d))
}

#' Sequential rank test for independence, with uniform marginals, MLE
#' 
#' The parameters and output are as for srt_simple, except for the following 
#' two.
#' 
#' @param dmax test martingales will be computed for grid with cells that have
#'     widths 2^(-1),...,2^(-dmax).
#' @param control list of control parameters for gradient descent when 
#'     uniform_marginals is TRUE. Has no influence if uniform_marginals is
#'     FALSE. Gradient descent is stopped when the relative decrease in 
#'     log-likelihood between two iterations is less than eps; delta is a 
#'     tolerance to determine if there is an improvement in the step size 
#'     correction; maxit is the maximal number of iterations.
srt_constrained <- function(
    X,
    Y,
    dmax = 4,
    n0 = 1,
    control = list(eps = 1e-2, delta = 1e-5, maxit = 100)) {
  H <- matrix(1, nrow = 1, ncol = 1)
  h <- matrix(nrow = 2, ncol = 2, c(1, 1, 1, -1))
  basis_mat <- vector("list", dmax)
  c0 <- numeric(dmax)
  for (j in seq_len(dmax)) {
    d <- 2^j
    H <- (H %x% h) %x% h
    drop_cols <- c(seq_len(d), seq(d + 1, d^2, d))
    basis_mat[[j]] <- H[, -drop_cols, drop = FALSE]
    c0[j] <- 1/d^2
  }
  out <- sequential_histogram_constrained(
    X = X,
    Y = Y,
    basis_mat = basis_mat,
    c0 = c0,
    dmax = dmax,
    n0 = n0,
    eps = control$eps,
    delta = control$delta,
    maxit = control$maxit
  )
  structure(out, type = "srt_constrained", nd = dmax)
}

#' Sequential rank test for independence, with uniform marginals, Sinkhorn
#' 
#' The parameters and output are as for srt_simple, except for the following.
#'
#' @param control list of control parameters for row and column normalization.
#'    The tolerance criterion is eps, a number > 1. The algorithm stops if
#'    the column and row sums are all in the interval (1/eps, eps), or if
#'    maxit iterations are reached.
srt_sinkhorn <- function(
    X,
    Y,
    d = c(2, 4, 8, 16),
    n0 = 1,
    control = list(eps = 1.001, maxit = 20)) {
  out <- sequential_histogram_sinkhorn(
    X = X,
    Y = Y,
    d = d,
    n0 = n0,
    eps = control$eps,
    maxit = control$maxit
  )
  structure(out, type = "srt_sinkhorn", nd = length(d))
}

#' Sequential rank test for independence, uniform marginals, not derandomized
#' 
#' The parameters and output are as for srt_sinkhorn. This test generates
#' runif() random variables to really randomize the ranks, and does not do any
#' derandomization
srt_sinkhorn_random <- function(
    X,
    Y,
    d = c(2, 4, 8, 16),
    n0 = 1,
    control = list(eps = 1.001, maxit = 20)) {
  epsX <- diff(sort(X))
  epsX <- min(epsX[epsX > 0]) / 4
  epsY <- diff(sort(Y))
  epsY <- min(epsY[epsY > 0]) / 4
  n <- length(X)
  X <- X + runif(n, -epsX, epsX)
  Y <- Y + runif(n, -epsY, epsY)
  out <- sequential_histogram_sinkhorn_rand(
    X = X,
    Y = Y,
    d = d,
    n0 = n0,
    eps = control$eps,
    maxit = control$maxit
  )
  structure(out, type = "srt_sinkhorn_random", nd = length(d))
}

#' Sequential version of the BET
#' 
#' The parameters and output are as for srt_constrained.
#'
#' @return 
#' A list containing the multiplicative increments of the test martingales
#' (there are (2^(dmax) - 1)^2 of them) and the depth of each interaction 
#' (integer vector, one entry for each test martingale).
srt_bet <- function(
    X,
    Y,
    dmax = 4,
    n0 = 1) {
  zero_one <- vector("list", (2^dmax - 1)^2)
  lwr <- 1
  H <- matrix(1, nrow = 1, ncol = 1)
  h <- matrix(nrow = 2, ncol = 2, c(1, 1, 1, -1))
  for (j in seq_len(dmax)) {
    d <- 2^j
    H <- (H %x% h) %x% h
    drop_cols <- unique(c(
      c(seq_len(d), seq(d + 1, d^2, d)),
      seq(1, d - 1, 2) + rep((2 * d) * (0:(d/2 - 1)), each = d/2)
    ))
    basis_mat <- H[, -drop_cols, drop = FALSE]
    zero_one_j <- lapply(
      asplit(basis_mat, 2),
      function(x) matrix(nrow = d, ncol = d, replace(x, x == -1, 0))
    )
    zero_one[lwr:(lwr + length(zero_one_j) - 1)] <- zero_one_j
    lwr <- lwr + length(zero_one_j)
  }
  ls <- lengths(zero_one)
  d_groups <- c(0, c(which(ls[-1] > ls[-length(ls)]), length(ls)))
  out <- sequential_bet(
    X = X,
    Y = Y,
    n0 = n0,
    zero_one = zero_one,
    d_groups = d_groups
  )
  out <- lapply(asplit(out, 2), c)
  out <- c(out, list(log2(ls) / 2))
  structure(out, type = "srt_bet")
}

#' Sequential Kolmogorov-Smirnov test for independence
#' 
#' @param X X variable (numeric vector with pairwise distinct entries).
#' @param Y Y variable (numeric vector with pairwise distinct entries).
#' @param xgrid grid for Y variable (sorted numeric vector).
#' @param ygrid grid for Y variable (sorted numeric vector).
#' 
#' @return 
#' 
seq_ks <- function(X, Y, xgrid, ygrid) {
  rx <- range(X)
  ry <- range(Y)
  nxg <- length(xgrid)
  nyg <- length(ygrid)
  xgrid <- c(min(xgrid[1], rx[1]) - 1, xgrid, max(xgrid[nxg], rx[2]) + 1)
  ygrid <- c(min(ygrid[1], ry[1]) - 1, ygrid, max(ygrid[nyg], ry[2]) + 1)
  out <- sequential_ks(X, Y, xgrid, ygrid)
  out[[3]] <- out[[3]][-c(1, nxg + 2), -c(1, nyg + 2)]
  structure(out, type = "sequential_ks")
}

#' Construct test martingale from function output
#'
#' Takes the output from one of the functions above and combines it into a 
#' single test martingale.
#'
#' @param output output of one of the functions above in this file.
#' @param combination combination method for the methods based on binning. There
#'     are three options: "mean_product" takes the average over the test 
#'     martingales for the different grid sizes (eta = 1 in the article).
#'     "product_mean" (eta = 0 in the article) first averages the factors of
#'     each martingale and then computes the product. A constant martingale (1)
#'     is included.
#'
#' @return 
#' A numeric vector, the test martingale.
get_martingale <- function(
    output,
    combination = "mean_product") {
  type <- attributes(output)$type
  if (type == "sequential_ks") {
    out <- rep(0.5 * cumprod(output[[1]]) + 0.5 * cumprod(output[[2]]), each = 2)
  } else if (type == "srt_bet") {
    groups <- output[[length(output)]]
    ngroups <- groups[length(groups)]
    w <- 1 / (length(output) - 1)
    out <- cumprod(output[[1]]) * w
    if (length(output) > 1) {
      for (j in 2:(length(output) - 1)) out <- out + cumprod(output[[j]]) * w
    }
  } else {
    nd <- attributes(output)$nd
    if (identical(combination, "mean_product")) {
      w0 <- 1 / nd
      out <- cumprod(output[[1]]) * w0
      if (nd > 1) {
        for (j in 2:nd) out <- out + cumprod(output[[j]]) * w0
      }
    } else if (identical(combination, "product_mean")) {
      w0 <- 1 / (nd + 1)
      out <- rep(1, length(output[[1]])) * w0
      for (j in seq_len(nd)) out <- out + output[[j]] * w0
      out <- cumprod(out)
    }
  }
  out
}