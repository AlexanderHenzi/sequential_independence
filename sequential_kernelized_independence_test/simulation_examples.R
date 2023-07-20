sim_null <- function(n, l) {
  x <- runif(n)
  y <- runif(n)
  data.frame(x = x, y = y)
}

sim_linear <- function(n, l) {
  x <- runif(n)
  y <- x + 6 * rnorm(n, sd = l/40)
  data.frame(x = x, y = y)
}

sim_parabolic <- function(n, l) {
  x <- runif(n)
  y <- (x - 0.5)^2 + 1.5 * rnorm(n, sd = l/40)
  data.frame(x = x, y = y)
}

sim_circular <- function(n, l) {
  theta <- runif(n, -2*pi, 2*pi)
  x <- cos(theta) + 2.5 * rnorm(n, sd = l/40)
  y <- sin(theta) + 2.5 * rnorm(n, sd = l/40)
  data.frame(x = x, y = y)
}

sim_sine <- function(n, l) {
  x <- runif(n)
  y <- sin(4*pi*x) + 8 * rnorm(n, sd = l/40)
  data.frame(x = x, y = y)
}

sim_checkerboard <- function(n, l) {
  w <- sample(x = seq_len(3), size = n, replace = TRUE)
  x <- w + rnorm(n, sd = l/40)
  v1 <- sample(x = c(2, 4), size = n, replace = TRUE) + 4 * rnorm(n, sd = l/40)
  v2 <- sample(x = c(1, 3, 5), size = n, replace = TRUE) + 4 * rnorm(n, sd = l/40)
  y <- v2
  y[w == 2] <- v1[w == 2]
  data.frame(x = x, y = y)
}

sim_local <- function(n, l) {
  g1 <- rnorm(n, sd = 0.5)
  g2 <- rnorm(n, sd = 0.5)
  x <- g1
  y <- numeric(n)
  y <- x + rnorm(n, sd = l/40)
  ind <- (g1 < 0 | g1 > 1 | g2 < 0 | g2 > 1)
  y[ind] <- g2[ind]
  data.frame(x = x, y = y)
}

sim_all <- function(n, l) {
  list(
    linear = sim_linear(n, l),
    parabolic = sim_parabolic(n, l),
    circular = sim_circular(n, l),
    sine = sim_sine(n, l),
    checkerboard = sim_checkerboard(n, l),
    local = sim_local(n, l),
    null = sim_null(n, l)
  )
}
