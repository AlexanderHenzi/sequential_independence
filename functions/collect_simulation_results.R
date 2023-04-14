# packages
library(tidyverse)

# functions
Rcpp::sourceCpp("functions/sequential_tests.cpp")

# parameters
M <- 10000
alpha <- c(1e-4, 1e-3, 1e-2, 0.05)

# get results for stopping times
pb <- txtProgressBar(max = M)
results <- vector("list", M)
for (i in seq_len(M)) {
  setTxtProgressBar(pb, i)
  load(paste0("independence_", i, ".rda"))
  tmp <- do.call(rbind, out) %>%
    unnest(cols = c(pm, mp)) %>%
    group_by(method, sim, l) %>%
    summarise(
      rej = list(tibble(
        alpha = alpha,
        pm = sapply(alpha, function(a) get_first_rejection(pm, 1/a)),
        mp = sapply(alpha, function(a) get_first_rejection(mp, 1/a)),
      )),
      .groups = "drop"
    ) %>%
    unnest(rej)
  results[[i]] <- tmp
}
close(pb)
results <- do.call(rbind, results)
save(list = "results", file = "simulation_results.rda")

# get results for comparison to non-sequential BET
ns <- c(128, 256, 512)
tt <- c(9.17, 13.8, 16.9)

pb <- txtProgressBar(max = M)
results <- vector("list", M)
for (i in seq_len(M)) {
  setTxtProgressBar(pb, i)
  load(paste0("independence_", i, ".rda"))
  tmp <- do.call(rbind, out) %>%
    unnest(cols = c(pm, mp)) %>%
    filter(method == "sinkhorn") %>%
    select(id, l, sim, pm) %>%
    group_by(sim, l) %>%
    summarise(
      rej = list(tibble(
        n = ns,
        rej = mapply(
          function(n, t) get_first_rejection(pm[seq_len(n)], t),
          n = ns,
          t = tt
        )
      )),
      .groups = "drop"
    ) %>%
    unnest(rej)
  results[[i]] <- tmp
}
close(pb)
results <- do.call(rbind, results)
save(list = "results", file = "simulations_comparison_bet.rda")
