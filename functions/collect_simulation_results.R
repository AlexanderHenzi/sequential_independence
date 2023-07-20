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
  
  # for methods implemented in R
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
  
  # for kernelized test
  out1 <- read.csv(
    paste0("simulations_aGRAPA_20_", i, ".csv"),
    sep = ",",
    header = TRUE
  )
  out2 <- read.csv(
    paste0("simulations_ONS_20_", i, ".csv"),
    sep = ",",
    header = TRUE
  )
  
  # find rejection times; multiply by 2 because observations are in pairs
  out <- rbind(out1, out2) %>%
    group_by(method, sim, l) %>%
    summarise(
      rej = list(tibble(
        alpha = alpha,
        pm = sapply(alpha, function(a) get_first_rejection(martingale, 1/a)) * 2,
        mp = sapply(alpha, function(a) get_first_rejection(martingale, 1/a)) * 2,
      )),
      .groups = "drop"
    ) %>%
    unnest(rej)
  results[[i]] <- rbind(out, tmp)
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
