# load packages
require(Rcpp)

# load functions
source("functions/simulation_examples.R")
source("functions/r_wrappers.R")
Rcpp::sourceCpp("functions/sequential_tests.cpp")
 
# parameters
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
l <- ceiling(id / 1000)
nmax <- 10000

# simulation
set.seed(id)
out <- sim_all(nmax, l = l)
for (j in seq_along(out)) {
  df <- out[[j]]
  xgrid <- seq(min(df$x), max(df$x), 0.025)
  ygrid <- seq(min(df$y), max(df$y), 0.025)
  outputs <- list(
    simple = srt_simple(df$x, df$y),
    sinkhorn = srt_sinkhorn(df$x, df$y),
    bet = srt_bet(df$x, df$y),
    ks = seq_ks(df$x, df$y, xgrid, ygrid)
  )
  martingales_pm <- lapply(
    outputs,
    function(out) get_martingale(out, combination = "product_mean")
  )
  martingales_mp <- lapply(
    outputs,
    function(out) get_martingale(out, combination = "mean_product")
  )
  out[[j]] <- data.frame(
    id = id,
    l = l,
    method = names(outputs),
    sim = names(out)[j],
    pm = I(martingales_pm),
    mp = I(martingales_mp)
  )
}

# export results
save(list = "out", file = paste0("independence_", id, ".rda"))

