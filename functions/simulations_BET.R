# load packages
require(Rcpp)

# load functions
source("functions/simulation_examples.R")
source("~functions/BET_functions_JASA0223.r")

# parameters
ns <- c(64, 128, 256, 512)
alpha <- 0.05
nsim <- 10000
nmax <- 10000
results <- vector("list", 10000)

pb <- txtProgressBar(max = nsim)
for (id in seq_len(nsim)) {
  l <- ceiling(id / 1000)
  
  # simulation
  set.seed(id)
  out <- sim_all(nmax, l = l)
  for (j in seq_along(out)) {
    df <- out[[j]]
    out[[j]] <- data.frame(
      n = ns,
      reject = sapply(
        ns,
        function(n) 
          as.integer(BETs(
            ecdf(df$x[seq_len(n)])(df$x[seq_len(n)]),
            ecdf(df$y[seq_len(n)])(df$y[seq_len(n)]),
          )$bet.s.pvalue <= alpha)
      ),
      sim = names(out)[j],
      id = id,
      l = l
    )
  }
  results[[id]] <- do.call(rbind, out)
  setTxtProgressBar(pb, id)
}
close(pb)

results <- do.call(rbind, results)
save(list = "results", file = "simulations_bet_results.rda")
