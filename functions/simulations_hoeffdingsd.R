# load packages
require(independence)

# load functions
source("functions/simulation_examples.R")

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
          as.integer(hoeffding.D.test(
            df$x[seq_len(n)],
            df$y[seq_len(n)],
            precision = 1e-4
          )$p.value <= alpha)
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
save(list = "results", file = "simulations_hoeffdingsd_results.rda")
