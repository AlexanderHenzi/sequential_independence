# load functions
source("simulation_examples.R")

# parameters
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
l <- ceiling(id / 1000)
nmax <- 10000

# simulation
set.seed(id)
dfs <- sim_all(nmax, l = l)
dfs <- do.call(
  rbind,
  mapply(function(d, n) cbind(d, sim = n), dfs, names(dfs), SIMPLIFY = FALSE)
)

# export data to read it in python
write.table(
  x = dfs,
  file = paste0("simulation_examples_", id, ".csv"),
  row.names = FALSE,
  sep = ","
)