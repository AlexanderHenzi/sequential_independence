# packages
library(tidyverse)

# functions
source("functions/simulation_examples.R")

# plot options
theme_set(theme_bw(base_size = 12))

# parameters
n <- 100
l <- 1

# generate data
set.seed(123)
sims <- sim_all(n, l)
sims <- sims[names(sims) != "null"]
sims <- mapply(
  function(data, name) {
    cbind(sim = name, data)
  },
  data = sims,
  name = names(sims),
  SIMPLIFY = FALSE
)
sims <- do.call(rbind, sims)
rownames(sims) <- NULL

# plot

## for main part of paper
sim_illustration_main <- sims %>%
  group_by(sim) %>%
  mutate(x = ecdf(x)(x), y = ecdf(y)(y)) %>%
  filter(sim %in% c("linear", "local", "circular")) %>%
  mutate(sim = str_to_title(sim)) %>%
  ggplot() +
  geom_point(aes(x = x, y = y)) +
  facet_grid(cols = vars(sim)) +
  labs(x = "X", y = "Y")


## for appendix
sim_illustration_suppl <- sims %>%
  group_by(sim) %>%
  mutate(x = ecdf(x)(x), y = ecdf(y)(y)) %>%
  filter(!sim %in% c("linear", "local", "circular")) %>%
  mutate(sim = str_to_title(sim)) %>%
  ggplot() +
  geom_point(aes(x = x, y = y)) +
  facet_grid(cols = vars(sim)) +
  labs(x = "X", y = "Y")

pdf(file = "sim_illustration_main.pdf", width = 8, height = 3)
sim_illustration_main
dev.off()

pdf(file = "sim_illustration_suppl.pdf", width = 8, height = 3)
sim_illustration_suppl
dev.off()


