# packages
library(tidyverse)
library(ggthemes)

# functions
source("functions/simulation_examples.R")

# plot options
theme_set(theme_bw(base_size = 12))

# colors
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# parameters
n <- 100
l1 <- 1
l2 <- 9

# generate data
set.seed(123)
sims <- sim_all(n, l1)
sims <- sims[names(sims) != "null"]
sims <- mapply(
  function(data, name) {
    cbind(sim = name, data)
  },
  data = sims,
  name = names(sims),
  SIMPLIFY = FALSE
)
sims1 <- do.call(rbind, sims)
sims1$l <- l1
sims <- sim_all(n, l2)
sims <- sims[names(sims) != "null"]
sims <- mapply(
  function(data, name) {
    cbind(sim = name, data)
  },
  data = sims,
  name = names(sims),
  SIMPLIFY = FALSE
)
sims2 <- do.call(rbind, sims)
sims2$l <- l2
sims <- rbind(sims1, sims2)
rownames(sims) <- NULL

# plot

## for main part of paper
sim_illustration_main <- sims %>%
  group_by(sim, l) %>%
  mutate(
    x = sapply(seq_len(n), function(k) rank(x[seq_len(k)])[k] / k - runif(1) / k),
    y = sapply(seq_len(n), function(k) rank(y[seq_len(k)])[k] / k - runif(1) / k),
    l = factor(l)
  ) %>%
  filter(sim %in% c("linear", "local", "circular")) %>%
  mutate(sim = str_to_title(sim)) %>%
  ggplot() +
  geom_point(aes(x = x, y = y, shape = l, color = l, alpha = l)) +
  scale_color_manual(values = colpal[c(2, 7)]) +
  scale_alpha_manual(values = c(1, 0.4)) +
  facet_grid(cols = vars(sim)) +
  labs(x = expression(R[n]), y = expression(S[n])) +
  theme(legend.position = "none")


## for appendix
sim_illustration_suppl <- sims %>%
  group_by(sim, l) %>%
  mutate(
    x = sapply(seq_len(n), function(k) rank(x[seq_len(k)])[k] / k - runif(1) / k),
    y = sapply(seq_len(n), function(k) rank(y[seq_len(k)])[k] / k - runif(1) / k),
    l = factor(l)
  ) %>%
  filter(!sim %in% c("linear", "local", "circular")) %>%
  mutate(sim = str_to_title(sim)) %>%
  ggplot() +
  geom_point(aes(x = x, y = y, shape = l, color = l, alpha = l)) +
  scale_color_manual(values = colpal[c(2, 7)]) +
  scale_alpha_manual(values = c(1, 0.4)) +
  facet_grid(cols = vars(sim)) +
  labs(x = expression(R[n]), y = expression(S[n])) +
  theme(legend.position = "none")

pdf(file = "sim_illustration_main_revised.pdf", width = 8, height = 3)
sim_illustration_main
dev.off()

pdf(file = "sim_illustration_suppl_revised.pdf", width = 8, height = 3)
sim_illustration_suppl
dev.off()


