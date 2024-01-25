# load packages
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(Rcpp)

# load functions
source("functions/r_wrappers.R")
source("functions/simulation_examples.R")
sourceCpp("functions/sequential_tests.cpp")

# ggplot2 settings
theme_set(theme_bw(base_size = 10))
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# parameters
theme_set(theme_bw(base_size = 12))
nsim <- 100000
d <- c(2, 4, 8, 16)
ns <- c(128, 256, 512, 1024, 2048, 4096)
nn <- length(ns)
nmax <- max(ns)
alpha <- 0.05

# run simulations
sinkhorn <- matrix(nrow = nsim, ncol = nn)
pb <- txtProgressBar(max = nsim)
for (i in seq_len(nsim)) {
  set.seed(i)
  x <- runif(nmax)
  y <- runif(nmax)
  m_sinkhorn <- get_martingale(srt_sinkhorn(x, y, d), "product_mean")
  sinkhorn[i, ] <- 1 / sapply(ns, function(k) max(m_sinkhorn[seq_len(k)]))
  setTxtProgressBar(pb, i)
}
close(pb)

# export results
save(list = "sinkhorn", file = "ville_gap.rda")

# plot
load(file = "ville_gap.rda")
colnames(sinkhorn) <- ns
data <- as_tibble(sinkhorn) %>%
  gather(key = "n", value = "p") %>%
  mutate(n = factor(n, levels = sort(ns), labels = sort(ns), ordered = TRUE))

grd <- seq(0, 1, 0.001)
bounds <- data %>%
  group_by(n) %>%
  summarise(
    bnds = list(tibble(
      grd = grd,
      cdf = ecdf(p)(grd)
    ))
  ) %>%
  unnest(cols = bnds) %>%
  mutate(n = factor(n, levels = sort(ns), labels = sort(ns), ordered = TRUE))

ville_illustration <- ggplot() +
  geom_abline(intercept = 0, slope = 1, lty = 5, color = "darkgray") +
  stat_ecdf(data = data, aes(x = p, group = n, color = n)) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.1)) +
  scale_color_manual(values = colpal[1:6]) +
  scale_fill_manual(values = colpal[1:6]) +
  labs(
    x = expression(p[N]),
    y = "Distribution function",
    color = "N",
    fill = "N"
  ) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))

pdf(file = "ville_illustration.pdf", width = 8, height = 3)
print(ville_illustration)
dev.off()

# numbers for text in article
freqs <- data %>%
  group_by(n) %>%
  summarise(p005 = mean(p <= 0.05), t005 = 1/quantile(p, 0.05, type = 1))
freqs
freqs$p005

# for illustration on real data
size <- 128
sinkhorn <- numeric(nsim)
pb <- txtProgressBar(max = nsim)
for (i in seq_len(nsim)) {
  set.seed(i)
  x <- runif(size)
  y <- runif(size)
  m_sinkhorn <- get_martingale(
    srt_sinkhorn_random(x, y, d = c(2, 4)),
    "product_mean"
  )
  sinkhorn[i] <- max(m_sinkhorn)
  setTxtProgressBar(pb, i)
}
close(pb)
quantile(sinkhorn, 0.95, type = 1)
