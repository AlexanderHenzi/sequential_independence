# load packages
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(Rcpp)

# load functions
source("functions/r_wrappers.R")
source("functions/simulation_examples.R")
sourceCpp("functions/sequential_tests.cpp")

# global parameters
theme_set(theme_bw(base_size = 12))
l <- 5
d <- c(2, 4, 8, 16)
K <- length(d)
n <- 10000
n0 <- 1
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# data simulation
set.seed(1991)
data <- sim_linear(n = n, l = l)
x <- data$x
y <- data$y

# test martingales
simple <- srt_simple(X = x, Y = y, d = d, n0 = n0)
sinkhorn <- srt_sinkhorn(X = x, Y = y, d = d, n0 = n0)

# data for plot
df_d <- bind_cols(
    as_tibble(set_names(simple[1:4], paste0("simple_", d))),
    as_tibble(set_names(sinkhorn[1:4], paste0("sinkhorn_", d))),
    tibble(n = seq_len(n))
  ) %>%
  pivot_longer(
    cols = starts_with("s"),
    names_to = c("method", "d"),
    names_sep = "_",
    values_to = "martingale"
  ) %>%
  group_by(method, d) %>%
  arrange(n) %>%
  mutate(martingale = cumprod(martingale)) %>%
  ungroup() %>%
  mutate(parameter = paste0("d == ", d))

df_combined <- tibble(
    simple_0 = get_martingale(simple, combination = "product_mean"),
    simple_1 = get_martingale(simple, combination = "mean_product"),
    sinkhorn_0 = get_martingale(sinkhorn, combination = "product_mean"),
    sinkhorn_1 = get_martingale(sinkhorn, combination = "mean_product"),
    n = seq_len(n)
  ) %>%
  pivot_longer(
    cols = starts_with("s"),
    names_to = c("method", "eta"),
    names_sep = "_",
    values_to = "martingale"
  ) %>%
  mutate(parameter = paste0("eta == ", eta))

df <- bind_rows(select(df_d, -d), select(df_combined, -eta)) %>%
  mutate(
    is_combined = grepl("eta", parameter),
    parameter = factor(
      parameter,
      levels = c(paste0("d == ", 2^(1:4)), c("eta == 0", "eta == 1")),
      labels = c(paste0("d == ", 2^(1:4)), c("eta == 0", "eta == 1")),
      ordered = TRUE
    )
  )

# plot
small_n <- 501
small_simple <- ggplot() +
  geom_hline(yintercept = 1, col = "darkgray", lty = 5) +
  geom_line(
    data = filter(df, n < small_n & method == "simple"),
    aes(x = n, y = martingale, color = parameter, lty = parameter)
  ) +
  scale_linetype_manual(values = c(5, 5, 5, 5, 1, 1), labels = scales::parse_format()) +
  scale_y_log10() +
  scale_color_manual(values = colpal, labels = scales::parse_format()) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  labs(x = "Sample size", y = "Test martingales")
large_simple <- ggplot() +
  geom_hline(yintercept = 1, col = "darkgray", lty = 5) +
  geom_line(
    data = filter(df, method == "simple"),
    aes(x = n, y = martingale, color = parameter, lty = parameter)
  ) +
  scale_linetype_manual(values = c(5, 5, 5, 5, 1, 1), labels = scales::parse_format()) +
  scale_y_log10() +
  scale_color_manual(values = colpal, labels = scales::parse_format()) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  labs(x = "Sample size", y = "Test martingales")
small_large_n <- ggarrange(
  small_simple,
  large_simple,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

pdf(file = "illustration_small_large_n.pdf", width = 8, height = 3.5)
print(small_large_n)
dev.off()

simple_sinkhorn <- ggplot() +
  geom_hline(yintercept = 1, col = "darkgray", lty = 5) +
  geom_line(
    data = mutate(
      filter(df, n < 1000 & !grepl("eta", parameter)),
      method = str_replace(method, "s", "S")
    ),
    aes(x = n, y = martingale, color = parameter, lty = method, group = interaction(parameter, method))
  ) +
  scale_linetype_manual(values = c(5, 1), labels = scales::parse_format()) +
  scale_shape_manual(values = c(16, 18)) +
  scale_y_log10() +
  scale_color_manual(values = colpal, labels = scales::parse_format()) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1), lty = guide_legend(nrow = 1)) +
  labs(x = "Sample size", y = "Test martingales")

pdf("illustration_simple_sinkhorn.pdf", width = 8, height = 3.5)
print(simple_sinkhorn)
dev.off()

# density plots
nn <- 500
simple <- srt_simple(X = x[seq_len(nn)], Y = y[seq_len(nn)], d = d, n0 = n0)
sinkhorn <- srt_sinkhorn(X = x[seq_len(nn)], Y = y[seq_len(nn)], d = d, n0 = n0)

dat_simple <- dat_sinkhorn <- vector("list", K + 1)
dens_simple <- dens_sinkhorn <- vector("list", K + 1)
for (j in seq_len(K)) {
  grd <- expand.grid(r = seq(2^(-j), 1, 2^(-j)), s = seq(2^(-j), 1, 2^(-j)))
  dat <- tibble(grd, density = 2^{2*j} * c(simple[[length(d) + j]]), d = 2^j)
  dens_simple[[j]] <- ggplot() +
    geom_raster(
      data = dat,
      aes(x = r, y = s, fill = density),
      hjust = 0,
      vjust = 0
    ) +
    scale_fill_gradientn(colors = c("white", "black"), limits = c(0, 3)) +
    theme(legend.position = "bottom") +
    labs(x = element_blank(), y = element_blank(), fill = "Density")
  dat <- tibble(grd, density = 2^{2*j} * c(sinkhorn[[length(d) + j]]), d = 2^j)
  dens_sinkhorn[[j]] <- ggplot() +
    geom_raster(
      data = dat,
      aes(x = r, y = s, fill = density),
      hjust = 0,
      vjust = 0
    ) +
    scale_fill_gradientn(colors = c("white", "black"), limits = c(0, 3)) +
    theme(legend.position = "bottom") +
    labs(x = element_blank(), y = element_blank(), fill = "Density")
  if (j > 1) {
    dens_simple[[j]] <- dens_simple[[j]] + theme(axis.text.y = element_blank())
    dens_sinkhorn[[j]] <- dens_sinkhorn[[j]] + theme(axis.text.y = element_blank())
  }
}
dens_simple_combined <- simple[[2 * K]] + 2^(-2*K)
dens_sinkhorn_combined <- sinkhorn[[2 * K]] + 2^(-2*K)

for (j in seq_len(K - 1)) {
  b <- 2^(K - j)
  for (k in seq_len(2^j)) {
    for (l in seq_len(2^j)) {
      ks <- (1 + b * (k - 1)):(b * k)
      ls <- (1 + b * (l - 1)):(b * l)
      dens_simple_combined[ks, ls] <- dens_simple_combined[ks, ls] +
        simple[[K + j]][k, l] / b^2
      dens_sinkhorn_combined[ks, ls] <- dens_sinkhorn_combined[ks, ls] +
        sinkhorn[[K + j]][k, l] / b^2
    }
  }
}
dens_simple_combined <- dens_simple_combined / (K + 1)
dens_sinkhorn_combined <- dens_sinkhorn_combined / (K + 1)

dat <- tibble(grd, density = 2^(2*K) * c(dens_simple_combined))
dens_simple[[K+1]] <- ggplot() +
  geom_raster(
    data = dat,
    aes(x = r, y = s, fill = density),
    hjust = 0,
    vjust = 0
  ) +
  scale_fill_gradientn(colors = c("white", "black"), limits = c(0, 3)) +
  theme(legend.position = "bottom", axis.text.y = element_blank()) +
  labs(x = element_blank(), y = element_blank(), fill = "Density")
dat <- tibble(grd, density = 2^(2*K) * c(dens_sinkhorn_combined))
dens_sinkhorn[[K+1]] <- ggplot() +
  geom_raster(
    data = dat,
    aes(x = r, y = s, fill = density),
    hjust = 0,
    vjust = 0
  ) +
  scale_fill_gradientn(colors = c("white", "black"), limits = c(0, 3)) +
  theme(legend.position = "bottom", axis.text.y = element_blank()) +
  labs(x = element_blank(), y = element_blank(), fill = "Density")

densities <- ggarrange(
  plotlist = dens_simple,
  ncol = 5,
  widths = c(1.15, 1, 1, 1, 1),
  common.legend = TRUE,
  legend = "bottom"
)

pdf(file = "illustration_densities.pdf", width = 11, height = 3)
print(densities)
dev.off()

# illustration of binning and counting
set.seed(123)
k <- 6
R <- pmin(0.9, pmax(0.1, runif(k)))
S <- pmin(0.9, pmax(0.1, runif(k)))
points_data <- tibble(
  ds = rep(c("d = 2", "d = 4"), each = k),
  R = rep(R, 2),
  S = rep(S, 2),
  time_order = rep(factor(seq_len(k)), 2),
  labels2 = rep(paste0("paste('(',", paste0("R[", seq_len(k), "])")), 2),
  labels3 = rep("       , ", 2 * k),
  labels4 = rep(paste0("~~~~~~~~~~~~~paste(S[", seq_len(k), "], ')')"), 2),
)

count_grid_cells <- ggplot(points_data) +
  geom_text(
    aes(x = R, y = S, color = time_order, label = labels2),
    cex = 3,
    parse = TRUE
  ) +
  geom_text(aes(x = R, y = S, color = time_order, label = labels3), cex = 3) +
  geom_text(
    aes(x = R, y = S, color = time_order, label = labels4),
    cex = 3,
    parse = TRUE
  ) +
  facet_grid(cols = vars(ds)) +
  geom_segment(
    data = tibble(
      ds = c("d = 2", rep("d = 4", 3)),
      yend = c(0.5, 0.25, 0.5, 0.75),
      y = yend,
      x = rep(0, 4),
      xend = rep(1, 4)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "darkgray"
  ) +
  geom_segment(
    data = tibble(
      ds = c("d = 2", rep("d = 4", 3)),
      xend = c(0.5, 0.25, 0.5, 0.75),
      x = xend,
      y = rep(0, 4),
      yend = rep(1, 4)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "darkgray"
  ) +
  geom_point(aes(x = R + 0.025, y = S + 0.0375, color = time_order)) +
  geom_rect(
    data = tibble(
      xmin = c(0, 0),
      xmax = c(1, 1),
      ymin = c(0, 0),
      ymax = c(1, 1)
    ),
    aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
    fill = "transparent",
    color = "black"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = colpal) +
  labs(x = element_blank(), y = element_blank()) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

pdf("count_grid_cells.pdf", width = 8, height = 4.25)
print(count_grid_cells)
dev.off()

# different discretizations and kl divergence
set.seed(1991)
n <- 1e4
ls <- c(1, 3, 5)
ds <- c(2, 4, 8, 16)
sims <- vector("list", length(ls))
for (i in seq_along(ls)) {
  sims_tmp <- sim_linear(n, ls[i])
  sims_tmp$l <- ls[i]
  sims[[i]] <- sims_tmp
}
sims <- do.call(rbind, sims)
kl_estimates <- vector("list", length(ds))
for (i in seq_along(ds)) {
  grd <- seq(0, 1, length.out = ds[i] + 1)
  kl_estimates[[i]] <- sims %>%
    group_by(l) %>%
    mutate(
      x = ecdf(x)(x),
      y = ecdf(y)(y),
      x = cut(x, grd, include.lowest = TRUE),
      y = cut(y, grd, include.lowest = TRUE)
    ) %>%
    group_by(l, x, y) %>%
    summarise(q = length(x) / n) %>%
    group_by(l) %>%
    summarise(kl = sum(q * log(q * ds[i]^2), na.rm = TRUE)) %>%
    mutate(d = ds[i]) %>%
    arrange(l)
}
kl_estimates <- do.call(rbind, kl_estimates)

n <- 10000
mndl <- vector("list", length(ls))
for (i in seq_along(ls)) {
  data <- sim_linear(n = n, l = ls[i])
  x <- data$x
  y <- data$y
  simple <- srt_simple(X = x, Y = y, d = d, n0 = n0)
  fs <- unlist(simple[seq_along(ds)])
  mndl[[i]] <- data.frame(
    fnd = fs,
    d = rep(ds, each = n),
    l = ls[i],
    index = rep(seq_len(n), length(ds))
  )
}
mndl <- do.call(rbind, mndl)
mndl <- mndl %>%
  group_by(d, l) %>%
  arrange(index) %>%
  mutate(lfnd = cummean(log(fnd)))

kl_growth_rate <- ggplot() +
  geom_line(
    data = mutate(mndl, d = factor(d), l = factor(l)),
    aes(x = index, y = lfnd, color = d, group = d, linetype = d)
  ) +
  geom_hline(
    data = mutate(kl_estimates, d = factor(d), l = factor(l)),
    aes(yintercept = kl, color = d, group = d, linetype = d),
    alpha = 0.75
  ) +
  scale_color_manual(values = colpal, labels = scales::parse_format()) +
  scale_linetype_manual(values = c(5, 1, 4, 3)) +
  facet_wrap(.~l, scales = "free_y") +
  labs(
    x = expression("Index"~italic(n)),
    y = "KL-divergence, empirical growth rate"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") 

pdf("kl_growth_rate.pdf", width = 8, height = 3)
print(kl_growth_rate)
dev.off()
