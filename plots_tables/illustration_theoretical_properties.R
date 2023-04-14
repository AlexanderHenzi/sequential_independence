# load packages
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(Rcpp)

# load functions
source("final_code/r_wrappers.R")
source("final_code/simulation_examples.R")
sourceCpp("final_code/sequential_tests.cpp")

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