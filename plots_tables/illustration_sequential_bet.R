# packages
library(tidyverse)
library(ggpubr)
library(ggthemes)

# ggplot2 settings
theme_set(theme_bw(base_size = 10))
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

dmax <- 2
zero_one <- vector("list", (2^dmax - 1)^2)
lwr <- 1
H <- matrix(1, nrow = 1, ncol = 1)
h <- matrix(nrow = 2, ncol = 2, c(1, 1, 1, -1))
for (j in seq_len(dmax)) {
  d <- 2^j
  H <- (H %x% h) %x% h
  drop_cols <- unique(c(
    c(seq_len(d), seq(d + 1, d^2, d)),
    seq(1, d - 1, 2) + rep((2 * d) * (0:(d/2 - 1)), each = d/2)
  ))
  basis_mat <- H[, -drop_cols, drop = FALSE]
  zero_one_j <- lapply(
    asplit(basis_mat, 2),
    function(x) matrix(nrow = d, ncol = d, replace(x, x == -1, 0))
  )
  zero_one[lwr:(lwr + length(zero_one_j) - 1)] <- zero_one_j
  lwr <- lwr + length(zero_one_j)
}
zero_one[[1]] <- 
  cbind(c(1, 1, 0, 0), c(1, 1, 0, 0), c(0, 0, 1, 1), c(0, 0, 1, 1))

grd <- expand.grid(x = seq(1/d, 1, 1/d), y = seq(1/d, 1, 1/d))
zero_one <- do.call(rbind, imap(zero_one, ~cbind(grd, col = c(.x), ind = .y)))

illustration_seq_bet <- ggplot() +
  geom_raster(
    data = zero_one,
    mapping = aes(x = x, y = y, fill = factor(col))
  ) +
  facet_wrap(.~ind, nrow = d - 1, ncol = d - 1) +
  scale_fill_manual(values = colpal[1:2]) +
  labs(
    x = element_blank(),
    y = element_blank()
  ) +
  theme(legend.position = "none")

# export figure
pdf(file = "illustration_seq_bet.pdf", width = 8, height = 8)
print(illustration_seq_bet)
dev.off()
