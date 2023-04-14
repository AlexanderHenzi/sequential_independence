# packages
library(tidyverse)
library(ggpubr)
library(ggthemes)

# ggplot2 settings
theme_set(theme_bw(base_size = 10))
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# load data, functions
load("simulation_results.rda")

# parameters
alpha0 <- 0.05
n <- 500
set.seed(1994)

# pre-format data
data <- results %>%
  gather(key = "combination", value = "time", pm, mp) %>%
  mutate(time = time + 1e9 * (time == 0))
plot_data <- data %>%
  filter(
    sim != "null" &
    alpha == alpha0 &
    ((method %in% c("bet", "ks") & combination == "pm") | 
       (method %in% c("simple", "sinkhorn") & combination != "bo"))
  ) %>%
  mutate(
    method = paste0(method, "_",  combination),
    method = factor(
      method,
      levels = c(
        "ks_pm",
        "bet_pm",
        "simple_mp",
        "simple_pm",
        "sinkhorn_mp",
        "sinkhorn_pm"
      ),
      labels = c(
        "SR",
        "Sequential~BET",
        "Simple~(eta == 1)",
        "Simple~(eta == 0)",
        "Sinkhorn~(eta == 1)",
        "Sinkhorn~(eta == 0)"
      ),
      ordered = TRUE
    ),
    sim = str_replace(sim, "^\\w{1}", toupper)
  ) %>%
  select(-combination)

## Figure in main part of paper
xlims <- seq(150, 2150, 200)
plotlist <- vector("list", 10)
for (l0 in seq_len(10)) {
  plotlist[[l0]] <- ggplot() +
    stat_ecdf(
      data = filter(
        plot_data,
        l == l0 &
        sim %in% c("Circular", "Linear", "Local") &
        method != "Sequential~BET"
      ),
      mapping = aes(x = time, color = method, group = method, lty = method),
      lwd = 0.8
    ) +
    scale_color_manual(values = colpal, label = scales::parse_format()) +
    scale_linetype_manual(values = 1:5, label = scales::parse_format()) +
    coord_cartesian(xlim = c(0, xlims[l0])) +
    facet_grid(cols = vars(sim), rows = vars(l)) +
    theme(legend.position = "bottom", axis.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 1))
  if (l0 != 9) {
    plotlist[[l0]] <- plotlist[[l0]] +
      labs(
        x = element_blank(),
        y = "Rejection rate",
        color = element_blank(),
        linetype = element_blank()
      )
  } else {
    plotlist[[l0]] <- plotlist[[l0]] +
      labs(
        x = "Sample size",
        y = "Rejection rate",
        color = element_blank(),
        linetype = element_blank()
      )
  }
}

rejection_rates <- ggarrange(
  plotlist = plotlist[c(1,5,9)],
  nrow = 3,
  heights = c(1, 1, 1.1),
  common.legend = TRUE,
  legend = "bottom"
)

pdf(file = "rejection_rates.pdf", width = 8, height = 6)
print(rejection_rates)
dev.off()

## Figures in the supplement
xlims <- seq(150, 2150, 200)
plotlist <- vector("list", 10)
for (l0 in seq_len(10)) {
  plotlist[[l0]] <- ggplot() +
    stat_ecdf(
      data = filter(
        plot_data,
        l == l0
      ),
      mapping = aes(x = time, color = method, group = method, lty = method),
      lwd = 0.8
    ) +
    scale_color_manual(
      values = colpal[c(1, 6, 2:5)],
      label = scales::parse_format()
    ) +
    scale_linetype_manual(values = c(1, 6, 2:5), label = scales::parse_format()) +
    coord_cartesian(xlim = c(0, xlims[l0])) +
    facet_grid(cols = vars(sim), rows = vars(l)) +
    theme(legend.position = "bottom", axis.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 1))
  if (l0 != 9) {
    plotlist[[l0]] <- plotlist[[l0]] +
      labs(
        x = element_blank(),
        y = "Rejection rate",
        color = element_blank(),
        linetype = element_blank()
      )
  } else {
    plotlist[[l0]] <- plotlist[[l0]] +
      labs(
        x = "Sample size",
        y = "Rejection rate",
        color = element_blank(),
        linetype = element_blank()
      )
  }
}

rejection_rates <- ggarrange(
  plotlist = plotlist[c(1, 3, 5, 7, 9)],
  nrow = 5,
  heights = c(1, 1, 1, 1, 1.1),
  common.legend = TRUE,
  legend = "bottom"
)

pdf(file = "rejection_rates_supplement.pdf", width = 8, height = 8)
print(rejection_rates)
dev.off()
