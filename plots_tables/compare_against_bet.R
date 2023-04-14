# packages
library(tidyverse)
library(ggpubr)
library(ggthemes)

# ggplot2 settings
theme_set(theme_bw(base_size = 10))
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# load data, functions
ll <- c(1, 3, 5, 7, 9)

load("data/simulations_comparison_bet.rda")
srt <- results %>%
  mutate(method = "sinkhorn") %>%
  rename(reject = rej)
load("data/simulations_bet_results.rda")
bet <- results %>%
  select(-id) %>%
  mutate(method = "bet")

# check validity of srt under independence
srt %>%
  filter(sim == "null") %>%
  group_by(l, n) %>%
  summarise(rej = mean(reject > 0), se = sqrt(rej * (1 - rej) / 1000)) %>%
  mutate(too_high = rej - se > 0.05) %>%
  print(n = 30)

srt %>%
  filter(sim == "null") %>%
  summarise(rej = mean(reject > 0))

# table comparing bet and srt
bet <- bet %>%
  filter(l %in% ll) %>%
  group_by(sim, l, n) %>%
  summarise(power = mean(reject)) %>%
  spread(key = "n", value = "power") %>%
  ungroup()

srt <- srt %>%
  filter(n == 512 & l %in% ll) %>%
  mutate(rejected = reject > 0) %>%
  mutate(reject = n * (reject == 0) + reject * (reject > 0)) %>%
  group_by(n, sim, l) %>%
  summarise(power = mean(rejected), mean = mean(reject)) %>%
  ungroup() %>%
  select(-n)

# table in main text of paper
merge(bet, srt, by = c("sim", "l")) %>%
  filter(sim %in% c("local", "linear", "circular")) %>%
  gather(key = "method", value = "pwr", `64`, `128`, `256`, `512`, power) %>%
  mutate(
    pwr = as.character(round(pwr, 2)),
    pwr = sapply(
      pwr,
      function(x) {
        if (nchar(x) > 1)
          return(paste0(x, paste0(rep(0, 4 - nchar(x)), collapse = "")))
        paste0(x, ".", paste0(rep(0, 3 - nchar(x)), collapse = ""))
      }
    )
  ) %>%
  spread(key = "method", value = "pwr") %>%
  select(sim, l, `64`, `128`, `256`, `512`, `power`, `mean`) %>%
  mutate(sim = str_to_title(sim)) %>%
  arrange(sim, l) %>%
  mutate(
    mean = round(mean),
    sim = replace(sim, setdiff(seq_along(sim), seq(3, length(sim), 5)), ""),
  ) %>%
  write.table(sep = " & ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)

# table in supplement
merge(bet, srt, by = c("sim", "l")) %>%
  filter(!sim %in% c("local", "linear", "circular", "null")) %>%
  gather(key = "method", value = "pwr", `64`, `128`, `256`, `512`, power) %>%
  mutate(
    pwr = as.character(round(pwr, 2)),
    pwr = sapply(
      pwr,
      function(x) {
        if (nchar(x) > 1)
          return(paste0(x, paste0(rep(0, 4 - nchar(x)), collapse = "")))
        paste0(x, ".", paste0(rep(0, 3 - nchar(x)), collapse = ""))
      }
    )
  ) %>%
  spread(key = "method", value = "pwr") %>%
  select(sim, l, `64`, `128`, `256`, `512`, `power`, `mean`) %>%
  mutate(sim = str_to_title(sim)) %>%
  arrange(sim, l) %>%
  mutate(
    mean = round(mean),
    sim = replace(sim, setdiff(seq_along(sim), seq(3, length(sim), 5)), ""),
  ) %>%
  write.table(sep = " & ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)
