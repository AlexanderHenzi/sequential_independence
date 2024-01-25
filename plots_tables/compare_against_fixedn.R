# packages
library(tidyverse)
library(ggpubr)
library(ggthemes)

# ggplot2 settings
theme_set(theme_bw(base_size = 10))
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# load data, functions
source("functions/simulation_examples.R")

load("data/simulations_comparison_bet.rda")
srt <- results %>%
  mutate(method = "sinkhorn") %>%
  rename(reject = rej)

load("data/simulations_bet_results.rda")
bet <- results %>%
  select(-id) %>%
  mutate(method = "bet")

load("data/simulations_xicor_results.rda")
xicor <- results %>%
  select(-id) %>%
  mutate(method = "xicor")

load("data/simulations_hoeffdingsd_results.rda")
hoeffdingsd <- results %>%
  select(-id) %>%
  mutate(method = "hoeddfingsd")

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

# table comparing bet, xicor, srt
ll <- c(1, 3, 5, 7, 9)

bet <- bet %>%
  filter(l %in% ll) %>%
  group_by(sim, l, n) %>%
  summarise(power = mean(reject)) %>%
  spread(key = "n", value = "power") %>%
  ungroup()

xicor <- xicor %>%
  filter(l %in% ll) %>%
  group_by(sim, l, n) %>%
  summarise(power = mean(reject)) %>%
  spread(key = "n", value = "power") %>%
  ungroup()

hoeffdingsd <- hoeffdingsd %>%
  filter(l %in% ll) %>%
  group_by(sim, l, n) %>%
  summarise(power = mean(reject)) %>%
  spread(key = "n", value = "power") %>%
  ungroup()

srt_256 <- srt %>%
  filter(n == 256 & l %in% ll) %>%
  mutate(rejected = reject > 0) %>%
  mutate(reject = n * (reject == 0) + reject * (reject > 0)) %>%
  group_by(n, sim, l) %>%
  summarise(power = mean(rejected), mean = mean(reject)) %>%
  ungroup() %>%
  select(-n)

srt <- srt %>%
  filter(n == 512 & l %in% ll) %>%
  mutate(rejected = reject > 0) %>%
  mutate(reject = n * (reject == 0) + reject * (reject > 0)) %>%
  group_by(n, sim, l) %>%
  summarise(power = mean(rejected), mean = mean(reject)) %>%
  ungroup() %>%
  select(-n)

# table in main text of paper

## estimate KL on simulated data
set.seed(123)
n <- 1e6
ls <- seq(1, 9, 2)
sims <- vector("list", length(ls))
for (i in seq_along(ls)) {
  sims_tmp <- sim_all(n, ls[i])
  sims_tmp <- sims_tmp[names(sims_tmp) != "null"]
  sims_tmp <- mapply(
    function(data, name) {
      cbind(sim = name, data)
    },
    data = sims_tmp,
    name = names(sims_tmp),
    SIMPLIFY = FALSE
  )
  sims_tmp <- do.call(rbind, sims_tmp)
  sims_tmp$l <- ls[i]
  sims[[i]] <- sims_tmp
}
sims <- do.call(rbind, sims)

ds <- c(2, 8)
kl_estimates <- vector("list", length(ds))
for (i in seq_along(ds)) {
  grd <- seq(0, 1, length.out = ds[i] + 1)
  kl_estimates[[i]] <- sims %>%
    group_by(sim, l) %>%
    mutate(
      x = ecdf(x)(x),
      y = ecdf(y)(y),
      x = cut(x, grd, include.lowest = TRUE),
      y = cut(y, grd, include.lowest = TRUE)
    ) %>%
    group_by(sim, l, x, y) %>%
    summarise(q = length(x) / n) %>%
    group_by(sim, l) %>%
    summarise(kl = sum(q * log(q * ds[i]^2), na.rm = TRUE)) %>%
    mutate(d = ds[i]) %>%
    arrange(sim, l)
}
kl_estimates <- do.call(rbind, kl_estimates)
kl <- kl_estimates %>%
  mutate(kl = sprintf("%.2f", kl)) %>%
  spread(key = "d", value = "kl")

## table
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
  merge(
    kl,
    by = c("sim", "l"),
    all.x = TRUE,
    all.y = FALSE
  ) %>%
  select(sim, l, `2`, `8`, `64`, `128`, `256`, `512`, `power`, `mean`) %>%
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
  merge(
    kl,
    by = c("sim", "l"),
    all.x = TRUE,
    all.y = FALSE
  ) %>%
  select(sim, l, `2`, `8`, `64`, `128`, `256`, `512`, `power`, `mean`) %>%
  mutate(sim = str_to_title(sim)) %>%
  arrange(sim, l) %>%
  mutate(
    mean = round(mean),
    sim = replace(sim, setdiff(seq_along(sim), seq(3, length(sim), 5)), ""),
  ) %>%
  write.table(sep = " & ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)

# tables for xicor, hoeffdingsd, sample size 256 in supplement
merge(xicor, srt, by = c("sim", "l")) %>%
  filter(sim != "null") %>%
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

merge(hoeffdingsd, srt, by = c("sim", "l")) %>%
  filter(sim != "null") %>%
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

merge(bet, srt_256, by = c("sim", "l")) %>%
  filter(sim != "null") %>%
  gather(key = "method", value = "pwr", `64`, `128`, `256`, power) %>%
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
  select(sim, l, `64`, `128`, `256`, `power`, `mean`) %>%
  mutate(sim = str_to_title(sim)) %>%
  arrange(sim, l) %>%
  mutate(
    mean = round(mean),
    sim = replace(sim, setdiff(seq_along(sim), seq(3, length(sim), 5)), ""),
  ) %>%
  write.table(sep = " & ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)
