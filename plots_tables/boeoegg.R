# packages
library(tidyverse)
library(ggpubr)
library(ggthemes)

# ggplot2 settings
theme_set(theme_bw(base_size = 10))
colpal <- colorblind_pal()(8)[-5][c(2, 3, 4, 6, 1, 5)]

# functions
source("functions/r_wrappers.R")
Rcpp::sourceCpp("functions/sequential_tests.cpp")

# data
boeoegg <- read.csv("data/boeoegg_burn_duration.csv", header = TRUE, sep = ",")
weather <- read.csv("data/zuri_historical_weather.csv", header = TRUE, sep = ",")

## add year 2022, data from https://www.meteoswiss.admin.ch/services-and-publications/applications/ext/climate-tables-homogenized.html
boeoegg <- rbind(
  boeoegg,
  data.frame(
    year = c(2022, 2023),
    burn_duration_seconds = c(2279, 57 * 60)
  )
)

weather <- rbind(
  weather,
  data.frame(
    Year = rep(c(2022, 2023), each = 3),
    Month = rep(6:8, 2),
    Temperature = c(19.5, 21.3, 20.4, rep(NA, 3)),
    Precipitation = c(126.6, 44.8, 75.1, rep(NA, 3))
  )
)

## preprocess
weather <- weather %>%
  filter(Month %in% 6:8 & Year %in% boeoegg$year) %>%
  group_by(Year) %>%
  summarise(temp = mean(Temperature), precip = mean(Precipitation)) %>%
  rename(year = Year)

data <- (merge(boeoegg, weather, by = "year"))

# time series plot
data_plot <- data %>%
  rename(
    # `Temperature (°C)` = temp,
    # `Precipitation (mm)` = precip,
    # `Burning time (s)` = "burn_duration_seconds",
    `Temperature` = temp,
    `Precipitation` = precip,
    `Burning time` = "burn_duration_seconds",
    Year = year
  ) %>%
  # mutate(most_recent = is.na(`Temperature (°C)`)) %>%
  mutate(most_recent = is.na(`Temperature`)) %>%
  gather(key = "var", value = "val", -Year, -most_recent) %>%
  ggplot() +
  geom_line(aes(x = Year, y = val)) +
  geom_text(
    data = tibble(
      most_recent = rep(TRUE, 3),
      var = c("Burning time", "Precipitation", "Temperature"),
      # var = c("Burning time (s)", "Precipitation (mm)", "Temperature (°C)"),
      # pch = c(NA, "?", "?"),
      pch = c(NA, NA, NA),
      Year = rep(2023, 3),
      y = c(2000, 125, 19)
    ),
    aes(x = Year, y = y, label = pch, color = most_recent),
    size = 5
  ) +
  geom_line(
    data = tibble(
      Year = tail(data$year, 2),
      var = c("Burning time", "Burning time"),
      # var = c("Burning time (s)", "Burning time (s)"),
      val = tail(data$burn_duration_seconds, 2),
      most_recent = TRUE
    ),
    aes(x = Year, y = val, color = most_recent),
    lwd = 0.85
  ) +
  geom_point(aes(x = Year, y = val, color = most_recent)) +
  scale_color_manual(values = c("black", colpal[2])) +
  theme(legend.position = "none") +
  # theme(
  #   legend.position = "none",
  #   strip.placement = "outside",
  #   strip.background = element_blank(),
  #   strip.text.y = element_text(size = 6)
  # ) +
  facet_grid(rows = vars(var), scales = "free_y") +
  # facet_wrap(
  #   ~var,
  #   strip.position = "left",
  #   nrow = 3, scales = "free_y"
  # ) +
  labs(
    x = "Year",
    y = element_blank()
  )



# construct test martingales
set.seed(20230417) # date of most recent Sechselaeuten

data <- na.omit(data)
n <- nrow(data)

d <- c(2, 4)
mtemp <- srt_sinkhorn_random(data$burn_duration_seconds, data$temp, d = d)
mprecip <- srt_sinkhorn_random(data$burn_duration_seconds, data$precip, d = d)
mtemp <- get_martingale(mtemp, combination = "product_mean")
mprecip <- get_martingale(mprecip, combination = "product_mean")

mtemp[n]
mprecip[n]
(mtemp[n] + mprecip[n]) / 2

martingales <- data.frame(
  year = data$year,
  Temperature = mtemp,
  Precipitation = mprecip
)
  
martingales_plot <- martingales %>%
  gather(key = "var", value = "martingale", -year) %>%
  ggplot() +
  geom_point(
    aes(x = year, y = martingale, color = var, group = var, shape = var)
  ) +
  geom_line(
    aes(x = year, y = martingale, color = var, group = var, lty = var)
  ) +
  scale_color_manual(values = colpal[1:2]) +
  scale_linetype_manual(values = c(1, 5)) +
  # theme(legend.position = c(0.5, 0.85)) +
  theme(legend.position = "none") +
  labs(
    x = "Year",
    y = "Test Martingales",
    color = element_blank(),
    lty = element_blank(),
    shape = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 2.5))

pdf(file = "sechselaeuten.pdf", width = 8, height = 3)
ggarrange(data_plot, martingales_plot, ncol = 2)
dev.off()

# derandomized p-values
derandomized_test <- function(
    X,
    Y,
    d = c(2, 4, 8, 16),
    B = 10000,
    method = "am",
    combo = "mean_product"
  ) {
  n <- length(X)
  ms <- matrix(nrow = n, ncol = B)
  if (method == "am") {
    for (j in seq_len(B)) {
      ms[, j] <- 1 / 
        cummax(get_martingale(srt_sinkhorn_random(X, Y, d = d), combination = combo))
    }
    if (B > 1) {
      p <- 2 * rowMeans(ms)
    } else {
      p <- rowMeans(ms)
    }
  } else if (method == "gm") {
    for (j in seq_len(B)) {
      ms[, j] <- 
        log(1 / cummax(get_martingale(srt_sinkhorn_random(X, Y, d = d), combination = combo)))
    }
    p <- exp(1) * exp(rowMeans(ms))
  }
  p
}

derandomized_test(
  X = data$temp,
  Y = data$burn_duration_seconds,
  d = d,
  combo = "product_mean"
)
derandomized_test(
  X = data$precip,
  Y = data$burn_duration_seconds,
  d = d,
  combo = "product_mean"
)
derandomized_test(
  X = data$temp,
  Y = data$burn_duration_seconds,
  d = d,
  combo = "product_mean",
  method = "gm"
)
derandomized_test(
  X = data$precip,
  Y = data$burn_duration_seconds,
  d = d,
  combo = "product_mean",
  method = "gm"
)