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
weather_april <- read.csv("data/climate-reports-tables-homogenized_SMA.txt", skip = 25, sep = "")

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

weather_april <- rbind(
  weather,
  data.frame(
    Year = c(2022, 2023),
    Month = c(4, 4),
    Temperature = c(13, NA),
    Precipitation = c(56.8, NA)
  )
)

## preprocess
weather <- weather %>%
  filter(Month %in% 6:8 & Year %in% boeoegg$year) %>%
  group_by(Year) %>%
  summarise(temp = mean(Temperature), precip = mean(Precipitation)) %>%
  rename(year = Year) %>%
  mutate(Period = "Summer")

april_weather <- weather_april %>%
  filter(Month == 4 & Year %in% boeoegg$year) %>%
  group_by(Year) %>%
  summarise(temp = mean(Temperature), precip = mean(Precipitation)) %>%
  rename(year = Year) %>%
  mutate(Period = "April")

weather <- rbind(weather, april_weather)
data <- (merge(boeoegg, weather, by = "year"))
m <- nrow(boeoegg)

# time series plot
ts_plot_data <- data %>%
  filter(year != 2023) %>%
  bind_rows(data.frame(
    year = c(2023, 2023),
    temp = c(mean(c(20.0, 20.3, 19.9)), 13),
    precip = c(mean(c(37.7, 131.5, 162.8)), 56.8),
    burn_duration_seconds = 57 * 60,
    Period = c("Summer", "April")
  )) %>%
  rename(
    `Temperature` = temp,
    `Precipitation` = precip,
    `Burning time` = "burn_duration_seconds",
    Year = year
  ) %>%
  gather(key = "var", value = "val", Temperature, Precipitation, `Burning time`)

data_plot <- ts_plot_data %>%
  ggplot() +
  geom_line(aes(
    x = Year,
    y = val,
    group = Period,
    linetype = Period,
    color = Period)
  ) +
  geom_point(
    aes(x = Year, y = val, group = Period, color = Period, shape = var),
    cex = 1
  ) +
  geom_point(
    aes(x = Year, y = val, group = Period, color = Period, shape = var),
    cex = 1,
    data = filter(ts_plot_data, Period == "Summer")
  ) +
  scale_linetype_manual(values = c(5, 1)) +
  scale_color_manual(values = c(colpal[7], "black")) +
  scale_shape_manual(values = c(15, 16, 17)) +
  theme(legend.position = "none") +
  facet_grid(rows = vars(var), scales = "free_y") +
  labs(x = "Year", y = element_blank())

# construct test martingales
data_all <- na.omit(data)
d <- c(2, 4)

## seed: date of most resent sechselaeuten
set.seed(20230417) 

## Summer: first everything up to 2022
data <- filter(data_all, Period == "Summer")
n <- nrow(data)
mtemp <- srt_sinkhorn_random(data$burn_duration_seconds, data$temp, d = d)
mprecip <- srt_sinkhorn_random(data$burn_duration_seconds, data$precip, d = d)

## include year 2023
temp <- mean(c(20.0, 20.3, 19.9))
precip <- mean(c(37.7, 131.5, 162.8))
burn <- 57 * 60
data$year[n] <- 2023

r_burn <- rank(c(data$burn_duration_seconds[-n], burn))[n] / n - 
  runif(1, 0, 1 / n)
r_temp <- rank(c(data$temp[-n], temp))[n] / n - runif(1, 0, 1 / n)
r_precip <- rank(c(data$precip[-n], precip))[n] / n - runif(1, 0, 1 / n)
p_burn_2 <- findInterval(r_burn, c(0.5, 1)) + 1
p_burn_4 <- findInterval(r_burn, c(0.25, 0.5, 0.75, 1)) + 1
p_temp_2 <- findInterval(r_temp, c(0.5, 1)) + 1
p_temp_4 <- findInterval(r_temp, c(0.25, 0.5, 0.75, 1)) + 1
p_precip_2 <- findInterval(r_precip, c(0.5, 1)) + 1
p_precip_4 <- findInterval(r_precip, c(0.25, 0.5, 0.75, 1)) + 1

f_temp_2 <- 2 * mtemp[[3]][p_burn_2, p_temp_2]
f_temp_4 <- 2 * mtemp[[4]][p_burn_2, p_temp_4]
f_precip_2 <- 4 * mprecip[[3]][p_burn_2, p_precip_2]
f_precip_4 <- 4 * mprecip[[4]][p_burn_2, p_precip_4]

mtemp[[1]][n] <- f_temp_2
mtemp[[2]][n] <- f_temp_4
mtemp[[1]][n] <- f_precip_2
mtemp[[2]][n] <- f_precip_4

mtemp <- get_martingale(mtemp, combination = "product_mean")
mprecip <- get_martingale(mprecip, combination = "product_mean")

round(mtemp[n], 2)
round(mprecip[n], 2)
(round(mtemp[n], 2) + round(mprecip[n], 2)) / 2

martingales <- data.frame(
  year = data$year,
  Temperature = mtemp,
  Precipitation = mprecip,
  Period = "Summer"
)

## April
set.seed(20230417)

data <- filter(data_all, Period == "April")
mtemp <- srt_sinkhorn_random(data$burn_duration_seconds, data$temp, d = d)
mprecip <- srt_sinkhorn_random(data$burn_duration_seconds, data$precip, d = d)

## include year 2023
temp <- 6
precip <- 198.1
burn <- 57 * 60
data$year[n] <- 2023

r_burn <- rank(c(data$burn_duration_seconds[-n], burn))[n] / n - 
  runif(1, 0, 1 / n)
r_temp <- rank(c(data$temp[-n], temp))[n] / n - runif(1, 0, 1 / n)
r_precip <- rank(c(data$precip[-n], precip))[n] / n - runif(1, 0, 1 / n)
p_burn_2 <- findInterval(r_burn, c(0.5, 1)) + 1
p_burn_4 <- findInterval(r_burn, c(0.25, 0.5, 0.75, 1)) + 1
p_temp_2 <- findInterval(r_temp, c(0.5, 1)) + 1
p_temp_4 <- findInterval(r_temp, c(0.25, 0.5, 0.75, 1)) + 1
p_precip_2 <- findInterval(r_precip, c(0.5, 1)) + 1
p_precip_4 <- findInterval(r_precip, c(0.25, 0.5, 0.75, 1)) + 1

f_temp_2 <- 2 * mtemp[[3]][p_burn_2, p_temp_2]
f_temp_4 <- 2 * mtemp[[4]][p_burn_2, p_temp_4]
f_precip_2 <- 4 * mprecip[[3]][p_burn_2, p_precip_2]
f_precip_4 <- 4 * mprecip[[4]][p_burn_2, p_precip_4]

mtemp[[1]][n] <- f_temp_2
mtemp[[2]][n] <- f_temp_4
mtemp[[1]][n] <- f_precip_2
mtemp[[2]][n] <- f_precip_4

mtemp <- get_martingale(mtemp, combination = "product_mean")
mprecip <- get_martingale(mprecip, combination = "product_mean")

round(mprecip[n], 2)

martingales <- rbind(
  martingales,
  data.frame(
    year = data$year,
    Temperature = mtemp,
    Precipitation = mprecip,
    Period = "April"
  )
)
  
martingales_plot <- martingales %>%
  gather(key = "var", value = "martingale", Temperature, Precipitation) %>%
  ggplot() +
  geom_hline(yintercept = 1, col = "darkgray") +
  geom_point(
    aes(
      x = year,
      y = martingale,
      color = Period,
      group = interaction(var, Period),
      shape = var
    ),
    cex = 1
  ) +
  geom_line(aes(
    x = year,
    y = martingale,
    color = Period,
    group = interaction(var, Period),
    lty = Period
  )) +
  scale_color_manual(values = c(colpal[7], "black")) +
  scale_linetype_manual(values = c(5, 1)) +
  scale_shape_manual(values = c(16, 17)) +
  theme(legend.position = "none") +
  labs(
    x = "Year",
    y = "Test Martingales",
    color = element_blank(),
    lty = element_blank(),
    shape = element_blank()
  ) +
  scale_y_log10()

pdf(file = "sechselaeuten.pdf", width = 8, height = 2.75)
ggarrange(data_plot, martingales_plot, ncol = 2)
dev.off()

## test the method by Shekhar and Ramdas (2023) for April precipitation
x <- data$burn_duration_seconds[data$Period == "April"]
y <- data$precip[data$Period == "April"]

xgrid <- seq(min(x), max(x), length.out = 1001)
ygrid <- seq(min(y), max(y), length.out = 1001)

sr_martingale <- get_martingale(seq_ks(x, y, xgrid, ygrid))
sr_martingale
max(sr_martingale)
