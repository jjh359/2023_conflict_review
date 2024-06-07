library("ggplot2")
library("ggridges")
library("tidyr")
library("dplyr")
library("viridis")

setwd("D:/Documents/Research/conflict_review/rate_estimation/")

data <- read.csv(file = "rates_100trees_udmodel_v2.csv", header = TRUE)

data_long <- data %>% gather(run, rate, run1:run100) %>%
  arrange(rate) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

#Colour by rate
ggplot(data_long, aes(x = rate, y = Family, group = Family, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 1) +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Transition rates inferred across 100 credible trees (Unidirectional model)")

#Colour by family
ggplot(data_long, aes(x = rate, y = Family, group = Family, fill = Family)) +
  geom_density_ridges(scale = 1) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Transition rates inferred across 100 credible trees (Unidirectional model)")

ggsave("Fig2.png", device = "png", width = 10, height = 7, units = "in")
