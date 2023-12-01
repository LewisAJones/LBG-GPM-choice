# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 05_Diversity_figures.R
# Aim: Plot stage-level diversity in our palaeolatitudinal bins (Fig. 2 and S1)
#      and diversity heatmaps (Fig. 3 and S2).
# Last updated: 2023-12-01
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(palaeoverse)
library(deeptime)

# Setup common things for figures ---------------------------------------
source("./R/functions/theme_will.R")
GTS2020_periods <- time_bins(rank = "period") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
GTS2020_eras <- time_bins(rank = "era") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Source data -----------------------------------------------------------
source("./R/functions/data_for_div.R")

# Plot data -------------------------------------------------------------
# sqs plot
p <- ggplot(data = div_sqs_join, aes(x = mid,
                                     y = qD_norm1,
                                     colour = model)) +
  geom_point(aes(shape = Method), size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8,
                         limits = na.exclude(unique(div_sqs_join$model))) +
  scale_shape_discrete(NULL, limits = c("Extrapolation", "Rarefaction", "Observed")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = 0.1)) +
  facet_wrap(~factor(interval_name, levels = rev(stages$interval_name)), ncol = 9) +
  labs(y = "Normalised estimated genus richness",
       x = "Palaeolatitudinal bin") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top")

#Update strip colours
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
for (i in strip_t) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <-
    stages$colour[match(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label,
                        stages$interval_name)]
}
#save
ggsave("./figures/LBGs_sqs.png", g, width = 14.4, height = 17.6)

# raw plot
p <- ggplot(data = div_raw_join, aes(x = mid,
                                     y = n_genera_norm1,
                                     colour = model)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8,
                         limits = na.exclude(unique(div_sqs_join$model))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = 0.1)) +
  facet_wrap(~factor(interval_name, levels = rev(stages$interval_name)), ncol = 9) +
  labs(y = "Normalised raw genus richness",
       x = "Palaeolatitudinal bin") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top")

#Update strip colours
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
for (i in strip_t) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <-
    stages$colour[match(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label,
                        stages$interval_name)]
}
#save
ggsave("./figures/LBGs_raw.png", g, width = 14.4, height = 17.6)
