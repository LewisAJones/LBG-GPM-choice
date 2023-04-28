# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 04_diversity_figures.R
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(grid)
library(ggplot2)
library(dplyr)

# Data ------------------------------------------------------------------
# Load bins
lats <- readRDS("./data/lat_bins.RDS")
stages <- readRDS("./data/time_bins.RDS")

# Load data
div_sqs <- readRDS("./data/processed/interpolations.RDS")
div_raw <- readRDS("./data/processed/genus_counts.RDS")

# Wrangling -------------------------------------------------------------
# Set factor levels
stages$interval_name <- factor(x = stages$interval_name, levels = stages$interval_name)

# Join datasets
div_raw_global <- div_raw %>%
  filter(is.na(paleolat_bin)) %>%
  select(-paleolat_bin)
div_raw_join <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  inner_join(div_raw_global, by = c("model", "stage_bin"), suffix = c("", ".global")) %>%
  group_by(model, stage_bin) %>%
  mutate(n_genera_norm1 = n_genera / max(n_genera, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(n_genera_norm2 = n_genera / n_genera.global) %>%
  full_join(stages, by = c("stage_bin" = "bin")) %>%
  full_join(lats, by = c("paleolat_bin" = "bin"))

div_sqs_global <- div_sqs %>%
  filter(is.na(paleolat_bin)) %>%
  select(-paleolat_bin)
div_sqs_join <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  inner_join(div_sqs_global, by = c("model", "stage"), suffix = c("", ".global")) %>%
  group_by(model, stage) %>%
  mutate(qD_norm1 = qD / max(qD, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(qD_norm2 = qD / qD.global) %>%
  full_join(stages, by = c("stage" = "bin")) %>%
  full_join(lats, by = c("paleolat_bin" = "bin"))

# Plot data -------------------------------------------------------------
# sqs plot
p <- ggplot(data = div_sqs_join, aes(x = mid,
                                     y = qD_norm1,
                                     colour = model)) +
  geom_line(linewidth = 0.75, alpha = 1) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8,
                         limits = na.exclude(unique(div_sqs_join$model))) +
  facet_wrap(~factor(interval_name, levels = rev(stages$interval_name)), nrow = 10) +
  labs(y = "Proportional estimated genus richness",
       x = "Palaeolatitude (\u00B0)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "top")

#Update strip colours
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
for (i in strip_t) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <-
    stages$colour[match(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label,
                        stages$interval_name)]
}

ggsave("./figures/LBGs_sqs.pdf", g, width = 15, height = 15)

# raw plot
p <- ggplot(data = div_raw_join, aes(x = mid,
                                     y = n_genera_norm1,
                                     colour = model)) +
  geom_line(linewidth = 0.75, alpha = 1) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8,
                         limits = na.exclude(unique(div_sqs_join$model))) +
  facet_wrap(~factor(interval_name, levels = rev(stages$interval_name)), nrow = 10) +
  labs(y = "Proportional raw genus richness",
       x = "Palaeolatitude (\u00B0)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "top")

#Update strip colours
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
for (i in strip_t) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <-
    stages$colour[match(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label,
                        stages$interval_name)]
}

ggsave("./figures/LBGs_raw.pdf", g, width = 15, height = 15)
