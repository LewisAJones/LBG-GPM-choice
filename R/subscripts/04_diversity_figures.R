# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 04_diversity_figures.R
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(palaeoverse)
library(deeptime)

# Data ------------------------------------------------------------------
# Load bins
lats <- readRDS("./data/lat_bins.RDS")
stages <- readRDS("./data/time_bins.RDS")

# Load data
div_sqs <- readRDS("./data/processed/interpolations.RDS")
div_raw <- readRDS("./data/processed/genus_counts.RDS")

# Setup common things for figures ---------------------------------------
source("./R/functions/theme_will.R")
GTS2020_periods <- time_bins(rank = "period") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
GTS2020_eras <- time_bins(rank = "era") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

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
  complete(stage_bin = stages$bin, paleolat_bin, model) %>%
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
  complete(stage = stages$bin, paleolat_bin, model) %>%
  left_join(stages, by = c("stage" = "bin")) %>%
  left_join(lats, by = c("paleolat_bin" = "bin"))

# Plot data -------------------------------------------------------------
# sqs plot
p <- ggplot(data = div_sqs_join, aes(x = mid,
                                     y = qD_norm1,
                                     colour = model)) +
  geom_point(aes(shape = Method), size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8,
                         limits = na.exclude(unique(div_sqs_join$model))) +
  scale_shape_discrete(NULL, limits = c("Extrapolation", "Rarefaction")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = 0.1)) +
  facet_wrap(~factor(interval_name, levels = rev(stages$interval_name)), nrow = 10) +
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

ggsave("./figures/LBGs_sqs.png", g, width = 16, height = 16)

# raw plot
p <- ggplot(data = div_raw_join, aes(x = mid,
                                     y = n_genera_norm1,
                                     colour = model)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8,
                         limits = na.exclude(unique(div_sqs_join$model))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = 0.1)) +
  facet_wrap(~factor(interval_name, levels = rev(stages$interval_name)), nrow = 10) +
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

ggsave("./figures/LBGs_raw.png", g, width = 16, height = 16)

# Heatmap --------------------------------------------------------
# sqs plot
gg_heatmap_sqs <- ggplot(data = div_sqs_join) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = qD_norm1)) +
  geom_hline(yintercept = 3.5) +
  annotate(geom = "text", x = 538, y = 0.25, label = "S. Hemisphere", hjust = 0, size = 5) +
  annotate(geom = "text", x = 538, y = 6.9, label = "N. Hemisphere", hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                   expand = expansion(add = 1.25)) +
  scale_fill_viridis_c("Norm. est. genus richness", limits = c(0, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0))) +
  facet_wrap(~model, ncol = 1)
ggsave("./figures/heatmap_SQS.png", gg_heatmap_sqs, width = 13, height = 13.5)

# raw plot
gg_heatmap_raw <- ggplot(data = div_raw_join %>% mutate(n_genera_norm1 = ifelse(n_genera_norm1 == 0, NA, n_genera_norm1))) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = n_genera_norm1)) +
  geom_hline(yintercept = 3.5) +
  annotate(geom = "text", x = 538, y = 0.25, label = "S. Hemisphere", hjust = 0, size = 5) +
  annotate(geom = "text", x = 538, y = 6.9, label = "N. Hemisphere", hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                   expand = expansion(add = 1.25)) +
  scale_fill_viridis_c("Norm. raw genus richness", limits = c(0, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0))) +
  facet_wrap(~model, ncol = 1)
ggsave("./figures/heatmap_raw.png", gg_heatmap_raw, width = 13, height = 13.5)

# Difference heatmaps ####
diffs_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin), !is.na(qD)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  group_by(model, stage) %>%
  mutate(qD_norm = qD / max(qD, na.rm = TRUE)) %>%
  ungroup() %>%
  inner_join(., ., by = c("stage" = "stage", "paleolat_bin" = "paleolat_bin"),
             relationship = "many-to-many") %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  group_by(stage, paleolat_bin, model.x, model.y) %>%
  summarise(diff = qD_norm.x - qD_norm.y, .groups = "drop") %>%
  complete(stage = stages$bin, model.x, model.y) %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  mutate(models = paste(model.x, model.y, sep = " - ")) %>%
  left_join(stages, by = c("stage" = "bin")) %>%
  left_join(lats, by = c("paleolat_bin" = "bin"))

diffs_heatmap_sqs <- ggplot(data = diffs_sqs) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = diff)) +
  geom_hline(yintercept = 3.5) +
  annotate(geom = "text", x = 538, y = 0.25, label = "S. Hemisphere", hjust = 0, size = 5) +
  annotate(geom = "text", x = 538, y = 6.9, label = "N. Hemisphere", hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                   expand = expansion(add = 1.25)) +
  scale_fill_viridis_c("Norm. est. genus richness", limits = c(-1, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0))) +
  facet_wrap(~models, ncol = 1)
ggsave("./figures/diffs_heatmap_sqs.png", diffs_heatmap_sqs, width = 13, height = 13.5)

diffs_raw <- div_raw %>%
  filter(!is.na(paleolat_bin), !is.na(n_genera), n_genera > 0) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  group_by(model, stage_bin) %>%
  mutate(n_genera_norm = n_genera / max(n_genera, na.rm = TRUE)) %>%
  ungroup() %>%
  inner_join(., ., by = c("stage_bin" = "stage_bin", "paleolat_bin" = "paleolat_bin"),
             relationship = "many-to-many") %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  group_by(stage_bin, paleolat_bin, model.x, model.y) %>%
  summarise(diff = n_genera_norm.x - n_genera_norm.y, .groups = "drop") %>%
  complete(stage_bin = stages$bin, model.x, model.y) %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  mutate(models = paste(model.x, model.y, sep = " - ")) %>%
  left_join(stages, by = c("stage_bin" = "bin")) %>%
  left_join(lats, by = c("paleolat_bin" = "bin"))

diffs_heatmap_raw <- ggplot(data = diffs_raw) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = diff)) +
  geom_hline(yintercept = 3.5) +
  annotate(geom = "text", x = 538, y = 0.25, label = "S. Hemisphere", hjust = 0, size = 5) +
  annotate(geom = "text", x = 538, y = 6.9, label = "N. Hemisphere", hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                   expand = expansion(add = 1.25)) +
  scale_fill_viridis_c("Norm. raw genus richness", limits = c(-1, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0))) +
  facet_wrap(~models, ncol = 1)
ggsave("./figures/diffs_heatmap_raw.png", diffs_heatmap_raw, width = 13, height = 13.5)

