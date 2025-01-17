# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_Diversity_and_Difference_heatmaps.R
# Aim: Plot stage-level  heatmaps of diversity (Fig. 3 and S2) and pairwise 
#      difference in normalised diversity between models (Fig. 4 and S3).
# Last updated: 2023-12-01
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(palaeoverse)
library(deeptime)
library(ggh4x)

# Setup common things for figures ---------------------------------------
source("./R/functions/theme_will.R")
ics_periods <- time_bins(scale = "international periods") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
ics_eras <- time_bins(scale = "international eras") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Source data -----------------------------------------------------------
source("./R/functions/data_for_div.R")

# Heatmap ---------------------------------------------------------------
# sqs plot
gg_heatmap_sqs <- ggplot(data = div_sqs_join) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = qD_norm1)) +
  geom_hline(yintercept = 3.5) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0)) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High")) +
  scale_fill_viridis_c("Normalised estimated genus richness", limits = c(0, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  guides(y.sec = guide_axis_manual(breaks = factor(sort(lats$mid))[c(2, 5)], 
                                   labels = c("S. Hemi.", "N. Hemi."),
                                   angle = 90, label_size = 16)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0)),
             axis.ticks.length.y.right = unit(0, "in")) +
  facet_wrap(~model, ncol = 1)
ggsave("./figures/heatmap_SQS.png", gg_heatmap_sqs, width = 13, height = 13.5)

# raw plot
gg_heatmap_raw <- ggplot(data = div_raw_join %>% mutate(n_genera_norm1 = ifelse(n_genera_norm1 == 0, NA, n_genera_norm1))) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = n_genera_norm1)) +
  geom_hline(yintercept = 3.5) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0)) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High")) +
  scale_fill_viridis_c("Normalised raw genus richness", limits = c(0, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  guides(y.sec = guide_axis_manual(breaks = factor(sort(lats$mid))[c(2, 5)], 
                                   labels = c("S. Hemi.", "N. Hemi."),
                                   angle = 90, label_size = 16)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0)),
             axis.ticks.length.y.right = unit(0, "in")) +
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
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  group_by(stage, paleolat_bin, model.x, model.y) %>%
  summarise(diff = qD_norm.x - qD_norm.y, .groups = "drop") %>%
  complete(stage = stages$bin, model.x, model.y) %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  mutate(models = paste(model.x, model.y, sep = " - ")) %>%
  left_join(stages, by = c("stage" = "bin")) %>%
  left_join(lats, by = c("paleolat_bin" = "bin"))

diffs_heatmap_sqs <- ggplot(data = diffs_sqs) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = diff)) +
  geom_hline(yintercept = 3.5) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0)) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High")) +
  scale_fill_viridis_c("Difference in normalised estimated genus richness", limits = c(-1, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  guides(y.sec = guide_axis_manual(breaks = factor(sort(lats$mid))[c(2, 5)], 
                                   labels = c("S. Hemi.", "N. Hemi."),
                                   angle = 90, label_size = 16)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0)),
             axis.ticks.length.y.right = unit(0, "in")) +
  facet_wrap(~models, ncol = 1)
ggsave("./figures/diffs_heatmap_sqs.png", diffs_heatmap_sqs, width = 13, height = 18)

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
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  group_by(stage_bin, paleolat_bin, model.x, model.y) %>%
  summarise(diff = n_genera_norm.x - n_genera_norm.y, .groups = "drop") %>%
  complete(stage_bin = stages$bin, model.x, model.y) %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  mutate(models = paste(model.x, model.y, sep = " - ")) %>%
  left_join(stages, by = c("stage_bin" = "bin")) %>%
  left_join(lats, by = c("paleolat_bin" = "bin"))

diffs_heatmap_raw <- ggplot(data = diffs_raw) +
  geom_tile(aes(x = mid_ma, y = factor(mid),
                width = duration_myr, height = 1, fill = diff)) +
  geom_hline(yintercept = 3.5) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0)) +
  scale_y_discrete("Palaeolatitudinal bin",
                   limits = factor(sort(lats$mid)),
                   labels = c("High", "Middle", "Low", "Low", "Middle", "High")) +
  scale_fill_viridis_c("Difference in normalised raw genus richness", limits = c(-1, 1),
                       option = "plasma", end = .8, na.value = "grey80",
                       guide = guide_colorbar(barwidth = 15)) +
  coord_geo(list("bottom", "bottom"), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  guides(y.sec = guide_axis_manual(breaks = factor(sort(lats$mid))[c(2, 5)], 
                                   labels = c("S. Hemi.", "N. Hemi."),
                                   angle = 90, label_size = 16)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             legend.title = element_text(margin = margin(0, 15, 0, 0)),
             axis.ticks.length.y.right = unit(0, "in")) +
  facet_wrap(~models, ncol = 1)
ggsave("./figures/diffs_heatmap_raw.png", diffs_heatmap_raw, width = 13, height = 18)

