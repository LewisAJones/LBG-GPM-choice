# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 07_metrics_figures.R
# Aim: Assess and plot max diversity bin through time per model (Fig.5 and S4)
#      and normalised averaged rank order differences in reconstructed palaeodiversity
#      between GPMs through time (Fig.6 and S4).
# Last updated: 2024-12-17
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(dplyr)
library(grid)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(broom)
library(deeptime)
library(palaeoverse)
library(AICcmodavg)
library(mgcv)
source("./R/options.R")

# Load data -------------------------------------------------------------
div_sqs <- readRDS("./data/processed/interpolations.RDS")
div_raw <- readRDS("./data/processed/genus_counts.RDS")
lat_bins <- readRDS("./data/lat_bins.RDS")
time_bins <- readRDS("./data/time_bins.RDS")

# Setup common things for figures
source("./R/functions/theme_will.R")
ics_periods <- time_bins(scale = "international periods") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
ics_eras <- time_bins(scale = "international eras") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Metric #1: max diversity bin  -----------------------------------------
met1_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  mutate(hemi = ifelse(mid > 0, "north", "south")) %>%
  group_by(model, stage, hemi) %>%
  filter(sum(!is.na(qD)) > 0) %>%
  # find lat bin with highest diversity
  summarise(max_bin = paleolat_bin[which.max(qD)], .groups = "drop") %>%
  complete(model, stage = time_bins$bin, hemi) %>%
  left_join(lat_bins, by = c("max_bin" = "bin")) %>%
  left_join(time_bins, by = c("stage" = "bin"))

met1_sqs_rects <- met1_sqs %>%
  filter(!is.na(max_bin)) %>%
  group_by(stage, hemi, max_ma, min_ma) %>%
  # do all models have the same bin?
  summarise(all_same = length(unique(max_bin)) == 1, .groups = "drop")

met1_sqs_rects %>%
  mutate(era = cut(max_ma, c(ics_eras$max_age[1], ics_eras$min_age), rev(ics_eras$name))) %>%
  group_by(era, hemi, all_same) %>%
  count()

met1_sqs_rects <- met1_sqs_rects %>%
  filter(all_same) %>%
  mutate(ymin = ifelse(hemi == "north", 3.5, -Inf),
         ymax = ifelse(hemi == "north", Inf, 3.5))
saveRDS(met1_sqs_rects, "./results/max_lat_bin_sqs.RDS")

gg_met1_sqs <- ggplot(met1_sqs, aes(x = mid_ma, y = as.numeric(factor(mid)), color = model)) +
  geom_rect(data = met1_sqs_rects, inherit.aes = FALSE,
            aes(xmin = max_ma, xmax = min_ma, ymin = ymin, ymax = ymax), fill = "grey90") +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(aes(group = interaction(hemi, model)), linewidth = .75,
            position = position_dodge(width = 2)) +
  geom_hline(yintercept = 3.5) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0), expand = expansion()) +
  scale_y_continuous("Most diverse bin", breaks = 1:6,
                     labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                     expand = expansion(add = .5)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, ylim = c(1, 6), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  guides(y.sec = guide_axis_manual(breaks = c(2, 5), 
                                   labels = c("S. Hemi.", "N. Hemi."),
                                   angle = 90, label_size = 16)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             axis.ticks.length.y.right = unit(0, "in"))
ggsave("./figures/metric_1.png", gg_met1_sqs, width = 13, height = 6)
saveRDS(met1_sqs, "./results/max_lat_sqs.RDS")

met1_raw <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  mutate(hemi = ifelse(mid > 0, "north", "south")) %>%
  group_by(model, stage_bin, hemi) %>%
  filter(sum(!is.na(n_genera)) > 0) %>%
  # find lat bin with highest diversity
  summarise(max_bin = paleolat_bin[which.max(n_genera)], .groups = "drop") %>%
  complete(model, stage_bin = time_bins$bin, hemi) %>%
  left_join(lat_bins, by = c("max_bin" = "bin")) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

met1_raw_rects <- met1_raw %>%
  filter(!is.na(max_bin)) %>%
  group_by(stage_bin, hemi, max_ma, min_ma) %>%
  # do all models have the same bin?
  summarise(all_same = length(unique(max_bin)) == 1, .groups = "drop")

met1_raw_rects %>%
  mutate(era = cut(max_ma, c(ics_eras$max_age[1], ics_eras$min_age), rev(ics_eras$name))) %>%
  group_by(era, hemi, all_same) %>%
  count()

met1_raw_rects <- met1_raw_rects %>%
  filter(all_same) %>%
  mutate(ymin = ifelse(hemi == "north", 3.5, -Inf),
         ymax = ifelse(hemi == "north", Inf, 3.5))

gg_met1_raw <- ggplot(met1_raw, aes(x = mid_ma, y = as.numeric(factor(mid)), color = model)) +
  geom_rect(data = met1_raw_rects, inherit.aes = FALSE,
            aes(xmin = max_ma, xmax = min_ma, ymin = ymin, ymax = ymax), fill = "grey90") +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(aes(group = interaction(hemi, model)), linewidth = .75,
            position = position_dodge(width = 2)) +
  geom_hline(yintercept = 3.5) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0), expand = expansion()) +
  scale_y_continuous("Most diverse bin", breaks = 1:6,
                     labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                     expand = expansion(add = .5)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, ylim = c(1, 6), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  guides(y.sec = guide_axis_manual(breaks = c(2, 5), 
                                   labels = c("S. Hemi.", "N. Hemi."),
                                   angle = 90, label_size = 16)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5),
             axis.ticks.length.y.right = unit(0, "in"))
ggsave("./figures/metric_1_raw.png", gg_met1_raw, width = 13, height = 6)

# Metric #2: rank order diffs ---------------------------------------
met2_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin), !is.na(qD)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  inner_join(., ., by = c("stage" = "stage", "paleolat_bin" = "paleolat_bin"),
             relationship = "many-to-many") %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  group_by(stage, model.x, model.y) %>%
  # get ranks for "x" model
  arrange(-qD.x) %>%
  mutate(rank.x = row_number()) %>%
  # get ranks for "y" model
  arrange(-qD.y) %>%
  mutate(rank.y = row_number()) %>%
  # calculate average rank difference
  summarise(avg = mean(abs(rank.x - rank.y)), n_bins = n(), .groups = "drop") %>%
  # remove time bins with 1 or fewer lat bins
  filter(n_bins > 1) %>%
  # normalize to maximum possible average difference
  mutate(avg_norm = avg / c(0, 1, 4/3, 2, 2.4, 3)[n_bins]) %>%
  complete(stage = time_bins$bin, model.x, model.y) %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  mutate(models = paste(model.x, model.y, sep = "/")) %>%
  left_join(time_bins, by = c("stage" = "bin"))

saveRDS(object = met2_sqs, file = "results/rank_order_sqs.RDS")

gg_met2_sqs <- ggplot(met2_sqs, aes(x = mid_ma, y = avg_norm, group = models)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0), expand = expansion()) +
  scale_y_continuous("Norm. avg. rank order diff.", limits = c(0, 1)) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5)) +
  facet_wrap(vars(models), ncol = 1)
ggsave("./figures/metric_2.png", gg_met2_sqs, width = 13, height = 18)

met2_sqs %>%
  nest_by(models) %>%
  mutate(mod = list(lm(avg_norm ~ mid_ma, data = data))) %>% 
  summarize(glance(mod))

ggsave("./figures/metric_2_reg.png",
       gg_met2_sqs + geom_smooth(method = "lm"),
       width = 13, height = 6)

met2_raw <- div_raw %>%
  filter(!is.na(paleolat_bin), !is.na(n_genera), n_genera > 0) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  inner_join(., ., by = c("stage_bin" = "stage_bin", "paleolat_bin" = "paleolat_bin"),
             relationship = "many-to-many") %>%
  # remove duplicates (with reversed x and y)
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  group_by(stage_bin, model.x, model.y) %>%
  # get ranks for "x" model
  arrange(-n_genera.x) %>%
  mutate(rank.x = row_number()) %>%
  # get ranks for "y" model
  arrange(-n_genera.y) %>%
  mutate(rank.y = row_number()) %>%
  # calculate average rank difference
  summarise(avg = mean(abs(rank.x - rank.y)), n_bins = n(), .groups = "drop") %>%
  # remove time bins with 1 or fewer lat bins
  filter(n_bins > 1) %>%
  # normalize to maximum possible average difference
  mutate(avg_norm = avg / c(0, 1, 4/3, 2, 2.4, 3)[n_bins]) %>%
  complete(stage_bin = time_bins$bin, model.x, model.y) %>%
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "MERDITH2021") |
           (model.x == "TorsvikCocks2017" & model.y == "PALEOMAP") |
           (model.x == "TorsvikCocks2017" & model.y == "GOLONKA")) %>%
  mutate(models = paste(model.x, model.y, sep = "/")) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

gg_met2_raw <- ggplot(met2_raw, aes(x = mid_ma, y = avg_norm, group = models)) +
  geom_point(size = 1.5) +
  geom_line(linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(538.800, 0), expand = expansion()) +
  scale_y_continuous("Norm. avg. rank order diff.", limits = c(0, 1)) +
  coord_geo(list("bottom", "bottom"), expand = TRUE, dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5)) +
  facet_wrap(vars(models), ncol = 1)
ggsave("./figures/metric_2_raw.png", gg_met2_raw, width = 13, height = 6)

met2_raw %>%
  nest_by(models) %>%
  mutate(mod = list(lm(avg_norm ~ mid_ma, data = data))) %>% 
  summarize(glance(mod))

ggsave("./figures/metric_2_raw_reg.png",
       gg_met2_raw + geom_smooth(method = "lm"),
       width = 13, height = 6)

# Combine #1 and #2 -------------------------------------------------
ggsave("./figures/metrics_1_2_sqs.png",
       ggarrange2(gg_met1_sqs, gg_met2_sqs, ncol = 1, draw = FALSE,
                  labels = c("(a)", "(b)"), label.args = list(gp = gpar(font = 2, cex = 2))),
       width = 13, height = 12)

ggsave("./figures/metrics_1_2_raw.png",
       ggarrange2(gg_met1_raw, gg_met2_raw, ncol = 1, draw = FALSE,
                  labels = c("(a)", "(b)"), label.args = list(gp = gpar(font = 2, cex = 2))),
       width = 13, height = 12)
