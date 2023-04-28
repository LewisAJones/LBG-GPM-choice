# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 05_metrics_figures.R
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
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
theme_will <- function(...) {
  theme(axis.ticks = element_line(color = "black", linewidth = 1),
        axis.line = element_blank(),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(linetype = "solid", colour = "black",
                                    fill = NA, linewidth = 2),
        ...)
}
GTS2020_periods <- time_bins(rank = "period") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Metric #1: max diversity bin  -----------------------------------------
met1_sqs <- div_sqs %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  mutate(hemi = ifelse(max > 0, "north", "south")) %>%
  group_by(model, stage, hemi) %>%
  filter(sum(!is.na(qD)) > 0) %>%
  summarise(max_bin = paleolat_bin[which.max(qD)], .groups = "drop") %>%
  complete(model, stage, hemi) %>%
  left_join(lat_bins, by = c("max_bin" = "bin")) %>%
  left_join(time_bins, by = c("stage" = "bin"))

met1_sqs_rects <- met1_sqs %>%
  group_by(stage, hemi, max_ma, min_ma) %>%
  summarise(all_same = length(unique(max_bin)) == 1, .groups = "drop") %>%
  filter(all_same) %>%
  mutate(ymin = ifelse(hemi == "north", 0, -Inf),
         ymax = ifelse(hemi == "north", Inf, 0))

gg_met1_sqs <- ggplot(met1_sqs) +
  geom_rect(data = met1_sqs_rects, aes(xmin = max_ma, xmax = min_ma,
                                       ymin = ymin, ymax = ymax), fill = "grey90") +
  geom_line(aes(x = mid_ma, y = mid, color = model,
                group = interaction(hemi, model)), linewidth = .75) +
  geom_hline(yintercept = 0) +
  annotate(geom = "text", x = 538, y = c(-77, 77),
           label = c("South. Hemisphere", "North. Hemisphere"),
           hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous(expression('Most Diverse Bin ('*degree*')'),
                     breaks = seq(-80, 80, 20), limits = c(-80, 80)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_1.pdf", gg_met1_sqs, width = 13, height = 4.5)

met1_raw <- div_raw %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  mutate(hemi = ifelse(paleolat_bin <= 6, "north", "south")) %>%
  group_by(model, stage_bin, hemi) %>%
  filter(sum(!is.na(n_genera)) > 0) %>%
  summarise(max_bin = paleolat_bin[which.max(n_genera)], .groups = "drop") %>%
  complete(model, stage_bin, hemi) %>%
  left_join(lat_bins, by = c("max_bin" = "bin")) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

met1_raw_rects <- met1_raw %>%
  group_by(stage_bin, hemi, max_ma, min_ma) %>%
  summarise(all_same = length(unique(max_bin)) == 1, .groups = "drop") %>%
  filter(all_same) %>%
  mutate(ymin = ifelse(hemi == "north", 0, -Inf),
         ymax = ifelse(hemi == "north", Inf, 0))

gg_met1_raw <- ggplot(met1_raw) +
  geom_rect(data = met1_raw_rects, aes(xmin = max_ma, xmax = min_ma,
                                       ymin = ymin, ymax = ymax), fill = "grey90") +
  geom_line(aes(x = mid_ma, y = mid, color = model,
                group = interaction(hemi, model)), linewidth = .75) +
  #geom_point(aes(x = mid_ma, y = mid, color = model, shape = model,
  #               group = interaction(hemi, model))) +
  geom_hline(yintercept = 0) +
  annotate(geom = "text", x = 538, y = c(-77, 77),
           label = c("South. Hemisphere", "North. Hemisphere"),
           hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous(expression('Most Diverse Bin ('*degree*')'),
                     breaks = seq(-80, 80, 20), limits = c(-80, 80)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_1_raw.pdf", gg_met1_raw, width = 13, height = 4.5)

# Metric #2: rank order diffs ---------------------------------------
# what do we do when one model has data for a bin but the other doesn't?
# if we don't include that bin, then the sum will be biased to be smaller
met2_sqs <- div_sqs %>%
  select(stage, paleolat_bin, model, qD) %>%
  group_by(model, stage) %>%
  arrange(-qD) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(stage, paleolat_bin), names_from = model, values_from = rank) %>%
  group_by(stage) %>%
  summarise(p_g = sum(abs(PALEOMAP - GOLONKA)),
            p_m = sum(abs(PALEOMAP - MERDITH2021)),
            g_m = sum(abs(GOLONKA - MERDITH2021)), .groups = "drop") %>%
  pivot_longer(cols = c(p_g, p_m, g_m)) %>%
  left_join(time_bins, by = c("stage" = "bin")) %>%
  complete(stage, name)

gg_met2_sqs <- ggplot(met2_sqs) +
  geom_line(aes(x = mid_ma, y = value, color = name, group = name), linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Sum of Rank Order Differences") +
  scale_colour_viridis_d(NULL, end = .9, labels = c("GOLONKA/MERDITH2021",
                                                    "PALEOMAP/GOLONKA",
                                                    "PALEOMAP/MERDITH2021")) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_2.pdf", gg_met2_sqs, width = 13, height = 4.5)

met2_raw <- div_raw %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  group_by(model, stage_bin) %>%
  arrange(-n_genera) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(stage_bin, paleolat_bin), names_from = model, values_from = rank) %>%
  group_by(stage_bin) %>%
  summarise(p_g = sum(abs(PALEOMAP - GOLONKA)),
            p_m = sum(abs(PALEOMAP - MERDITH2021)),
            g_m = sum(abs(GOLONKA - MERDITH2021)), .groups = "drop") %>%
  pivot_longer(cols = c(p_g, p_m, g_m)) %>%
  left_join(time_bins, by = c("stage_bin" = "bin")) %>%
  complete(stage_bin, name)

gg_met2_raw <- ggplot(met2_raw) +
  geom_line(aes(x = mid_ma, y = value, color = name, group = name), linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Sum of Rank Order Differences") +
  scale_colour_viridis_d(NULL, end = .9, labels = c("GOLONKA/MERDITH2021",
                                                    "PALEOMAP/GOLONKA",
                                                    "PALEOMAP/MERDITH2021")) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_2_raw.pdf", gg_met2_raw, width = 13, height = 4.5)

# Metric #3: sum of squares -----------------------------------------
met3_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  pivot_wider(id_cols = c(stage, paleolat_bin), names_from = model, values_from = qD) %>%
  group_by(stage) %>%
  summarise(pg_sum = sum(abs(PALEOMAP - GOLONKA), na.rm = TRUE),
            pg_n = sum(!is.na(PALEOMAP) & !is.na(GOLONKA)),
            pm_sum = sum(abs(PALEOMAP - MERDITH2021), na.rm = TRUE),
            pm_n = sum(!is.na(PALEOMAP) & !is.na(MERDITH2021)),
            gm_sum = sum(abs(GOLONKA - MERDITH2021), na.rm = TRUE),
            gm_n = sum(!is.na(GOLONKA) & !is.na(MERDITH2021)),
            .groups = "drop") %>%
  pivot_longer(cols = c(pg_sum, pg_n, pm_sum, pm_n, gm_sum, gm_n)) %>%
  separate_wider_delim(name, "_", names = c("models", "name")) %>%
  pivot_wider(id_cols = c(stage, models)) %>%
  mutate(avg = sqrt(sum)/n) %>%
  left_join(time_bins, by = c("stage" = "bin"))

gg_met3_sqs <- ggplot(met3_sqs) +
  geom_line(aes(x = mid_ma, y = avg, color = models, group = models), linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Average Diversity Estimate Difference", trans = "log10") +
  scale_colour_viridis_d(NULL, end = .9, labels = c("GOLONKA/MERDITH2021",
                                                    "PALEOMAP/GOLONKA",
                                                    "PALEOMAP/MERDITH2021")) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_3.pdf", gg_met3_sqs, width = 13, height = 4.5)

met3_raw <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  pivot_wider(id_cols = c(stage_bin, paleolat_bin), names_from = model, values_from = n_genera) %>%
  group_by(stage_bin) %>%
  summarise(pg_sum = sum(abs(PALEOMAP - GOLONKA), na.rm = TRUE),
            pg_n = sum(!is.na(PALEOMAP) & !is.na(GOLONKA)),
            pm_sum = sum(abs(PALEOMAP - MERDITH2021), na.rm = TRUE),
            pm_n = sum(!is.na(PALEOMAP) & !is.na(MERDITH2021)),
            gm_sum = sum(abs(GOLONKA - MERDITH2021), na.rm = TRUE),
            gm_n = sum(!is.na(GOLONKA) & !is.na(MERDITH2021)),
            .groups = "drop") %>%
  pivot_longer(cols = c(pg_sum, pg_n, pm_sum, pm_n, gm_sum, gm_n)) %>%
  separate_wider_delim(name, "_", names = c("models", "name")) %>%
  pivot_wider(id_cols = c(stage_bin, models)) %>%
  mutate(avg = sqrt(sum)/n) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

gg_met3_raw <- ggplot(met3_raw) +
  geom_line(aes(x = mid_ma, y = avg, color = models, group = models), linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Average Raw Diversity Difference", trans = "log10") +
  scale_colour_viridis_d(NULL, end = .9, labels = c("GOLONKA/MERDITH2021",
                                                    "PALEOMAP/GOLONKA",
                                                    "PALEOMAP/MERDITH2021")) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_3_raw.pdf", gg_met3_raw, width = 13, height = 4.5)

# Metric #4: best fit model -----------------------------------------
met4_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  group_by(stage, model) %>%
  summarise(poly = which.min(sapply(1:5, function(i) {
    tryCatch(AICc(glm(qD ~ poly(mid, i))),
             error = function(e) NA)
  }))) %>%
  left_join(time_bins, by = c("stage" = "bin"))

gg_met4_sqs <- ggplot(met4_sqs) +
  geom_tile(aes(x = mid_ma, y = model, width = max_ma - min_ma, height = 1, fill = as.factor(poly))) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0)) +
  scale_y_discrete(NULL) +
  scale_fill_viridis_d(NULL, end = .9) +
  coord_geo(dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_4.pdf", gg_met4_sqs, width = 13, height = 4.5)

met4_raw <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  group_by(stage_bin, model) %>%
  summarise(poly = which.min(sapply(1:5, function(i) {
    tryCatch(AICc(glm(n_genera ~ poly(mid, i))),
             error = function(e) NA)
  }))) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

gg_met4_raw <- ggplot(met4_raw) +
  geom_tile(aes(x = mid_ma, y = model, width = max_ma - min_ma, height = 1, fill = as.factor(poly))) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0)) +
  scale_y_discrete(NULL) +
  scale_fill_viridis_d(NULL, end = .9) +
  coord_geo(dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_4_raw.pdf", gg_met4_raw, width = 13, height = 4.5)

# Metric #4b: GAM knots -----------------------------------------
# met4b_sqs <- div_sqs %>%
#   filter(!is.na(paleolat_bin)) %>%
#   select(stage, paleolat_bin, model, qD) %>%
#   left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
#   group_by(stage, model) %>%
#   summarise(edf = 
#     tryCatch(summary(gam(qD ~ s(mid)))$edf,
#              error = function(e) NA)
#   ) %>%
#   left_join(time_bins, by = c("stage" = "bin"))
# 
# gg_met4b_sqs <- ggplot(met4b_sqs) +
#   geom_tile(aes(x = mid_ma, y = model, width = max_ma - min_ma, height = 1, fill = edf)) +
#   scale_x_reverse("Time (Ma)", limits = c(541, 0)) +
#   scale_y_discrete(NULL) +
#   scale_fill_viridis_c(NULL, end = .9) +
#   coord_geo(dat = GTS2020_periods, lwd = 1,
#             bord = c("left", "right", "bottom")) +
#   theme_classic(base_size = 14) +
#   theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
# ggsave("./figures/metric_4b.pdf", gg_met4b_sqs, width = 13, height = 4.5)
