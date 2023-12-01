# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 05b_additional_metrics_figures.R
# Aim: Assess and display additional proxy metrics for the GPM differences
#      on reconstructed LBGs, that are sum of the squares of the estimated 
#      diversities between models (metric 3) and best-fit model for the
#      diversity curves between flat, unimodal and bimodal.
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(dplyr)
library(grid)
library(tidyr)
library(ggplot2)
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
  mutate(avg = sum/n) %>%
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
ggsave("./figures/metric_3.png", gg_met3_sqs, width = 13, height = 4.5)

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
ggsave("./figures/metric_3_raw.png", gg_met3_raw, width = 13, height = 4.5)

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
ggsave("./figures/metric_4.png", gg_met4_sqs, width = 13, height = 4.5)

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
ggsave("./figures/metric_4_raw.png", gg_met4_raw, width = 13, height = 4.5)

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
# ggsave("./figures/metric_4b.png", gg_met4b_sqs, width = 13, height = 4.5)
