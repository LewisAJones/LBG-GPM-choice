# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_Diversity_and_Difference_heatmaps.R
# Aim: Data processing for 05_Diversity_figure.R and
#      06_Diversity_and_Difference_heatmaps.R
# Last updated: 2023-12-01
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(dplyr)
library(tidyr)

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
  complete(stage_bin, paleolat_bin, model) %>%
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
  complete(stage, paleolat_bin, model) %>%
  left_join(stages, by = c("stage" = "bin")) %>%
  left_join(lats, by = c("paleolat_bin" = "bin"))