# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: data_for_div.R
# Aim: Data processing for diversity figures
# Last updated: 2026-09-02
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(dplyr)
library(tidyr)

# Data ------------------------------------------------------------------
# Load bins
lats <- read.csv("data/lat_bins.csv")
stages <- read.csv("data/time_bins.csv")

# Load data
div_sqs <- read.csv("data/processed/interpolations.csv")
div_raw <- read.csv("data/processed/genus_counts.csv")

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