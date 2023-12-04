# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 09_summary_stats.R
# Aim: Generate summary stats for the manuscript
# Last updated: 2023-12-04
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(tidyverse)

# Palaeogeographic reconstructions --------------------------------------
# How many available bins for estimating richness (non-na)?
df <- readRDS(file = "./data/processed/interpolations.RDS")
df %<>% 
  group_by(stage, model) %>%
  tally(!is.na(qD)) %>%
  group_by(stage) %>%
  summarise(max = max(n), min = min(n), diff = abs(max(n) - min(n)))
# Counts
table(df$diff)
# Percentages
(table(df$diff) / nrow(df)) * 100

# How many collections with palaeocoordinates?
df <- readRDS("./data/processed/pbdb_data.RDS")
# Reduce to collections
df <- df[, c("collection_no", "bin_midpoint",
             "p_lng_MERDITH2021", "p_lat_MERDITH2021", "MERDITH2021_bin",
             "p_lng_PALEOMAP", "p_lat_PALEOMAP", "PALEOMAP_bin",
             "p_lng_GOLONKA", "p_lat_GOLONKA", "GOLONKA_bin")]
df <- unique(df)
# Subset data and count (Tortonian)
tortonian <- subset(df, bin_midpoint == 9.438)
nrow(tortonian) - sum(is.na(tortonian$PALEOMAP_bin))
nrow(tortonian) - sum(is.na(tortonian$GOLONKA_bin))
nrow(tortonian) - sum(is.na(tortonian$MERDITH2021_bin))
# Subset data and count (Anisian)
anisian <- subset(df, bin_midpoint == 244.60000)
nrow(anisian) - sum(is.na(anisian$PALEOMAP_bin))
nrow(anisian) - sum(is.na(anisian$GOLONKA_bin))
nrow(anisian) - sum(is.na(anisian$MERDITH2021_bin))

# Palaeogeographic differences per period
## plat
plat <- readRDS("./results/plat_diff.RDS")
plat <- plat[order(plat$bin_midpoint), ]
cenozoic <- subset(plat, bin_midpoint < 66)
mesozoic <- subset(plat, bin_midpoint <= 251.9020 & bin_midpoint > 66)
palaeozoic <- subset(plat, bin_midpoint > 251.9020)
mean(cenozoic$med)
mean(mesozoic$med)
mean(palaeozoic$med)
## geodes
geodes <- readRDS("./results/geodes_diff.RDS")
geodes <- geodes[order(geodes$time), ]
cenozoic <- subset(geodes, geodes$time < 66)
mesozoic <- subset(geodes, geodes$time <= 251.9020 & geodes$time > 66)
palaeozoic <- subset(geodes, geodes$time > 251.9020)
mean(cenozoic$GDD)
mean(mesozoic$GDD)
mean(palaeozoic$GDD)

# Latitudinal biodiversity gradient reconstructions ---------------------



