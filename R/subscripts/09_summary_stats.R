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
plat <- readRDS("./results/plat_diff_degrees.RDS")
plat <- plat[order(plat$bin_midpoint), ]
cenozoic <- subset(plat, bin_midpoint < 66)
mesozoic <- subset(plat, bin_midpoint <= 251.9020 & bin_midpoint > 66)
palaeozoic <- subset(plat, bin_midpoint > 251.9020)
mean(cenozoic$med)
mean(mesozoic$med)
mean(palaeozoic$med)
## gcd
gcd <- readRDS("./results/plat_diff.RDS")
gcd <- gcd[order(gcd$time), ]
cenozoic <- subset(gcd, gcd$time < 66)
mesozoic <- subset(gcd, gcd$time <= 251.9020 & gcd$time > 66)
palaeozoic <- subset(gcd, gcd$time > 251.9020)
mean(cenozoic$med_latD)
mean(mesozoic$med_latD)
mean(palaeozoic$med_latD)
# How many fossil collections more than 25º palaeolatitudinal difference 
plat <- readRDS("./results/coll_plat_diff.RDS")
plat <- plat[order(plat$bin_midpoint), ]
jurassic <- subset(plat, bin_midpoint <= 201.3000 & bin_midpoint > 145.0000)
triassic <- subset(plat, bin_midpoint <= 251.9020 & bin_midpoint > 201.3000)
permian <- subset(plat, bin_midpoint <= 298.9000 & bin_midpoint > 251.9020)
carboniferous <- subset(plat, bin_midpoint <= 358.9000 & bin_midpoint > 298.9000)
devonian <- subset(plat, bin_midpoint <= 419.2000 & bin_midpoint > 358.9000)
silurian <- subset(plat, bin_midpoint <= 443.8000 & bin_midpoint > 419.2000)
ordovician <- subset(plat, bin_midpoint <= 485.4000 & bin_midpoint > 443.8000)
cambrian <- subset(plat, bin_midpoint > 485.4000)
val <- 25
length(which(jurassic$median_pair_lat_diff > val))
length(which(triassic$median_pair_lat_diff > val))
length(which(permian$median_pair_lat_diff > val))
length(which(carboniferous$median_pair_lat_diff > val))
length(which(devonian$median_pair_lat_diff > val))
length(which(silurian$median_pair_lat_diff > val))
length(which(ordovician$median_pair_lat_diff > val))
length(which(cambrian$median_pair_lat_diff > val))

# Latitudinal biodiversity gradient reconstructions ---------------------
maxlat <- readRDS("./results/max_lat_sqs.RDS")
# Which stages are tropical peaks in diversity under any model?
n_peak_trop <- unique(maxlat$stage[which(maxlat$mid == "9.74")])
s_peak_trop <- unique(maxlat$stage[which(maxlat$mid == "-9.74")])
# Which stages are temperate peaks in diversity under any model?
n_peak_tem <- unique(maxlat$stage[which(maxlat$mid == "30.655")])
s_peak_tem <- unique(maxlat$stage[which(maxlat$mid == "-30.655")])
# Which stages are polar peaks in diversity under any model?
n_peak_pol <- unique(maxlat$stage[which(maxlat$mid == "65.915")])
s_peak_pol <- unique(maxlat$stage[which(maxlat$mid == "-65.915")])
# How many stages match?
sum(!is.na(match(x = n_peak_trop, table = n_peak_tem)))
sum(!is.na(match(x = s_peak_trop, table = s_peak_tem)))
# Which stages have the same peak in diversity?
df <- readRDS(file = "./data/processed/interpolations.RDS")
n_df <- subset(df, paleolat_bin <= 3)
s_df <- subset(df, paleolat_bin >= 4)
df %<>% 
  group_by(stage, model) %>%
  summarise(n = which.max(qD)) %>%
  group_by(stage) %>%
  summarise(n = length(unique(n)))
# Palaeozoic
palaeozoic <- subset(df, stage <= 48)
length(which(palaeozoic$n == 1)) / nrow(palaeozoic)
# Mesozoic
mesozoic <- subset(df, stage >= 49 & stage <= 78)
length(which(mesozoic$n == 1)) / nrow(mesozoic)
#Cenozoic
cenozoic <- subset(df, stage >= 79)
length(which(cenozoic$n == 1)) / nrow(cenozoic)
# per hemisphere
n_df %<>% 
  group_by(stage, model) %>%
  summarise(n = which.max(qD)) %>%
  group_by(stage) %>%
  summarise(n = length(unique(n)))
length(which(n_df$n == 1)) / nrow(n_df)

s_df %<>% 
  group_by(stage, model) %>%
  summarise(n = which.max(qD)) %>%
  group_by(stage) %>%
  summarise(n = length(unique(n)))
length(which(s_df$n == 1)) / nrow(s_df)
