# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 09_summary_stats.R
# Aim: Generate summary stats for the manuscript
# Last updated: 2026-09-03
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(tidyverse)
time_bins <- read.csv("data/time_bins.csv")
lat_bins <- read.csv("data/lat_bins.csv")

# Palaeogeographic reconstructions --------------------------------------
# How many collections with palaeocoordinates?
df <- readRDS("data/processed/pbdb_data.RDS")
# Reduce to collections
df <- df[, c("bin_midpoint", "lng","lat", 
             "p_lng_MERDITH2021", "p_lat_MERDITH2021", "MERDITH2021_bin",
             "p_lng_PALEOMAP", "p_lat_PALEOMAP", "PALEOMAP_bin",
             "p_lng_GOLONKA", "p_lat_GOLONKA", "GOLONKA_bin",
             "p_lng_TorsvikCocks2017", "p_lat_TorsvikCocks2017", "TorsvikCocks2017_bin")]
df <- unique(df)
# Subset data and count (Lutetian)
lutetian <- subset(df, bin_midpoint == 44.55000)
nrow(lutetian) - sum(is.na(lutetian$PALEOMAP_bin))
nrow(lutetian) - sum(is.na(lutetian$GOLONKA_bin))
nrow(lutetian) - sum(is.na(lutetian$MERDITH2021_bin))
nrow(lutetian) - sum(is.na(lutetian$TorsvikCocks2017_bin))
# Subset data and count (Anisian)
anisian <- subset(df, bin_midpoint == 244.08200)
nrow(anisian) - sum(is.na(anisian$PALEOMAP_bin))
nrow(anisian) - sum(is.na(anisian$GOLONKA_bin))
nrow(anisian) - sum(is.na(anisian$MERDITH2021_bin))
nrow(anisian) - sum(is.na(anisian$TorsvikCocks2017_bin))

# How many available bins for estimating richness (non-na)?
df <- read.csv("data/processed/interpolations.csv")
df %<>% 
  group_by(stage, model) %>%
  tally(!is.na(qD)) %>%
  group_by(stage) %>%
  summarise(max = max(n), min = min(n), diff = abs(max(n) - min(n)))
# Counts
table(df$diff)
# Percentages
(table(df$diff) / nrow(df)) * 100

# Palaeogeographic differences per period
## plat
plat <- read.csv("results/plat_diff_degrees.csv")
plat <- plat[order(plat$time), ]
cenozoic <- subset(plat, time < 66)
mesozoic <- subset(plat, time <= 251.902 & time > 66)
palaeozoic <- subset(plat, time > 251.902)
mean(cenozoic$med)
mean(mesozoic$med)
mean(palaeozoic$med)
# How many fossil collections more than 25º palaeolatitudinal difference 
plat <- read.csv("results/coll_plat_diff.csv")
plat <- plat[, c("bin_midpoint", "collection_no", "median_pair_lat_diff",
             "p_lng_MERDITH2021", "p_lat_MERDITH2021", "MERDITH2021_bin",
             "p_lng_PALEOMAP", "p_lat_PALEOMAP", "PALEOMAP_bin",
             "p_lng_GOLONKA", "p_lat_GOLONKA", "GOLONKA_bin",
             "p_lng_TorsvikCocks2017", "p_lat_TorsvikCocks2017", "TorsvikCocks2017_bin")]
jurassic <- subset(plat, bin_midpoint <= 201.4 & bin_midpoint > 149.2000)
triassic <- subset(plat, bin_midpoint <= 251.902 & bin_midpoint > 201.4)
permian <- subset(plat, bin_midpoint <= 298.9 & bin_midpoint > 251.902)
carboniferous <- subset(plat, bin_midpoint <= 358.8600 & bin_midpoint > 298.9)
devonian <- subset(plat, bin_midpoint <= 419.6200 & bin_midpoint > 358.8600)
silurian <- subset(plat, bin_midpoint <= 443.100 & bin_midpoint > 419.6200)
ordovician <- subset(plat, bin_midpoint <= 486.850 & bin_midpoint > 443.100)
cambrian <- subset(plat, bin_midpoint <= 538.800 & bin_midpoint > 486.850)
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
maxlat <- read.csv("results/max_lat_sqs.csv")
# Create supplementary table
tbl <- maxlat %>%
  group_by(model, max_bin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = max_bin, values_from = count) %>%
  select(-`NA`) %>%
  rename_with(~c("Global Plate Model", "High (N)", "Middle (N)", "Low (N)", 
                 "Low (S)", "Middle (S)", "High (S)"))
write.csv(tbl, "results/bin_counts_tbl.csv", row.names = FALSE)

# Which stages are tropical peaks in diversity under any model?
n_peak_trop <- unique(maxlat$stage[which(maxlat$max_bin == 3)])
s_peak_trop <- unique(maxlat$stage[which(maxlat$max_bin == 4)])
# Which stages are temperate peaks in diversity under any model?
n_peak_tem <- unique(maxlat$stage[which(maxlat$max_bin == 2)])
s_peak_tem <- unique(maxlat$stage[which(maxlat$max_bin == 5)])
# Which stages are polar peaks in diversity under any model?
n_peak_pol <- unique(maxlat$stage[which(maxlat$max_bin == 1)])
s_peak_pol <- unique(maxlat$stage[which(maxlat$max_bin == 6)])
# How many stages match?
sum(!is.na(match(x = n_peak_trop, table = n_peak_tem)))
sum(!is.na(match(x = s_peak_trop, table = s_peak_tem)))
# Tropical and Polar peaks
sum(!is.na(match(x = n_peak_trop, table = n_peak_pol)))
sum(!is.na(match(x = s_peak_trop, table = s_peak_pol)))
# Tropical and Temperate peaks
sum(!is.na(match(x = n_peak_trop, table = n_peak_tem)))
sum(!is.na(match(x = s_peak_trop, table = s_peak_tem)))
    

# Which stages have the same peak in diversity?
df <- read.csv("data/processed/interpolations.csv")
df <- subset(df, !is.na(paleolat_bin))
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
  filter(!is.na(qD)) %>%
  summarise(n = which.max(qD)) %>%
  group_by(stage) %>%
  summarise(n = length(unique(n)))
length(which(n_df$n == 1)) / nrow(n_df)

s_df %<>% 
  group_by(stage, model) %>%
  filter(!is.na(qD)) %>%
  summarise(n = which.max(qD)) %>%
  group_by(stage) %>%
  summarise(n = length(unique(n)))
length(which(s_df$n == 1)) / nrow(s_df)

# Rank order differences
df <- read.csv("results/rank_order_sqs.csv")
# Palaeozoic
palaeozoic <- subset(df, stage <= 48)
max(palaeozoic$avg_norm, na.rm = TRUE)
mean(palaeozoic$avg_norm, na.rm = TRUE)
# Mesozoic
mesozoic <- subset(df, stage >= 49 & stage <= 78)
max(mesozoic$avg_norm, na.rm = TRUE)
mean(mesozoic$avg_norm, na.rm = TRUE)
#Cenozoic
cenozoic <- subset(df, stage >= 79)
max(cenozoic$avg_norm, na.rm = TRUE)
mean(cenozoic$avg_norm, na.rm = TRUE)
