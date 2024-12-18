# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 09_summary_stats.R
# Aim: Generate summary stats for the manuscript
# Last updated: 2024-12-17
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(tidyverse)
time_bins <- readRDS("data/time_bins.RDS")
lat_bins <- readRDS("data/lat_bins.RDS")

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
lutetian <- subset(df, bin_midpoint == 44.5)
nrow(lutetian) - sum(is.na(lutetian$PALEOMAP_bin))
nrow(lutetian) - sum(is.na(lutetian$GOLONKA_bin))
nrow(lutetian) - sum(is.na(lutetian$MERDITH2021_bin))
nrow(lutetian) - sum(is.na(lutetian$TorsvikCocks2017_bin))
# Subset data and count (Anisian)
anisian <- subset(df, bin_midpoint == 244.6)
nrow(anisian) - sum(is.na(anisian$PALEOMAP_bin))
nrow(anisian) - sum(is.na(anisian$GOLONKA_bin))
nrow(anisian) - sum(is.na(anisian$MERDITH2021_bin))
nrow(anisian) - sum(is.na(anisian$TorsvikCocks2017_bin))

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

# Palaeogeographic differences per period
## plat
plat <- readRDS("./results/plat_diff_degrees.RDS")
plat <- plat[order(plat$time), ]
cenozoic <- subset(plat, time < 66)
mesozoic <- subset(plat, time <= 251.902 & time > 66)
palaeozoic <- subset(plat, time > 251.902)
mean(cenozoic$med)
mean(mesozoic$med)
mean(palaeozoic$med)
## gcd
gcd <- readRDS("./results/plat_diff_km.RDS")
gcd <- gcd[order(gcd$time), ]
cenozoic <- subset(gcd, time < 66)
mesozoic <- subset(gcd, time <= 251.902 & time > 66)
palaeozoic <- subset(gcd, time > 251.902)
mean(cenozoic$med)
mean(mesozoic$med)
mean(palaeozoic$med)
# How many fossil collections more than 25º palaeolatitudinal difference 
plat <- readRDS("./results/coll_plat_diff.RDS")
plat <- plat[, c("bin_midpoint", "collection_no", "median_pair_lat_diff",
             "p_lng_MERDITH2021", "p_lat_MERDITH2021", "MERDITH2021_bin",
             "p_lng_PALEOMAP", "p_lat_PALEOMAP", "PALEOMAP_bin",
             "p_lng_GOLONKA", "p_lat_GOLONKA", "GOLONKA_bin",
             "p_lng_TorsvikCocks2017", "p_lat_TorsvikCocks2017", "TorsvikCocks2017_bin")]
jurassic <- subset(plat, bin_midpoint <= 201.4 & bin_midpoint > 145)
triassic <- subset(plat, bin_midpoint <= 251.902 & bin_midpoint > 201.4)
permian <- subset(plat, bin_midpoint <= 298.9 & bin_midpoint > 251.902)
carboniferous <- subset(plat, bin_midpoint <= 358.9 & bin_midpoint > 298.9)
devonian <- subset(plat, bin_midpoint <= 419.2 & bin_midpoint > 358.9)
silurian <- subset(plat, bin_midpoint <= 443.8 & bin_midpoint > 419.2)
ordovician <- subset(plat, bin_midpoint <= 485.4 & bin_midpoint > 443.8)
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
# Create supplementary table
tbl <- maxlat %>%
  group_by(model, max_bin) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = max_bin, values_from = count) %>%
  select(-`NA`) %>%
  rename_with(~c("Global Plate Model", "High (N)", "Middle (N)", "Low (N)", 
                 "Low (S)", "Middle (S)", "High (S)"))
saveRDS(tbl, "results/bin_counts_tbl.RDS")

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
df <- readRDS(file = "./data/processed/interpolations.RDS")
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

# Rank order differences
df <- readRDS("./results/rank_order_sqs.RDS")
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
