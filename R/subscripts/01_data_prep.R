# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 01_data_prep.R
# Aim: Prepare fossil occurrence data for downstream analyses.
# Last updated: 2024-12-16
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Load libraries --------------------------------------------------------
library(dplyr)
library(palaeoverse)
# Load project options
source("./R/options.R")
# Get data --------------------------------------------------------------
if (params$download || !file.exists("./data/raw/pbdb_data.RDS")) {
  # Use for fresh downloads
  library(RCurl)
  RCurl::curlSetOpt(3000)
  # Read data 
  occdf <- RCurl::getURL(url = params$pbdb_api, ssl.verifypeer = FALSE)
  occdf <- read.csv(textConnection(occdf))
  #save raw data
  saveRDS(occdf, "./data/raw/pbdb_data.RDS")
} else {
  occdf <- readRDS("./data/raw/pbdb_data.RDS")
}
# Process data ----------------------------------------------------------
# Collapse subgenera (remove characters from space onwards)
if (params$collapse_subgenera == TRUE) {
  occdf$genus <- sub(" .*", "", occdf$genus)
}
# Round off coordinates to stack collections 
occdf$lng <- round(occdf$lng, digits = params$n_decs)
occdf$lat <- round(occdf$lat, digits = params$n_decs)
# Bin data by collection for faster processing
colldf <- occdf[, c("collection_no", "lng", "lat",
                    "early_interval", "late_interval",
                    "max_ma", "min_ma")]
colldf <- unique(colldf)

## Set up time bins -----------------------------------------------------
bins <- time_bins(interval = "Phanerozoic",
                  rank = params$rank,
                  scale = params$GTS)
# Collapse Holocene equivalent bins
vec <- which(bins$interval_name == "Greenlandian")
bins$interval_name[vec] <- "Holocene"
bins$abbr[vec] <- "H"
# Update min_ma
bins$min_ma[vec] <- 0.0000
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Drop rows
bins <- bins[-which(bins$interval_name %in% c("Meghalayan", "Northgrippian")), ]
# Collapse Pleistocene equivalent bins
# Drop bins
pleis <- c("Late Pleistocene", "Chibanian", "Calabrian")
bins <- bins[-which(bins$interval_name %in% pleis), ]
# update Gelasian to be all of the Pleistocene
vec <- which(bins$interval_name == "Gelasian")
bins$interval_name[vec] <- "Pleistocene"
# Update min_ma
bins$min_ma[vec] <- bins[which(bins$interval_name == "Holocene"), "max_ma"]
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Update bin numbers
bins$bin <- 1:nrow(bins)
row.names(bins) <- 1:nrow(bins)
# Save time bins
saveRDS(object = bins, file = "./data/time_bins.RDS")

## Clean up age assignments ---------------------------------------------
# Remove any regular prefixes
# Early
colldf$early_interval <- gsub(pattern = "Early ",
                              replacement = "",
                              x = colldf$early_interval)
colldf$late_interval <- gsub(pattern = "Early ",
                              replacement = "",
                              x = colldf$late_interval)
# Middle
colldf$early_interval <- gsub(pattern = "Middle ",
                              replacement = "",
                              x = colldf$early_interval)
colldf$late_interval <- gsub(pattern = "Middle ",
                             replacement = "",
                             x = colldf$late_interval)
# Late
colldf$early_interval <- gsub(pattern = "Late ",
                              replacement = "",
                              x = colldf$early_interval)
colldf$late_interval <- gsub(pattern = "Late ",
                             replacement = "",
                             x = colldf$late_interval)
# Look up ages for intervals names using GTS2023/08
bins$early_stage <- bins$interval_name
bins$late_stage <- bins$interval_name
# Add ages
colldf <- look_up(occdf = colldf,
                  early_interval = "early_interval",
                  late_interval = "late_interval",
                  int_key = bins,
                  assign_with_GTS = FALSE)
# Which intervals could not be looked up?
vec_max <- which(is.na(colldf$interval_max_ma))
vec_min <- which(is.na(colldf$interval_min_ma))
# Use original input ages
colldf$interval_max_ma[vec_max] <- colldf$max_ma[vec_max]
colldf$interval_min_ma[vec_min] <- colldf$min_ma[vec_min]
colldf$early_stage[vec_max] <- colldf$early_interval[vec_max]
colldf$late_stage[vec_min] <- colldf$late_interval[vec_min]
# Replace max_ma and min_ma ages
colldf$max_ma <- colldf$interval_max_ma
colldf$min_ma <- colldf$interval_min_ma
# Calculate interval_mid_ma
colldf$interval_mid_ma <- (colldf$interval_max_ma + colldf$interval_min_ma) / 2
# Remove any collections with a large age range (> 50 Myr)
colldf <- colldf[-which(abs(colldf$max_ma - colldf$min_ma) > 50), ]

## Temporal binning -----------------------------------------------------
# Use the majority method
colldf <- bin_time(occdf = colldf, bins = bins, method = params$method)

# Remove data which do not hit the majority threshold (params$threshold)
colldf <- colldf[-which(colldf$overlap_percentage < params$threshold), ]

## Palaeorotate collections ---------------------------------------------
colldf <- palaeorotate(occdf = colldf,
                       lng = params$lng,
                       lat = params$lat,
                       age = params$age,
                       model = params$models,
                       method = "point",
                       uncertainty = FALSE,
                       round = NULL)
# Create equal area latitudinal bins
bins <- lat_bins_area(n = params$lat_bin)
# Save bins
saveRDS(object = bins, file = "./data/lat_bins.RDS")
# Create empty columns
colldf[, paste0(params$models, "_bin")] <- NA
# Bin collections using palaeolatitudes
# Run across models
for (i in params$models) {
  cnme <- paste0(params$p_lat, "_", i)
  # Run across bins
  for (j in seq_len(nrow(bins))) {
    vec <- which(colldf[, cnme] < bins$max[j] &
                 colldf[, cnme] > bins$min[j])
    colldf[vec, paste0(i, "_bin")] <- bins$bin[j]
  }
}
# Match up datasets -----------------------------------------------------
# Retain collections present in colldf
occdf <- occdf[which(occdf$collection_no %in% colldf$collection_no), ]
# Join datasets
m <- match(x = occdf$collection_no, table = colldf$collection_no)
# Add data
occdf[, colnames(colldf)] <- colldf[m, colnames(colldf)]
# Filter for unique occurrences from stacked collections
occdf <- distinct(occdf, lat, lng, family, genus, bin_assignment,
                  .keep_all = TRUE)
# Save processed data
saveRDS(object = occdf, file = "./data/processed/pbdb_data.RDS")
# Notify
if (params$notify) {
  beepr::beep(4)
}
