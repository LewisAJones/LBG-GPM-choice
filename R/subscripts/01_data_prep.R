# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 01_data_prep.R
# Last updated: 2023-03-06
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Load libraries --------------------------------------------------------
library(palaeoverse)
# Load project options
source("./R/options.R")
# Load custom functions
source("./R/functions/equal_area_lat_bins.R")
# Get data --------------------------------------------------------------
if (params$download) {
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
# Bin data by collection for faster processing
colldf <- occdf[, c("collection_no", "lng", "lat", "min_ma", "max_ma")]
colldf <- unique(colldf)
# Create time bins
bins <- time_bins(interval = "Phanerozoic",
                  rank = params$rank,
                  scale = params$GTS)
# Save time bins
saveRDS(object = bins, file = "./data/time_bins.RDS")
# Temporal binning (using the majority method)
colldf <- bin_time(occdf = colldf,
                   bins = bins,
                   method = params$method)
# Remove data which do not hit the majority threshold (params$threshold)
colldf <- colldf[-which(colldf$overlap_percentage < params$threshold), ]
# Palaeorotate collections (this can take some time, put on a coffee)
colldf <- palaeorotate(occdf = colldf,
                       lng = params$lng,
                       lat = params$lat,
                       age = params$age,
                       model = params$models,
                       method = "point",
                       uncertainty = FALSE,
                       round = NULL)
# Create equal area latitudinal bins
bins <- equal_area_lat_bins(n = params$lat_bin, by = 0.01)
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
# Join datasets (for loop used to prevent vector memory exhaustion)
m <- match(x = occdf$collection_no, table = colldf$collection_no)
# Which columns are the palaeocoords?
coord_cols <- sort(c(grep("p_lng_", x = colnames(colldf)),
                     grep("p_lat_", x = colnames(colldf))))
# Which columns are not already present in occdf?
cols <- !colnames(colldf) %in% colnames(occdf)
# Bind data
occdf <- cbind.data.frame(occdf, colldf[m, cols])
# Save processed data
saveRDS(object = occdf, file = "./data/processed/pbdb_data.RDS")
# Notify
if (params$notify) {
  beepr::beep(4)
}