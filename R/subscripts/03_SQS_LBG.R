# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 03_SQS_LBG.R
# Last updated: 2023-03-13
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Load libraries --------------------------------------------------------
library(palaeoverse)
library(dplyr)
library(rlang)
library(iNEXT)
# Load project options
source("./R/options.R")
# Get data --------------------------------------------------------------
occdf <- readRDS("./data/processed/pbdb_data.RDS")
# Process data ----------------------------------------------------------
# Collapse subgenera (remove characters from space onwards)
if (params$collapse_subgenera = TRUE) {
  occdf$genus <- sub(" .*", "", occdf$genus)
}
# Filter for unique occurrences from each collection
occdf <- distinct(occdf, collection_no, family, genus, bin_assignment,
                  .keep_all = TRUE)
# Filter for unique occurrences from stacked collections
# (those with same lat/lng, usually beds from the same section)
occdf$lng <- round(occdf$lng, digits = params$n_decs)
occdf$lat <- round(occdf$lat, digits = params$n_decs)
occdf <- distinct(occdf, lat, lng, family, genus, bin_assignment,
                  .keep_all = TRUE)
# Streamline to necessary variables
genus_occs <- occdf[, c("collection_no", "phylum", "family", "genus",
                        "bin_assignment")]
genus_occs <- cbind.data.frame(genus_occs, occdf[, c(
  grep("p_lng_", x = colnames(occdf)),
  grep("p_lat_", x = colnames(occdf)),
  grep("_bin", x = colnames(occdf)))])
# Create abundance strings ----------------------------------------------
# Count unique genera in each spatio-temporal bin for each rotation model
genus_counts <- data.frame()
count_NAs <- data.frame()
for (i in 1:length(params$models)){
  # Specify column name using rotation model name
  column_name <- data_sym(paste0(params$models[i], "_bin"))
  # Generate a table of sample sizes
  counts <- group_by(genus_occs, bin_assignment, !!column_name) %>%
    count()
  # Remove NA values into a separate data frame
  keep_NAs <- filter(counts, is.na(!!column_name))
  counts <- filter(counts, !is.na(!!column_name))
  # Rename columns
  colnames(counts) <- c("stage_bin", "paleolat_bin", "n_genera")
  colnames(keep_NAs) <- c("stage_bin", "paleolat_bin", "n_genera")
  # Add model label
  keep_NAs$model <- params$models[i]
  # Record NAs in separate dataframe
  count_NAs <- rbind.data.frame(count_NAs, keep_NAs[, -2])
}
# Save interpolated estimates and NA counts
# saveRDS(object = genus_counts, file = "./data/processed/genus_counts.RDS")
saveRDS(object = count_NAs, file = "./data/processed/NA_counts.RDS")