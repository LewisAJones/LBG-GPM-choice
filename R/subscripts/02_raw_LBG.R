# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 02_raw_LBG.R
# Aim: Estimate diversity through time using each GPM with raw counts.
# Last updated: 2023-03-13
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Load libraries --------------------------------------------------------
library(palaeoverse)
library(dplyr)
library(rlang)
# Load project options
source("./R/options.R")
# Get data --------------------------------------------------------------
occdf <- readRDS("./data/processed/pbdb_data.RDS")
# Process data ----------------------------------------------------------
# Collapse subgenera (remove characters from space onwards)
if (params$collapse_subgenera == TRUE) {
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
# Count genera ----------------------------------------------------------
# Count unique genera in each spatio-temporal bin for each rotation model
genus_counts <- data.frame()
for (i in 1:length(params$models)) {
  # Specify column name using rotation model name
  column_name <- data_sym(paste0(params$models[i], "_bin"))
  # Filter genera in each spatio-temporal bin, so each row is a unique
  # genus, according to the relevant plate model spatial bins
  one_model <- distinct(genus_occs, family, genus, bin_assignment,
                         !!column_name, .keep_all = TRUE)
  # Generate genus counts per time bin per palaeolatitude bin
  counts <- group_by(one_model, bin_assignment, !!column_name) %>%
    count()
  # Remove NA values
  counts <- filter(counts, !is.na(!!column_name))
  # Rename columns
  colnames(counts) <- c("stage_bin", "paleolat_bin", "n_genera")
  # Create vector of stages
  stages <- sort(unique(genus_occs$bin_assignment))
  # Fill in any gaps
  counts <- as.data.frame(counts)
  for (j in 1:length(stages)) {
    one_stage <- filter(counts, stage_bin == stages[j])
    for (k in 1:params$lat_bin) {
      if ((k %in% one_stage$paleolat_bin) == FALSE) {
        counts <- rbind(counts, c(stages[j], k, 0))
      }
    }
  }
  # Filter genera in each temporal bin, so each row is a unique
  # genus, according to the relevant plate model spatial bins
  one_model <- distinct(one_model, family, genus, bin_assignment,
                        .keep_all = TRUE)
  # Generate genus counts globally per time bin, removing NAs
  global_counts <- filter(one_model, !is.na(!!column_name)) %>%
    group_by(bin_assignment) %>% count()
  # Add 'global' label
  global_counts$paleolat_bin <- NA
  global_counts <- select(global_counts, bin_assignment, paleolat_bin, n)
  colnames(global_counts) <- c("stage_bin", "paleolat_bin", "n_genera")
  # Append to other table
  counts <- rbind.data.frame(counts, global_counts)
  counts <- arrange(counts, stage_bin, paleolat_bin)
  # Add model label
  counts$model <- params$models[i]
  # Add counts to overall dataframe
  genus_counts <- rbind.data.frame(genus_counts, counts)
}
# Save genus counts
saveRDS(object = genus_counts, file = "./data/processed/genus_counts.RDS")