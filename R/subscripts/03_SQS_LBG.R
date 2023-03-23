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
# Interpolate each spatial bin ------------------------------------------
# Count unique genera in each spatio-temporal bin for each rotation model
interpolations <- data.frame()
count_NAs <- data.frame()
for (i in 1:length(params$models)){
  # Specify column name using rotation model name
  column_name <- data_sym(paste0(params$models[i], "_bin"))
  # Generate a table of sample sizes
  counts <- group_by(genus_occs, bin_assignment, !!column_name) %>%
    count()
  # Keep NA values
  keep_NAs <- filter(counts, is.na(!!column_name))
  # Rename columns
  colnames(keep_NAs) <- c("stage_bin", "paleolat_bin", "n_genera")
  # Add model label
  keep_NAs$model <- params$models[i]
  # Record NAs in separate dataframe
  count_NAs <- rbind.data.frame(count_NAs, keep_NAs[, -2])
  # Create vector of stages
  stages <- sort(unique(genus_occs$bin_assignment))
  # Filter occurrences to one stage
  for (l in 1:length(stages)) {
    one_stage <- filter(genus_occs, bin_assignment == stages[l])
    # Generate list of frequencies by stage, starting with number of 
    # samples, as needed by iNEXT
    lat_freq <- list()
    for (m in 1:12) {
      one_bin <- filter(one_stage, !!column_name == m)
      lat_list <- count(one_bin, genus) %>% arrange(desc(n)) %>%
        add_row(n = length(unique(one_bin$collection_no)),
                .before = 1) %>% select(n)
      lat_list <- unlist(lat_list, use.names = F)
      if(lat_list[1] < 3){lat_list <- NA}
      if(length(lat_list) < 4){lat_list <- NA}
      lat_freq[[m]] <- lat_list
    }
    # Name lists
    names(lat_freq) <- seq(1, 12, 1)
    # Filter out empty lists
    lat_freq <- lat_freq[!is.na(lat_freq)]
    if (length(lat_freq) > 0) {
      # Estimate D using estimateD in iNEXT
      estD <- estimateD(lat_freq, q = 0, datatype = "incidence_freq",
                        base = "coverage", level = params$quorum_level)
      # Add sample size in additional column (from first value in lists)
      estD$reference_t <- unlist(lapply(lat_freq, '[[', 1))
      # Remove values when t is more than two times the sample size
      estD[which(estD$t >= 2 * estD$reference_t),
           c("qD", "qD.LCL", "qD.UCL")] <- rep(NA, 3)
      # Add latitude bin labels
      estD$paleolat_bin <- names(lat_freq)
      # Fill in any gaps
      for (n in 1:12) {
        if ((n %in% estD$paleolat_bin) == FALSE) {
          estD <- rbind.data.frame(estD, c(rep(NA, 9), n))
        }
      }
      # Sort rows
      estD$paleolat_bin <- as.numeric(estD$paleolat_bin)
      estD <- arrange(estD, paleolat_bin)
      # Add stage label
      estD$stage <- stages[l]
      # Add model label
      estD$model <- params$models[i]
      # Add values to overall dataframe
      interpolations <- rbind.data.frame(interpolations, estD)
    }
  }
}
# Interpolate global values ---------------------------------------------
for (p in 1:length(params$models)) {
  # Specify column name using rotation model name
  column_name <- data_sym(paste0(params$models[p], "_bin"))
  # Generate list of frequencies by stage, starting with number of 
  # samples, as needed by iNEXT
  no_NA <- filter(genus_occs, !is.na(!!column_name))
  temp_freq <- list()
  for (q in 1:length(stages)){
    one_stage <- filter(no_NA, bin_assignment == stages[q])
    temp_list <- count(one_stage, genus) %>% arrange(desc(n)) %>%
                add_row(n = length(unique(one_stage$collection_no)),
                .before = 1) %>% select(n)
    temp_list <- unlist(temp_list, use.names = F)
    if(temp_list[1] < 3){temp_list <- NA}
    if(length(temp_list) < 4){temp_list <- NA}
    temp_freq[[q]] <- temp_list
  }
  # Name lists
  names(temp_freq) <- stages
  # Estimate D using estimateD in iNEXT
  estD <- estimateD(temp_freq, q = 0, datatype = "incidence_freq",
                    base = "coverage", level = params$quorum_level)
  # Add sample size in additional column (from first value in lists)
  estD$reference_t <- unlist(lapply(temp_freq, '[[', 1))
  # Remove values when t is more than two times the sample size
  estD[which(estD$t >= 2 * estD$reference_t),
       c("qD", "qD.LCL", "qD.UCL")] <- rep(NA, 3)
  # Add latitude bin NA label
  estD$paleolat_bin <- NA
  # Add stage label
  estD$stage <- stages
  # Add model label
  estD$model <- params$models[p]
  # Add values to overall dataframe
  interpolations <- rbind.data.frame(interpolations, estD)
}
# Tidy dataframe
interpolations <- select(interpolations, stage, paleolat_bin, model,
                         reference_t, t, Method, SC, qD, qD.LCL, qD.UCL)
interpolations <- arrange(interpolations, model, stage, paleolat_bin)
# Save interpolated estimates and NA counts
saveRDS(object = interpolations,
        file = "./data/processed/interpolations.RDS")
saveRDS(object = count_NAs, file = "./data/processed/NA_counts.RDS")
# Notify
if (params$notify) {
  beepr::beep(4)
}