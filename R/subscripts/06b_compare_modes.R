# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06b_compare_modes.R
# Description: assess the number of modes in the diversity curves obtained with the 3 GPMs
# Last updated: 2023-03-07
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Load libraries --------------------------------------------------------
library(multimode)
library(stats)
library(dplyr)
library(ggplot2)
library(deeptime)
source("./R/options.R")
time_bins <- readRDS("./data/time_bins.RDS")
## Load data -------------------------------------------------------------------
occdf <- readRDS("./data/processed/pbdb_data.RDS")[, c("genus", "max_ma", "min_ma", "bin_assignment", "p_lat_PALEOMAP", "p_lat_GOLONKA", "p_lat_MERDITH2021")]
#rename time bin column to avoid ambiguities with lat bin one
occdf <- occdf %>% rename(time_bin_assignment = "bin_assignment")
#add lat bin columns (Bin collections using palaeolatitudes) 
lat_bins <- readRDS("./data/lat_bins.RDS")
occdf[, paste0(params$models, "_bin")] <- NA
for (i in params$models) {
  cnme <- paste0(params$p_lat, "_", i)
  for (j in seq_len(nrow(lat_bins))) {
    vec <- which(occdf[, cnme] < lat_bins$max[j] &
                   occdf[, cnme] > lat_bins$min[j])
    occdf[vec, paste0(i, "_bin")] <- lat_bins$bin[j]
  }
}

## Seek number of modes in binned latitudinal diversity ------------------------
  #returns a random number comprised between (x-1) and x
assign_rand <- function(x){
  rand <- runif(n = 1, min = x-1, max = x)
  return(rand)
}
  #assess the number of modes with N independent repeats of adding the random compound
mode_gathering <- function(subset, N=10){
  MODES <- c()
  for(i in 1:N){
    #assign each occurrence to a random number between (k-1) and k, where k is its lat bin
    subset$rand <- unlist(lapply(X = subset[, paste0(mdl, "_bin")], FUN = assign_rand))
    #evaluate band width
    bwRT <- stats::bw.nrd0(subset$rand)
    #evaluate number of modes
    nmodes <- multimode::nmodes(data = subset$rand, bw = bwRT)
    MODES <- c(MODES, nmodes)
  }
  return(round(mean(MODES, na.rm = TRUE)))
}
  #Run across models
mode_df <- data.frame(time_bin = sort(unique(occdf$time_bin_assignment)))
for(mdl in params$models){
  print(mdl)
  mode_df[, paste0("Nmodes_", mdl)] <- NA
  for(t in mode_df$time_bin){
    subset <- occdf[which( (occdf$time_bin_assignment == t) & (is.na(occdf[,paste0(mdl, "_bin")]) == FALSE)),
                    c("genus", paste0(mdl, "_bin"))]
    if(nrow(subset) == 0){
      mode_df[which(mode_df$time_bin == t), paste0("Nmodes_", mdl)] <- NA
    }
    else{
      # Remove duplicates (same genera present several times within the same time bin)
      subset <- unique(subset)
      row.names(subset) <- 1:nrow(subset)
      mode_df[which(mode_df$time_bin == t), paste0("Nmodes_", mdl)] <- mode_gathering(subset = subset, N=50)
    }
  }
}

## Save ------------------------------------------------------------------------
saveRDS(object = mode_df,
        file = "./data/mode_counts.RDS")

## Informal plot ---------------------------------------------------------------
mode_df$Nmodes_PALEOMAP[3] <- NA
mode_df$Nmodes_MERDITH2021[3] <- NA
#time equivalent of mode_df bins
find_equiv <- function(bin){
  return(time_bins$mid_ma[which(time_bins$bin == bin)])
}

plot_df <- data.frame(Time = rep(unlist(lapply(X = mode_df$time_bin, FUN = find_equiv))),
                      N_modes = c(mode_df$Nmodes_PALEOMAP, mode_df$Nmodes_GOLONKA, mode_df$Nmodes_MERDITH2021),
                      model = c(rep("PALEOMAP", nrow(mode_df)), rep("GOLONKA", nrow(mode_df)), rep("MERDITH2021", nrow(mode_df))))

p <- ggplot(data = plot_df, aes(x = Time, y = N_modes, group = model, colour = model)) +
  geom_step(linewidth = 1) +
  scale_x_reverse() +
  scale_y_continuous(limits = c(0.75,3.25), breaks = 1:3) +
  scale_color_manual(values = c("#fe9929", "#dd3497", "#3690c0")) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
  deeptime::coord_geo(dat = "periods")
  #Save
ggsave(plot = p,
       filename = "./figures/mode_number_draft.png",
       height = 200,
       width = 300,
       units = "mm",
       dpi = 600)

## For visualisation, subset PALEOMAP in time bin 79 (Danian) ------------------
subset <- occdf[which( (occdf$time_bin_assignment == 79) & (is.na(occdf$PALEOMAP_bin) == FALSE)),
                c("genus", "PALEOMAP_bin", "p_lat_PALEOMAP")]
  # Remove duplicates (same genera present several times within the same time bin)
row.names(subset) <- 1:nrow(subset)
idx <- row.names(unique(subset[, c("genus", "PALEOMAP_bin")]))
idx <- unlist(lapply(X = idx, FUN = as.numeric))
subset <- subset[idx,]
  # Add the random number and plot
subset$rand <- unlist(lapply(X = 1:nrow(subset),
                             FUN = assign_rand))
par(mfrow = c(1,2))
plot(density(subset$p_lat_PALEOMAP))
plot(density(subset$rand))

