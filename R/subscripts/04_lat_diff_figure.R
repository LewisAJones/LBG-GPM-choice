# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 04_lat_diff_figure.R
# Aim: Plot the average pairwise latitudinal difference through time.
# Last updated: 2026-09-02
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(deeptime)
library(palaeoverse)
library(dplyr)
source("R/functions/theme_will.R")

ics_periods <- time_bins(scale = "international periods") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
ics_eras <- time_bins(scale = "international eras") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Load data --------------------------------------------------------------------
t_bins <- read.csv("data/time_bins.csv") #time bins
occdf <- readRDS("data/processed/pbdb_data.RDS") #occurrence data
colldf <- occdf %>% select(collection_no, 
                           bin_assignment,
                           bin_midpoint, 
                           matches("GOLONKA"),
                           matches("PALEOMAP"), 
                           matches("MERDITH2021"),
                           matches("TorsvikCocks2017"))
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)

################### PAIRWISE LATITUDINAL DIFFERENCE (degrees) ##################
# Create median pairwise lat. diff. for each occurrence --------------------------
colldf$median_pair_lat_diff <- sapply(X = 1:nrow(colldf),
                                      FUN = function(x){
                                        lat_diffGxT <- abs(colldf$p_lat_GOLONKA[x] - colldf$p_lat_TorsvikCocks2017[x])
                                        lat_diffGxP <- abs(colldf$p_lat_GOLONKA[x] - colldf$p_lat_PALEOMAP[x])
                                        lat_diffPxT <- abs(colldf$p_lat_PALEOMAP[x] - colldf$p_lat_TorsvikCocks2017[x])
                                        lat_diffMxT <- abs(colldf$p_lat_MERDITH2021[x] - colldf$p_lat_TorsvikCocks2017[x])
                                        lat_diffMxP <- abs(colldf$p_lat_MERDITH2021[x] - colldf$p_lat_PALEOMAP[x])
                                        lat_diffMxG <- abs(colldf$p_lat_MERDITH2021[x] - colldf$p_lat_GOLONKA[x])
                                        return(median(c(lat_diffGxP, lat_diffMxG, lat_diffMxP), na.rm = TRUE))
                                      })
colldf <- colldf %>% filter(!is.na(median_pair_lat_diff)) #remove NAs
write.csv(colldf, "results/coll_plat_diff.csv", row.names = FALSE)
# Assign 5%, median and 95% quantiles ------------------------------------------
assign_quantiles <- function(time_midpoint, level){
  pairwise_lat_vec <- colldf$median_pair_lat_diff[which(colldf$bin_midpoint == time_midpoint)]
  return(quantile(pairwise_lat_vec, probs = level, na.rm = TRUE))
}
med <- sapply(X = unique(colldf$bin_midpoint),
              FUN = assign_quantiles,
              level = 0.5)
low <- sapply(X = unique(colldf$bin_midpoint),
              FUN = assign_quantiles,
              level = 0.05)
up <- sapply(X = unique(colldf$bin_midpoint),
             FUN = assign_quantiles,
             level = 0.95)
# Wrap up in a table and plot --------------------------------------------------
plot_df <- data.frame(time = unique(colldf$bin_midpoint),
                      med = med,
                      lower = low,
                      upper = up,
                      metric = rep("Palaeolatitudinal Difference (°)", length(med)))
# Order by time
plot_df <- plot_df[order(plot_df$time), ]
# Save datafile
write.csv(plot_df, "results/plat_diff_degrees.csv", row.names = FALSE)

# Plot the two -----------------------------------------------------------------
p <- ggplot(data = plot_df, aes(x = time, y = med)) +
  scale_x_reverse(limits = c(538.800, 0),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, NA)) +
  geom_point(size = 2, aes(color = metric)) +
  geom_line(linewidth = 1, aes(color = metric)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = metric),
              alpha = 0.2) +
  scale_fill_manual(values = c("#e7298a")) +
  scale_color_manual(values = c("#e7298a")) +
  labs(x = "Time (Ma)",
       y = "Palaeolatitudinal Difference (°)",
       fill = NULL) +
  coord_geo(list("bottom", "bottom"), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE),
            size = 6) +
  theme_classic(base_size = 24) +
  theme_will(axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 18),
             axis.text.y = element_text(size = 18),
             legend.position = "none")
p
# Save plot
ggsave(filename = "figures/lat_diff.png", p,
       width = 13, height = 7)
