# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 07_lat_sd_and_mean_pairwise_GDD.R
# Last updated: 2023-03-07
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(deeptime)
library(geosphere)
library(dplyr)
library(Rmisc)
source("./R/functions/theme_will.R")

# Load data -------------------------------------------------------------
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf %>% select(collection_no, bin_assignment, bin_midpoint, matches("GOLONKA"), matches("PALEOMAP"), matches("MERDITH2021"))
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)
##### LAT SD #####
# Calculate palaeolatitudinal sd ----------------------------------------
colldf1 <- colldf
colldf1$sd <- apply(X = colldf1[, c("p_lat_GOLONKA", "p_lat_MERDITH2021", "p_lat_PALEOMAP")],
                   MARGIN = 1,
                   FUN = sd) 
# Function to jointly assess mean, 5% and 95% ---------------------------
medCi <- function(midpoint, Which=c("mean", "5%", "95%")){
  stdev <- colldf1$sd[which(colldf$bin_midpoint == midpoint)]
  stdev <- stdev[-which(is.na(stdev))]
  ci <- CI(stdev)
  if(Which == "mean"){
    return(ci[2])
  }
  else if(Which == "5%"){
    return(ci[3])
  }
  else if(Which == "95%"){
    return(ci[1])
  }
}
# Summary per bin -------------------------------------------------------
summ <- data.frame(bin_midpoint = unique(colldf1$bin_midpoint))
summ$mean <- unlist(lapply(X = summ$bin_midpoint, FUN = medCi, Which = "mean"))
summ$lower <- unlist(lapply(X = summ$bin_midpoint, FUN = medCi, Which = "5%"))
summ$upper <- unlist(lapply(X = summ$bin_midpoint, FUN = medCi, Which = "95%"))
summ <- summ %>% filter(is.na(mean) == FALSE)
# Plot data -------------------------------------------------------------
p1 <- ggplot(data = summ, aes(x = bin_midpoint, y = mean)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, 5),
                     labels = seq(0, 20, 5)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = "#e7298a",
              alpha = 0.2) +
  labs(x = "Time (Ma)",
       y = "Mean palaeolatitudinal standard deviation (º)") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             aspect.ratio = 0.5) +
  coord_geo(lwd = 1)

#### MEAN PAIRWISE GEODESIC DISTANCE ####
rm(colldf1)
# GDD function ----------------------------------------------------------
GDD <- function(x, Which=c("mean", "5%", "95%")) {
  gol <- as.numeric(colldf[x, c("p_lng_GOLONKA", "p_lat_GOLONKA")])
  pal <- as.numeric(colldf[x, c("p_lng_PALEOMAP", "p_lat_PALEOMAP")])
  mer <- as.numeric(colldf[x, c("p_lng_MERDITH2021", "p_lat_MERDITH2021")])
  tmp <-  matrix(c(gol, pal, mer),
                 ncol = 2,
                 byrow = TRUE)
  dist_mat <- distm(tmp)
  dist_mat <- as.numeric(dist_mat[upper.tri(dist_mat, diag = FALSE)])
  dist_mat <- as.matrix(dist_mat)
  return(mean(dist_mat, na.rm = TRUE) / 10**3)
}
# Assess mean pairwise geodesic distance per lat bin for each time bin --
geodes <- c()
lower <- c()
upper <- c()
for(t in unique(sort(colldf$bin_assignment, decreasing = FALSE))){
  coll_idx <- which(colldf$bin_assignment == t)
  #Vector of the mean pairwise geodesic distances per occurrence 
  gd_dist <- unlist(lapply(coll_idx, FUN=GDD, Which="mean"))
  gd_dist <- gd_dist[which(is.na(gd_dist) == FALSE)] #otherwise hampers CI function
  ci <- CI(gd_dist)
  #Mean
  geodes <- c(geodes, ci[2])
  #5%
  lower <- c(lower, ci[3])
  #95%
  upper <- c(upper, ci[1])
}
## Plot ------------------------------------------------------------------------
#Proper plot
t_bins <- readRDS("./data/time_bins.RDS")
plot_df <- data.frame(time = sort(t_bins$mid_ma[which(t_bins$bin %in% colldf$bin_assignment)], decreasing = TRUE),
                      GDD = geodes,
                      lower = lower,
                      upper = upper)
p2 <- ggplot(data = plot_df, aes(x = time, y = GDD)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 7600),
                     breaks = c(0, 1500, 3000, 4500, 6000, 7500),
                     labels = c(0, 1500, 3000, 4500, 6000, 7500)) +
  geom_point(size = 2, colour = "#e34a33") +
  geom_line(linewidth = 1, colour = "#e34a33") +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = "#e34a33",
              alpha = 0.2) +
  labs(x = "Time (Ma)",
       y = "Average geodesic distance (km)") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(lwd = 1)


#### Assemble the two plots ####
p <- ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("(a)", "(b)"))
# Save plot
ggsave(filename = "./figures/Lat_sd_and_GDD.png",
       height = 300,
       width = 300,
       units = "mm",
       dpi = 600)
