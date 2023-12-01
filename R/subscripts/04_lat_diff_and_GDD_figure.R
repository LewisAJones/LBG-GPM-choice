# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 04_lat_diff_and_GDD_figure.R
# Aim: Plot the average pairwise latitudinal difference and geodesic
#      distance through time.
# Last updated: 2023-12-01
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(grid)
library(deeptime)
library(geosphere)
library(palaeoverse)
library(dplyr)
source("./R/functions/theme_will.R")

GTS2020_periods <- time_bins(rank = "period") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
GTS2020_eras <- time_bins(rank = "era") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Load data --------------------------------------------------------------------
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf %>% select(collection_no, 
                           bin_assignment,
                           bin_midpoint, 
                           matches("GOLONKA"),
                           matches("PALEOMAP"), 
                           matches("MERDITH2021"))
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)
############################ LAT PAIRWISE DIFFERENCE ###########################
# Lat mean pairwise difference per row -----------------------------------------
colldf1 <- colldf
colldf1$mean_pair_lat_diff <- sapply(X = 1:nrow(colldf1),
                                     FUN = function(x){
                                       lat_diffGxP <- abs(colldf1$p_lat_GOLONKA[x] - colldf1$p_lat_PALEOMAP[x])
                                       lat_diffMxP <- abs(colldf1$p_lat_MERDITH2021[x] - colldf1$p_lat_PALEOMAP[x])
                                       lat_diffMxG <- abs(colldf1$p_lat_MERDITH2021[x] - colldf1$p_lat_GOLONKA[x])
                                       return(mean(c(lat_diffGxP, lat_diffMxG, lat_diffMxP)))
                                     })
colldf1 <- colldf1 %>% filter(!is.na(mean_pair_lat_diff)) #remove NAs
# Assign 5%, median and 95% quantiles ------------------------------------------
assign_quantiles <- function(time_midpoint, level){
  pairwise_lat_vec <- colldf1$mean_pair_lat_diff[which(colldf1$bin_midpoint == time_midpoint)]
  return(quantile(pairwise_lat_vec, probs = level, na.rm = TRUE))
}

med <- sapply(X = unique(colldf1$bin_midpoint),
              FUN = assign_quantiles,
              level = 0.5)
low <- sapply(X = unique(colldf1$bin_midpoint),
              FUN = assign_quantiles,
              level = 0.025)
up <- sapply(X = unique(colldf1$bin_midpoint),
             FUN = assign_quantiles,
             level = 0.975)

# Wrap up in a table and plot --------------------------------------------------
plot_df <- data.frame(bin_midpoint = unique(colldf1$bin_midpoint),
                      med = med,
                      lower = low,
                      upper = up)
p1 <- ggplot(data = plot_df, aes(x = bin_midpoint, y = med)) +
  scale_x_reverse(limits = c(542, 0),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 40, 5),
                     labels = seq(0, 40, 5)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  geom_ribbon(aes(ymin = low, ymax = upper), 
              fill = "#e7298a",
              alpha = 0.2) +
  labs(x = "Time (Ma)",
       y = "Palaeolatitudinal difference (º)") +
  theme_classic(base_size = 20) +
  theme_will() +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))

###################### MEAN PAIRWISE GEODESIC DISTANCE #########################
rm(colldf1)
# GDD function ----------------------------------------------------------
GDD <- function(x) {
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
  gd_dist <- unlist(lapply(coll_idx, FUN=GDD))
  gd_dist <- gd_dist[which(is.na(gd_dist) == FALSE)]
  #Median
  geodes <- c(geodes, quantile(gd_dist, probs = 0.5))
  #5%
  lower <- c(lower, quantile(gd_dist, probs = 0.025))
  #95%
  upper <- c(upper, quantile(gd_dist, probs = 0.975))
}
## Plot ------------------------------------------------------------------------
#Proper plot
t_bins <- readRDS("./data/time_bins.RDS")
plot_df2 <- data.frame(time = sort(t_bins$mid_ma[which(t_bins$bin %in% colldf$bin_assignment)], decreasing = TRUE),
                      GDD = geodes,
                      lower = lower,
                      upper = upper)
p2 <- ggplot(data = plot_df2, aes(x = time, y = GDD)) +
  scale_x_reverse(limits = c(542, 0),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = c(0, 1500, 3000, 4500, 6000, 7500),
                     labels = c(0, 1500, 3000, 4500, 6000, 7500)) +
  geom_point(size = 2, colour = "#e34a33") +
  geom_line(linewidth = 1, colour = "#e34a33") +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = "#e34a33",
              alpha = 0.2) +
  labs(x = "Time (Ma)",
       y = "Geodesic distance (km)") +
  theme_classic(base_size = 20) +
  theme_will() +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))

########################## Assemble the two plots ##############################
gg <- ggarrange2(p1, p2, ncol = 1, nrow = 2, labels = c("(a)", "(b)"),
                 label.args = list(gp = gpar(font = 2, cex = 2)), draw = FALSE)
# Save plot
ggsave(filename = "./figures/Lat_sd_and_GDD.png", gg,
       width = 13, height = 12)
