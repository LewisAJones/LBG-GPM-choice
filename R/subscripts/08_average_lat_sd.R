# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 08_average_lat_sd.R
# Last updated: 2023-03-07
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Load libraries --------------------------------------------------------
library(ggplot2)
library(deeptime)
library(dplyr)
source("./R/functions/theme_will.R")
# Load data -------------------------------------------------------------
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf[, c("collection_no", "bin_midpoint",
                    "p_lat_GOLONKA", "p_lat_MERDITH2021", "p_lat_PALEOMAP")]
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)
# Calculate palaeolatitudinal sd ----------------------------------------
colldf$sd <- apply(colldf[, c("p_lat_GOLONKA", "p_lat_MERDITH2021", "p_lat_PALEOMAP")], 1, sd) 
# Summary per bin -------------------------------------------------------
summ <- colldf %>%
  group_by(bin_midpoint) %>%
  summarise(mean = mean(sd, na.rm = TRUE))
# Plot data -------------------------------------------------------------
p <- ggplot(data = summ, aes(x = bin_midpoint, y = mean)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, 5),
                     labels = seq(0, 20, 5)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  labs(x = "Time (Ma)",
       y = "Mean palaeolatitudinal standard deviation (º)") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             aspect.ratio = 0.5) +
  coord_geo(lwd = 1)
ggsave("./figures/mean_sd.png", p, width = 10, height = 6)
