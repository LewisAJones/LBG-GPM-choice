###############################################################################################################
############### Calculate average geodesic distance between palaeorotations of each collections ###############
###############################################################################################################
## Load libraries --------------------------------------------------------------
library(dplyr)
library(geosphere)
library(ggplot2)
library(deeptime)
## Read palaeorotated cleaned occurrence dataframe and extract collections -----
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf %>% select(collection_no, bin_assignment, matches("GOLONKA"), matches("PALEOMAP"), matches("MERDITH2021"))
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)
## Assess mean pairwise geaodesic distance per lat bin for each time bin -------
geodes <- c()
for(t in unique(sort(colldf$bin_assignment, decreasing = FALSE))){
  coll_idx <- which(colldf$bin_assignment == t)
  gd_dist <- lapply(coll_idx, function(x) {
    gol <- as.numeric(colldf[x, c("p_lng_GOLONKA", "p_lat_GOLONKA")])
    pal <- as.numeric(colldf[x, c("p_lng_PALEOMAP", "p_lat_PALEOMAP")])
    mer <- as.numeric(colldf[x, c("p_lng_MERDITH2021", "p_lat_MERDITH2021")])
    tmp <-  matrix(c(gol, pal, mer),
                   ncol = 2,
                   byrow = TRUE)
    dist_mat <- distm(tmp)
    dist_mat <- as.numeric(dist_mat[upper.tri(dist_mat, diag = FALSE)])
    return(mean(dist_mat, na.rm = TRUE) / (10^3))
  })
  gd_dist <- do.call(rbind, gd_dist)
  geodes <- c(geodes, mean(gd_dist, na.rm = TRUE))
}
## Plot ------------------------------------------------------------------------
  #Harmonise plotting setup with the rest of the plots
theme_will <- function(...) {
  theme(axis.ticks = element_line(color = "black", linewidth = 1),
        axis.line = element_blank(),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(colour = "black", size = 12),
        panel.border = element_rect(linetype = "solid", colour = "black",
                                    fill = NA, linewidth = 2),
        panel.background = element_blank(),
        ...)
}
  #Proper plot
t_bins <- readRDS("./data/time_bins.RDS")
plot_df <- data.frame(time = sort(t_bins$mid_ma[which(t_bins$bin %in% colldf$bin_assignment)], decreasing = TRUE),
                      GDD = geodes)
p <- ggplot(data = plot_df, aes(x = time, y = GDD)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(600, 3400),
                     breaks = c(1000, 1500, 2000, 2500, 3000),
                  labels = c(1000, 1500, 2000, 2500, 3000)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  labs(x = "Time (Ma)",
       y = "Average geodesic distance (km)") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(lwd = 1)
ggsave("./figures/Average_GDD.pdf", p, width = 13, height = 6)
