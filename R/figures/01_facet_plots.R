# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 01_facet_plots.R
# Last updated: 2023-03-23
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Libraries -------------------------------------------------------------
library(grid)
library(ggplot2)
# Data ------------------------------------------------------------------
# Load bins
lats <- readRDS("./data/lat_bins.RDS")
stages <- readRDS("./data/time_bins.RDS")
# Load data
df <- readRDS("./data/processed/genus_counts.RDS")
# Wrangling -------------------------------------------------------------
# Set factor levels
stages$interval_name <- factor(x = stages$interval_name, levels = stages$interval_name)
# Rename columns for merging
colnames(stages)[which(colnames(stages) == "bin")] <- "stage_bin"
colnames(lats)[which(colnames(lats) == "bin")] <- "paleolat_bin"
# Join datasets
df <- merge(x = df, y = stages, by = c("stage_bin"), all = TRUE)
df <- merge(x = df, y = lats, by = c("paleolat_bin"), all = TRUE)
# Add unique column
df$unique <- paste0(df$stage_bin, df$model)
# Bins present in dataset
bins <- unique(df$unique)
# Add column for populating
df$normalised_richness <- NA
for (i in bins) {
  vec <- which(df$unique == i)
  df$normalised_richness[vec] <- df$n_genera[vec] / max(df$n_genera[vec])
}
# Plot data -------------------------------------------------------------
# number of rows desired
n <- 10
p <- ggplot(data = df, aes(x = mid,
                           y = normalised_richness,
                           colour = model,
                           strip)) +
            geom_line(linewidth = 0.75, alpha = 1) +
            scale_colour_viridis_d(option = "plasma") +
            facet_wrap(~interval_name, nrow = n) +
            labs(y = "Normalised genus richness (raw)",
                 x = "Palaeolatitude (\u00B0)") +
            theme_bw() +
            theme(aspect.ratio = 1,
                  legend.position = "top") 
#p
# Update strip colours (work in progress doesn't seem to work well!)
# g <- ggplot_gtable(ggplot_build(p))
# strip_t <- which(grepl('strip-t', g$layout$name))
# fills <- c("white", "white", "white", stages$colour)
# k <- 1
# for (i in strip_t) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   if (length(j) == 0) {
#     next
#   }
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k + 1
# }

#grid.draw(g)

# Save plot -------------------------------------------------------------
ggsave("./figures/raw_LBGs.png", plot = p,
       units = "mm", width = 250, height = 250, dpi = 300, scale = 3)
  

