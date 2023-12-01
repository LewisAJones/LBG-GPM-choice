# Header -----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_nb_collections_figure.R
# Aim: Plot number of collections throug time and binned in the different 
#      palaeolat bins per stages (Fig.S5).
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

## Load libraries -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(deeptime)
## 1. Overall plot ------------------------------------------------------
#Load data
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf %>% select(collection_no, bin_assignment)
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)
#Evaluate number of collection through time
n_col <- sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                FUN = function(x){
                  idx <- which(colldf$bin_assignment == x)
                  return(length(idx))
                  })

nb_coll.df <- data.frame(time_bin = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                         number_of_collections = n_col)
nb_coll.df$mid_time <- sapply(X = nb_coll.df$time_bin,
                              FUN = function(x){
                                return(time_bins$mid_ma[which(time_bins$bin == x)])
                                })

col_plot <- ggplot(data = nb_coll.df, aes(x = mid_time, y = number_of_collections)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 4100),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, 2000, 3000, 4000)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  labs(x = "Time (Ma)",
       y = "Number of collections") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
ggsave("./figures/Number_of_collections.png", col_plot, width = 13, height = 6)

## 2. Per-stage plot -------------------------------------------------
#Load data
occdf <- readRDS("./data/processed/pbdb_data.RDS")
#Retain features of interest and filter by collection
tmp <- occdf %>% select(collection_no, bin_assignment, PALEOMAP_bin, GOLONKA_bin, MERDITH2021_bin)
tmp <- unique(tmp)
#Create collection dataframe
colldf <- rbind(tmp[, 1:2], tmp[, 1:2], tmp[, 1:2])
colldf$model <- c(rep("GOLONKA", nrow(tmp)),
                  rep("PALEOMAP", nrow(tmp)),
                  rep("MERDITH2021", nrow(tmp)))
colldf$plat_bin <- c(tmp$GOLONKA_bin, 
                     tmp$PALEOMAP_bin, 
                     tmp$MERDITH2021_bin)
#Assess number of collections per lat bin through time
colldf1 <- colldf %>% 
  group_by(model, bin_assignment, plat_bin) %>% 
  count(plat_bin) %>% 
  rename("occ_number" = n) %>% 
  ungroup() %>%
  complete(bin_assignment = time_bins$bin, plat_bin, model) %>%
  full_join(time_bins, by = c("bin_assignment" = "bin")) %>% 
  full_join(lat_bins, by = c("plat_bin" = "bin")) %>%
  filter(!is.na(model)) %>%
  filter(!is.na(plat_bin))
#Plot
p <- ggplot(data = colldf1, aes(x = mid,
                                y = occ_number,
                                colour = model)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  facet_wrap(~factor(interval_name, levels = rev(time_bins$interval_name)), ncol = 9, scales = "free") +
  labs(y = "Number of collections",
       x = "Palaeolatitudinal bin") +
  theme_bw(base_size = 18) +
  theme(strip.text.x = element_text(size = 14),
        legend.position = "top")

#Update strip colours
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
for (i in strip_t) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <-
    time_bins$colour[match(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label,
                           time_bins$interval_name)]
}
#Save
ggsave("./figures/Nb_collections_per_stage.png", g, width = 18, height = 22)
