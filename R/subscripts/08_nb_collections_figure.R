# Header -----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_nb_collections_figure.R
# Aim: Plot number of collections through time and binned in the different 
#      palaeolat bins per stage (Fig.S5).
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries -------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(palaeoverse)
library(deeptime)
# 1. Overall plot ------------------------------------------------------
#Plotting accessories
source("./R/functions/theme_will.R")
GTS2020_periods <- time_bins(rank = "period") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
GTS2020_eras <- time_bins(rank = "era") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
#Load data
time_bins <- readRDS("./data/time_bins.RDS")
lat_bins <- readRDS("./data/lat_bins.RDS")
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf %>% 
  select(collection_no, bin_assignment, matches("p_lat"), bin_midpoint) %>% 
  group_by(collection_no, bin_assignment, p_lat_PALEOMAP, p_lat_GOLONKA, p_lat_MERDITH2021, bin_midpoint) %>% 
  #subset based on unique values of collection number
  distinct(collection_no)

#Total number of collection through time
n_col <- sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                FUN = function(x){
                  idx <- which(colldf$bin_assignment == x)
                  return(length(idx))
                  })

nb_coll.df <- data.frame(mid_time = unique(sort(colldf$bin_midpoint, decreasing = FALSE)), #increasing order as the oldest bin has the smallest bin assignment number
                         number_of_collections = n_col)

col_plot <- ggplot(data = nb_coll.df, aes(x = mid_time, y = number_of_collections)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 4100),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, 2000, 3000, 4000)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  ggtitle("Stage-level bins") +
  labs(x = "Time (Ma)",
       y = "Number of collections") +
  theme_will(plot.title = element_text(size = 18, hjust = .5),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))

# 2. Total number of collection available for each model ------------
avail_mdl <- function(bin, model, ds=colldf){
  idx <- which(ds$bin_assignment == bin)
  nas <- which(is.na(ds[idx, paste0("p_lat_", model)]))
  return((length(idx)-length(nas)))
}
coll_av.df <- data.frame(mid_time = rep(unique(sort(colldf$bin_midpoint, decreasing = FALSE)), 3),
                         model = c(rep("GOLONKA", length(unique(colldf$bin_assignment))),
                                   rep("MERDITH2021", length(unique(colldf$bin_assignment))),
                                   rep("PALEOMAP", length(unique(colldf$bin_assignment)))),
                         coll_av = c( sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                                             FUN = avail_mdl,
                                             model = "GOLONKA"),
                                      sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                                             FUN = avail_mdl,
                                             model = "MERDITH2021"),
                                      sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                                             FUN = avail_mdl,
                                             model = "PALEOMAP")))
col_av_plot <- ggplot(data = coll_av.df, aes(x = mid_time, y = coll_av, colour = model)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 4100),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, 2000, 3000, 4000)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  ggtitle("Stage-level bins") +
  labs(x = "Time (Ma)",
       y = "Number of available collections") +
  theme_will(plot.title = element_text(size = 18, hjust = .5),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.position = "top") +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))

# 3. Collections/myr through time -----------------------------------
rm(colldf)
colldf <- occdf %>% 
  select(collection_no, matches("p_lat"), min_ma, max_ma) %>% 
  group_by(collection_no, p_lat_PALEOMAP, p_lat_GOLONKA, p_lat_MERDITH2021, min_ma, max_ma) %>% 
  distinct(collection_no)
#create time bins with a 1Myr step
tb1 <- data.frame(bin = 1:541,
                  max_ma = seq(541, 1, -1),
                  mid_ma = seq(540.5, .5, -1),
                  min_ma = seq(540, 0, -1))
#bin fossil collections based 
cdf <- bin_time(occdf = colldf,
                bins = tb1,
                method = "mid")
#prepare plotting dataset
#for some reason, the dplyr::count function decided not to work anymore
count <- function(value){
  return(length(which(cdf$bin_midpoint == value)))
}

plot_df <- data.frame(bin_midpoint = unique(sort(cdf$bin_midpoint)),
                      n = sapply(X = unique(sort(cdf$bin_midpoint)),
                                 FUN = count))
plot_df$bin_assignment <- sapply(X = plot_df$bin_midpoint,
                                 FUN = function(x){
                                   return(tb1$bin[which(tb1$mid_ma == x)])
                                 })
## plot => all collections
plot_all <- ggplot(data = plot_df, aes(x = bin_midpoint, y = n)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, 2000, 3000, 4000)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  ggtitle("1Myr-long bins") +
  labs(x = "Time (Ma)",
       y = "Number of collections (1My bins)") +
  theme_will(plot.title = element_text(size = 18, hjust = .5),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
#save with total stage-level bins
p <- ggarrange(col_plot, plot_all, nrow = 2, 
               labels = c("(a)", "(b)"), font.label = list(size = 12))
ggsave("./figures/Number_of_collections_total.png", p, width = 13, height = 12)

## plot => available collections per GPM
coll_av.df1 <- data.frame(mid_time = rep(unique(sort(plot_df$bin_midpoint, decreasing = FALSE)), 3),
                          model = c(rep("GOLONKA", length(unique(plot_df$bin_assignment))),
                                   rep("MERDITH2021", length(unique(plot_df$bin_assignment))),
                                   rep("PALEOMAP", length(unique(plot_df$bin_assignment)))),
                          coll_av = c( sapply(X = unique(sort(plot_df$bin_assignment, decreasing = TRUE)),
                                             FUN = avail_mdl,
                                             model = "GOLONKA",
                                             ds = cdf),
                                      sapply(X = unique(sort(plot_df$bin_assignment, decreasing = TRUE)),
                                             FUN = avail_mdl,
                                             model = "MERDITH2021",
                                             ds = cdf),
                                      sapply(X = unique(sort(plot_df$bin_assignment, decreasing = TRUE)),
                                             FUN = avail_mdl,
                                             model = "PALEOMAP",
                                             ds = cdf)))
col_av_plot1 <- ggplot(data = coll_av.df1, aes(x = mid_time, y = coll_av, colour = model)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, 2000, 3000, 4000)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  ggtitle("1Myr-long bins") +
  labs(x = "Time (Ma)",
       y = "Number of available collections") +
  theme_will(plot.title = element_text(size = 18, hjust = .5),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.position = "top") +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
#save
p1 <- ggarrange(col_av_plot, col_av_plot1, nrow = 2, 
                labels = c("(a)", "(b)"), font.label = list(size = 12))
ggsave("./figures/Number_of_collections_per_model.png", p1, width = 13, height = 12)
# 4. Per-stage plot -------------------------------------------------
#Retain features of interest and filter by collection
tmp <- occdf %>% select(collection_no, bin_assignment, PALEOMAP_bin, GOLONKA_bin, MERDITH2021_bin)
tmp <- unique(tmp)
#Create collection dataframe
colldf2 <- rbind(tmp[, 1:2], tmp[, 1:2], tmp[, 1:2])
colldf2$model <- c(rep("GOLONKA", nrow(tmp)),
                  rep("PALEOMAP", nrow(tmp)),
                  rep("MERDITH2021", nrow(tmp)))
colldf2$plat_bin <- c(tmp$GOLONKA_bin, 
                      tmp$PALEOMAP_bin, 
                      tmp$MERDITH2021_bin)
#Assess number of collections per lat bin through time
colldf2 <- colldf2 %>% 
  group_by(model, bin_assignment, plat_bin) %>% 
  count(plat_bin) %>% 
  rename("occ_number" = n) %>% 
  ungroup() %>%
  complete(bin_assignment = time_bins$bin, plat_bin = lat_bins$bin, model, fill = list(occ_number = 0)) %>%
  full_join(time_bins, by = c("bin_assignment" = "bin")) %>% 
  full_join(lat_bins, by = c("plat_bin" = "bin")) %>%
  filter(!is.na(model)) %>%
  filter(!is.na(plat_bin))
#Plot
p <- ggplot(data = colldf2, aes(x = mid,
                                y = occ_number,
                                colour = model)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  # hack for the y-axis limits for the empty time bins
  geom_point(data = data.frame(x = 0, y = 4), aes(x = x, y = y), color = "transparent") +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  facet_wrap(~factor(interval_name, levels = rev(time_bins$interval_name)), ncol = 9, scales = "free") +
  coord_cartesian(ylim = c(0, NA)) +
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
