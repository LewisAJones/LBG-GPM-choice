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
#assess number of available collections per GPM relative to the stage length
nb_coll.df$stage_duration <- sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                                    FUN = function(x){
                                      return(time_bins$duration_myr[which(time_bins$bin == x)])
                                    })
nb_coll.df$coll_per_My <- nb_coll.df$number_of_collections / nb_coll.df$stage_duration

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
  theme_will(plot.title = element_text(size = 18, hjust = .5),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
#collections/Myr
col_pm_plot <- ggplot(data = nb_coll.df, aes(x = mid_time, y = coll_per_My)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 900),
                     breaks = c(0, 200, 400, 600, 800),
                     labels = c(0, 200, 400, 600, 800)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  labs(x = "Time (Ma)",
       y = "Nb. collections per My") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
#save
p <- ggarrange(col_plot, col_pm_plot, nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 20))
ggsave("./figures/Number_of_collections_total.png", p, width = 13, height = 12)

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
#plot
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
  labs(x = "Time (Ma)",
       y = "Number of available collections") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.position = "top") +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
#assess number of available collections per GPM relative to the stage length
coll_av.df$stage_duration <- sapply(X = rep(unique(sort(colldf$bin_assignment, decreasing = TRUE)), 3),
                                    FUN = function(x){
                                      return(time_bins$duration_myr[which(time_bins$bin == x)])
                                    })
coll_av.df$coll_av_pm <- coll_av.df$coll_av / coll_av.df$stage_duration
#plot
col_av_pm_plot <- ggplot(data = coll_av.df, aes(x = mid_time, y = coll_av_pm, colour = model)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 900),
                     breaks = c(0, 200, 400, 600, 800),
                     labels = c(0, 200, 400, 600, 800)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  labs(x = "Time (Ma)",
       y = "Nb. available collections per Myr") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.position = "top") +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE))
#save
p <- ggarrange(col_av_plot, col_av_pm_plot, nrow = 2, labels = c("(a)", "(b)"), font.label = list(size = 20))
ggsave("./figures/Number_of_collections_per_model.png", p, width = 13, height = 12)

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
