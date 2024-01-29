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
                         number_of_collections = n_col, group = rep("Number of collections", length(n_col)))
#Number of available collections per GPM relative to the stage length
stage_duration <- sapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                         FUN = function(x){
                         return(time_bins$duration_myr[which(time_bins$bin == x)])
                         })
coll_per_My <- nb_coll.df$number_of_collections / stage_duration
coll_per_My[1] <- 900 #restrict nb collections per million years of the Holocene for plotting issues
#Total dataframe
DF <- rbind.data.frame(nb_coll.df,
                       data.frame(mid_time = unique(sort(colldf$bin_midpoint, decreasing = FALSE)), #increasing order as the oldest bin has the smallest bin assignment number
                                  number_of_collections = coll_per_My, group = rep("Number of collections per My", length(n_col))))
#Plot
col_plot <- ggplot(data = DF, aes(x = mid_time, y = number_of_collections)) +
  scale_x_reverse(breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  labs(x = "Time (Ma)",
       y = NULL) +
  coord_geo(list("bottom", "bottom"), expand = FALSE, dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "none") +
  facet_wrap(~group, ncol = 1, scales = "free_y")
#save
ggsave("./figures/Number_of_collections_total.png", col_plot, width = 13, height = 12)

# 2. Total number of collections available for each model ----------------------
coll_av.df <- colldf %>%
  #pivot across the three palaeolat values
  pivot_longer(starts_with("p_lat_"), names_to = "model", names_prefix = "p_lat_") %>%
  group_by(bin_assignment, bin_midpoint, model) %>%
  #sum the number of available collections per time bins
  summarise(coll_av = sum(!is.na(value))) %>% 
  mutate(group = rep("TOTAL", 3))


#assess number of available collections per GPM relative to the stage length
coll_av.df2 <- coll_av.df %>% 
  mutate(group = rep("PER_My", 3))

stage_duration <- function(bin){
  return(time_bins$duration_myr[which(time_bins$bin == bin)])
}
coll_av.df2$duration <- sapply(X = coll_av.df2$bin_assignment, FUN = stage_duration)

coll_av.df2 <- coll_av.df2 %>% 
  mutate(coll_av = coll_av / duration) %>%
  select(-duration) %>% 
  bind_rows(coll_av.df)
  #restrict holocene for plotting issues
coll_av.df2$coll_av[which(coll_av.df2$bin_assignment == 97 & 
                            coll_av.df2$group == "PER_My")] <- 900
  #Customise group label
coll_av.df2$group <- factor(coll_av.df2$group, levels = c("TOTAL", "PER_My"),
                  labels = c("Total number of available collections","Number of available collections per million years")
)

#plot
coll_av_plot <- ggplot(data = coll_av.df2, aes(x = bin_midpoint, y = coll_av, colour = model)) +
  scale_x_reverse(breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  geom_point(size = 1.75, position = position_dodge(width = 2)) +
  geom_line(linewidth = 1, alpha = 1, position = position_dodge(width = 2)) +
  labs(x = "Time (Ma)",
       y = NULL) +
  theme_classic(base_size = 24) +
  theme_will(axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 18),
             axis.text.y = element_text(size = 18),
             legend.position = "top") +
  coord_geo(list("bottom", "bottom"), dat = list(GTS2020_eras, GTS2020_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE),
            size = 6) +
  facet_wrap(~group, ncol = 1, scales = "free_y")
#save
ggsave("./figures/Number_of_collections_per_model.png", coll_av_plot, width = 13, height = 12)

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
