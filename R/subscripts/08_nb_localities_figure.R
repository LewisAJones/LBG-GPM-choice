# Header -----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_nb_collections_figure.R
# Aim: Plot number of collections through time and binned in the different 
#      palaeolat bins per stage (Fig.S5).
# Last updated: 2024-12-17
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
ics_periods <- time_bins(scale = "international periods") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
ics_eras <- time_bins(scale = "international eras") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)
#Load data
time_bins <- readRDS("./data/time_bins.RDS")
lat_bins <- readRDS("./data/lat_bins.RDS")
occdf <- readRDS("./data/processed/pbdb_data.RDS")
pointdf <- occdf %>%
  distinct(lat, lng, bin_assignment, .keep_all = TRUE) %>%
  select(bin_assignment, p_lat_PALEOMAP, p_lat_GOLONKA, p_lat_MERDITH2021, p_lat_TorsvikCocks2017, bin_midpoint)

#Total number of collection through time
n_points <- pointdf %>% group_by(bin_assignment) %>% count() %>%
  arrange(desc(bin_assignment))
nb_points.df <- data.frame(mid_time = unique(sort(pointdf$bin_midpoint,
                                                decreasing = FALSE)), #increasing order as the oldest bin has the smallest bin assignment number
                         number_of_points = n_points$n,
                         group = rep("Number of localities", nrow(n_points)))
#Number of available collections per GPM relative to the stage length
stage_duration <- sapply(X = unique(sort(n_points$bin_assignment, decreasing = TRUE)),
                         FUN = function(x){
                         return(time_bins$duration_myr[which(time_bins$bin == x)])
                         })
points_per_My <- nb_points.df$number_of_points / stage_duration
points_per_My[1] <- 900 #restrict nb collections per million years of the Holocene for plotting issues
#Total dataframe
DF <- rbind.data.frame(nb_points.df,
                       data.frame(mid_time = unique(sort(pointdf$bin_midpoint,
                                                         decreasing = FALSE)), #increasing order as the oldest bin has the smallest bin assignment number
                                  number_of_points = points_per_My,
                                  group = rep("Number of localities per My",
                                              nrow(n_points))))
#Plot
point_plot <- ggplot(data = DF, aes(x = mid_time, y = number_of_points)) +
  scale_x_reverse(breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  geom_point(size = 2, colour = "#e7298a") +
  geom_line(linewidth = 1, colour = "#e7298a") +
  labs(x = "Time (Ma)",
       y = NULL) +
  coord_geo(list("bottom", "bottom"), expand = FALSE, dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE)) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "none") +
  facet_wrap(~group, ncol = 1, scales = "free_y")
#save
ggsave("./figures/Number_of_localities_total.png", point_plot, width = 13, height = 12)

# 2. Total number of collections available for each model ----------------------
point_av.df <- pointdf %>%
  #pivot across the three palaeolat values
  pivot_longer(starts_with("p_lat_"), names_to = "model", names_prefix = "p_lat_") %>%
  group_by(bin_assignment, bin_midpoint, model) %>%
  #sum the number of available collections per time bins
  summarise(point_av = sum(!is.na(value))) %>% 
  mutate(group = rep("TOTAL", 4))


#assess number of available collections per GPM relative to the stage length
point_av.df2 <- point_av.df %>% 
  mutate(group = rep("PER_My", 4))

stage_duration <- function(bin){
  return(time_bins$duration_myr[which(time_bins$bin == bin)])
}
point_av.df2$duration <- sapply(X = point_av.df2$bin_assignment,
                                FUN = stage_duration)

point_av.df2 <- point_av.df2 %>% 
  mutate(point_av = point_av / duration) %>%
  select(-duration) %>% 
  bind_rows(point_av.df)
  #restrict holocene for plotting issues
point_av.df2$point_av[which(point_av.df2$bin_assignment == 97 & 
                            point_av.df2$group == "PER_My")] <- 900
  #Customise group label
point_av.df2$group <- factor(point_av.df2$group, levels = c("TOTAL", "PER_My"),
                  labels = c("Total number of sampled localities",
                             "Number of sampled localities per million years")
)

#plot
point_av_plot <- ggplot(data = point_av.df2, aes(x = bin_midpoint, y = point_av,
                                                 colour = model)) +
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
  coord_geo(list("bottom", "bottom"), dat = list(ics_eras, ics_periods),
            lwd = 1, bord = c("left", "right", "bottom"), abbrv = list(FALSE, TRUE),
            size = 6) +
  facet_wrap(~group, ncol = 1, scales = "free_y")
#save
ggsave("./figures/Number_of_localities_per_model.png", point_av_plot, width = 13, height = 12)

# 3. Per-stage plot -------------------------------------------------
#Retain features of interest and filter by collection
tmp <- occdf %>% select(lng, lat, bin_assignment, PALEOMAP_bin,
                        GOLONKA_bin, MERDITH2021_bin, TorsvikCocks2017_bin)
tmp <- unique(tmp)
#Create collection dataframe
pointdf2 <- rbind(tmp[, 1:3], tmp[, 1:3], tmp[, 1:3], tmp[, 1:3])
pointdf2$model <- c(rep("GOLONKA", nrow(tmp)),
                  rep("PALEOMAP", nrow(tmp)),
                  rep("MERDITH2021", nrow(tmp)),
                  rep("TorsvikCocks2017", nrow(tmp)))
pointdf2$plat_bin <- c(tmp$GOLONKA_bin, 
                      tmp$PALEOMAP_bin, 
                      tmp$MERDITH2021_bin,
                      tmp$TorsvikCocks2017_bin)
#Assess number of collections per lat bin through time
pointdf2 <- pointdf2 %>% 
  group_by(model, bin_assignment, plat_bin) %>% 
  count(plat_bin) %>% 
  rename("loc_number" = n) %>% 
  ungroup() %>%
  complete(bin_assignment = time_bins$bin, plat_bin = lat_bins$bin, model,
           fill = list(loc_number = 0)) %>%
  full_join(time_bins, by = c("bin_assignment" = "bin")) %>% 
  full_join(lat_bins, by = c("plat_bin" = "bin")) %>%
  filter(!is.na(model)) %>%
  filter(!is.na(plat_bin))
#Plot
p <- ggplot(data = pointdf2, aes(x = mid,
                                y = loc_number,
                                colour = model)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  # hack for the y-axis limits for the empty time bins
  geom_point(data = data.frame(x = 0, y = 4), aes(x = x, y = y),
             color = "transparent") +
  geom_line(linewidth = 0.75, alpha = 1, position = position_dodge(width = 2)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  facet_wrap(~factor(interval_name, levels = rev(time_bins$interval_name)),
             ncol = 9, scales = "free") +
  coord_cartesian(ylim = c(0, NA)) +
  labs(y = "Number of localities",
       x = "Palaeolatitude (º)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        strip.text = element_text(size = 12,
                                  margin = margin(3.3, 4.4, 3.3, 4.4)),
        strip.background = element_rect(size = 1))

#Update strip colours
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
for (i in strip_t) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <-
    time_bins$colour[match(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label,
                           time_bins$interval_name)]
}
#Save
ggsave("./figures/Nb_localities_per_stage.png", g, width = 20, height = 22)
  