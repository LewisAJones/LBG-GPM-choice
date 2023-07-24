# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 05_metrics_figures.R
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Load libraries --------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(deeptime)
library(palaeoverse)
library(AICcmodavg)
library(mgcv)
source("./R/options.R")

# Load data -------------------------------------------------------------
div_sqs <- readRDS("./data/processed/interpolations.RDS")
div_raw <- readRDS("./data/processed/genus_counts.RDS")
lat_bins <- readRDS("./data/lat_bins.RDS")
time_bins <- readRDS("./data/time_bins.RDS")

# Setup common things for figures
theme_will <- function(...) {
  theme(axis.ticks = element_line(color = "black", linewidth = 1),
        axis.line = element_blank(),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(linetype = "solid", colour = "black",
                                    fill = NA, linewidth = 2),
        ...)
}
GTS2020_periods <- time_bins(rank = "period") %>%
  rename(name = interval_name, max_age = max_ma, min_age = min_ma,
         color = colour, lab_color = font)

# Metric #1: max diversity bin  -----------------------------------------
met1_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  mutate(hemi = ifelse(mid > 0, "north", "south")) %>%
  group_by(model, stage, hemi) %>%
  filter(sum(!is.na(qD)) > 0) %>%
  summarise(max_bin = paleolat_bin[which.max(qD)], .groups = "drop") %>%
  complete(model, stage, hemi) %>%
  left_join(lat_bins, by = c("max_bin" = "bin")) %>%
  left_join(time_bins, by = c("stage" = "bin"))

met1_sqs_rects <- met1_sqs %>%
  group_by(stage, hemi, max_ma, min_ma) %>%
  summarise(all_same = length(unique(max_bin)) == 1, .groups = "drop") %>%
  filter(all_same) %>%
  mutate(ymin = ifelse(hemi == "north", 3.5, -Inf),
         ymax = ifelse(hemi == "north", Inf, 3.5))

gg_met1_sqs <- ggplot(met1_sqs, aes(x = mid_ma, y = as.numeric(factor(mid)), color = model)) +
  geom_rect(data = met1_sqs_rects, inherit.aes = FALSE,
            aes(xmin = max_ma, xmax = min_ma, ymin = ymin, ymax = ymax), fill = "grey90") +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(aes(group = interaction(hemi, model)), linewidth = .75,
            position = position_dodge(width = 2)) +
  geom_hline(yintercept = 3.5) +
  annotate(geom = "text", x = 538, y = c(0.6, 6.4),
           label = c("S. Hemisphere", "N. Hemisphere"),
           hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Most Diverse Bin", breaks = 1:6,
                     labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                     expand = expansion(add = .75)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  coord_geo(expand = TRUE, ylim = c(1, 6), dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_1.png", gg_met1_sqs, width = 13, height = 6)

met1_raw <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  mutate(hemi = ifelse(mid > 0, "north", "south")) %>%
  group_by(model, stage_bin, hemi) %>%
  filter(sum(!is.na(n_genera)) > 0) %>%
  summarise(max_bin = paleolat_bin[which.max(n_genera)], .groups = "drop") %>%
  complete(model, stage_bin, hemi) %>%
  left_join(lat_bins, by = c("max_bin" = "bin")) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

met1_raw_rects <- met1_raw %>%
  group_by(stage_bin, hemi, max_ma, min_ma) %>%
  summarise(all_same = length(unique(max_bin)) == 1, .groups = "drop") %>%
  filter(all_same) %>%
  mutate(ymin = ifelse(hemi == "north", 3.5, -Inf),
         ymax = ifelse(hemi == "north", Inf, 3.5))

gg_met1_raw <- ggplot(met1_raw, aes(x = mid_ma, y = as.numeric(factor(mid)), color = model)) +
  geom_rect(data = met1_raw_rects, inherit.aes = FALSE,
            aes(xmin = max_ma, xmax = min_ma, ymin = ymin, ymax = ymax), fill = "grey90") +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(aes(group = interaction(hemi, model)), linewidth = .75,
            position = position_dodge(width = 2)) +
  geom_hline(yintercept = 3.5) +
  annotate(geom = "text", x = 538, y = c(0.6, 6.4),
           label = c("S. Hemisphere", "N. Hemisphere"),
           hjust = 0, size = 5) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Most Diverse Bin", breaks = 1:6,
                     labels = c("High", "Middle", "Low", "Low", "Middle", "High"),
                     expand = expansion(add = .75)) +
  scale_colour_viridis_d(NULL, option = "plasma", end = .8) +
  coord_geo(expand = TRUE, ylim = c(1, 6), dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_1_raw.png", gg_met1_raw, width = 13, height = 6)

# Metric #2: rank order diffs ---------------------------------------
met2_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin), !is.na(qD)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  inner_join(., ., by = c("stage" = "stage", "paleolat_bin" = "paleolat_bin"),
             relationship = "many-to-many") %>%
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  group_by(stage, model.x, model.y) %>%
  arrange(-qD.x) %>%
  mutate(rank.x = row_number()) %>%
  arrange(-qD.y) %>%
  mutate(rank.y = row_number()) %>%
  summarise(avg = mean(abs(rank.x - rank.y)), n_bins = n(), .groups = "drop") %>%
  filter(n_bins > 1) %>%
  mutate(avg_norm = avg / c(0, 1, 4/3, 2, 2.4, 3)[n_bins]) %>%
  complete(stage = time_bins$bin, model.x, model.y) %>%
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  mutate(models = paste(model.x, model.y, sep = "/")) %>%
  left_join(time_bins, by = c("stage" = "bin"))

gg_met2_sqs <- ggplot(met2_sqs, aes(x = mid_ma, y = avg_norm, color = models, group = models)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = .75, position = position_dodge(width = 2)) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Norm. Avg. Rank Order Diff.", limits = c(0, 1)) +
  scale_colour_viridis_d(NULL, end = .9) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_2.png", gg_met2_sqs, width = 13, height = 6)

met2_sqs %>%
  nest_by(models) %>%
  mutate(mod = list(lm(avg_norm ~ mid_ma, data = data))) %>% 
  summarize(glance(mod))

ggsave("./figures/metric_2_reg.png",
       gg_met2_sqs + geom_smooth(method = "lm"),
       width = 13, height = 6)

met2_raw <- div_raw %>%
  filter(!is.na(paleolat_bin), !is.na(n_genera), n_genera > 0) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  inner_join(., ., by = c("stage_bin" = "stage_bin", "paleolat_bin" = "paleolat_bin"),
             relationship = "many-to-many") %>%
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  group_by(stage_bin, model.x, model.y) %>%
  arrange(-n_genera.x) %>%
  mutate(rank.x = row_number()) %>%
  arrange(-n_genera.y) %>%
  mutate(rank.y = row_number()) %>%
  summarise(avg = mean(abs(rank.x - rank.y)), n_bins = n(), .groups = "drop") %>%
  filter(n_bins > 1) %>%
  mutate(avg_norm = avg / c(0, 1, 4/3, 2, 2.4, 3)[n_bins]) %>%
  complete(stage_bin = time_bins$bin, model.x, model.y) %>%
  filter((model.x == "PALEOMAP" & model.y == "GOLONKA") |
           (model.x == "PALEOMAP" & model.y == "MERDITH2021") |
           (model.x == "GOLONKA" & model.y == "MERDITH2021")) %>%
  mutate(models = paste(model.x, model.y, sep = "/")) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

gg_met2_raw <- ggplot(met2_raw, aes(x = mid_ma, y = avg_norm, color = models, group = models)) +
  geom_point(size = 1.5, position = position_dodge(width = 2)) +
  geom_line(linewidth = .75, position = position_dodge(width = 2)) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Norm. Avg. Rank Order Diff.", limits = c(0, 1)) +
  scale_colour_viridis_d(NULL, end = .9) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 20) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_2_raw.png", gg_met2_raw, width = 13, height = 6)

ggsave("./figures/metric_2_raw_reg.png",
       gg_met2_raw + geom_smooth(method = "lm"),
       width = 13, height = 6)

# Combine #1 and #2 -------------------------------------------------
ggsave("./figures/metrics_1_2_sqs.png",
       ggarrange2(gg_met1_sqs, gg_met2_sqs, ncol = 1, draw = FALSE), width = 13, height = 12)

ggsave("./figures/metrics_1_2_raw.png",
       ggarrange2(gg_met1_raw, gg_met2_raw, ncol = 1, draw = FALSE), width = 13, height = 12)

# Metric #3: sum of squares -----------------------------------------
met3_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  pivot_wider(id_cols = c(stage, paleolat_bin), names_from = model, values_from = qD) %>%
  group_by(stage) %>%
  summarise(pg_sum = sum(abs(PALEOMAP - GOLONKA), na.rm = TRUE),
            pg_n = sum(!is.na(PALEOMAP) & !is.na(GOLONKA)),
            pm_sum = sum(abs(PALEOMAP - MERDITH2021), na.rm = TRUE),
            pm_n = sum(!is.na(PALEOMAP) & !is.na(MERDITH2021)),
            gm_sum = sum(abs(GOLONKA - MERDITH2021), na.rm = TRUE),
            gm_n = sum(!is.na(GOLONKA) & !is.na(MERDITH2021)),
            .groups = "drop") %>%
  pivot_longer(cols = c(pg_sum, pg_n, pm_sum, pm_n, gm_sum, gm_n)) %>%
  separate_wider_delim(name, "_", names = c("models", "name")) %>%
  pivot_wider(id_cols = c(stage, models)) %>%
  mutate(avg = sum/n) %>%
  left_join(time_bins, by = c("stage" = "bin"))

gg_met3_sqs <- ggplot(met3_sqs) +
  geom_line(aes(x = mid_ma, y = avg, color = models, group = models), linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Average Diversity Estimate Difference", trans = "log10") +
  scale_colour_viridis_d(NULL, end = .9, labels = c("GOLONKA/MERDITH2021",
                                                    "PALEOMAP/GOLONKA",
                                                    "PALEOMAP/MERDITH2021")) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_3.png", gg_met3_sqs, width = 13, height = 4.5)

met3_raw <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  pivot_wider(id_cols = c(stage_bin, paleolat_bin), names_from = model, values_from = n_genera) %>%
  group_by(stage_bin) %>%
  summarise(pg_sum = sum(abs(PALEOMAP - GOLONKA), na.rm = TRUE),
            pg_n = sum(!is.na(PALEOMAP) & !is.na(GOLONKA)),
            pm_sum = sum(abs(PALEOMAP - MERDITH2021), na.rm = TRUE),
            pm_n = sum(!is.na(PALEOMAP) & !is.na(MERDITH2021)),
            gm_sum = sum(abs(GOLONKA - MERDITH2021), na.rm = TRUE),
            gm_n = sum(!is.na(GOLONKA) & !is.na(MERDITH2021)),
            .groups = "drop") %>%
  pivot_longer(cols = c(pg_sum, pg_n, pm_sum, pm_n, gm_sum, gm_n)) %>%
  separate_wider_delim(name, "_", names = c("models", "name")) %>%
  pivot_wider(id_cols = c(stage_bin, models)) %>%
  mutate(avg = sqrt(sum)/n) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

gg_met3_raw <- ggplot(met3_raw) +
  geom_line(aes(x = mid_ma, y = avg, color = models, group = models), linewidth = .75) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0), expand = expansion()) +
  scale_y_continuous("Average Raw Diversity Difference", trans = "log10") +
  scale_colour_viridis_d(NULL, end = .9, labels = c("GOLONKA/MERDITH2021",
                                                    "PALEOMAP/GOLONKA",
                                                    "PALEOMAP/MERDITH2021")) +
  coord_geo(expand = TRUE, dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_3_raw.png", gg_met3_raw, width = 13, height = 4.5)

# Metric #4: best fit model -----------------------------------------
met4_sqs <- div_sqs %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage, paleolat_bin, model, qD) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  group_by(stage, model) %>%
  summarise(poly = which.min(sapply(1:5, function(i) {
    tryCatch(AICc(glm(qD ~ poly(mid, i))),
             error = function(e) NA)
  }))) %>%
  left_join(time_bins, by = c("stage" = "bin"))

gg_met4_sqs <- ggplot(met4_sqs) +
  geom_tile(aes(x = mid_ma, y = model, width = max_ma - min_ma, height = 1, fill = as.factor(poly))) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0)) +
  scale_y_discrete(NULL) +
  scale_fill_viridis_d(NULL, end = .9) +
  coord_geo(dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_4.png", gg_met4_sqs, width = 13, height = 4.5)

met4_raw <- div_raw %>%
  filter(!is.na(paleolat_bin)) %>%
  select(stage_bin, paleolat_bin, model, n_genera) %>%
  left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
  group_by(stage_bin, model) %>%
  summarise(poly = which.min(sapply(1:5, function(i) {
    tryCatch(AICc(glm(n_genera ~ poly(mid, i))),
             error = function(e) NA)
  }))) %>%
  left_join(time_bins, by = c("stage_bin" = "bin"))

gg_met4_raw <- ggplot(met4_raw) +
  geom_tile(aes(x = mid_ma, y = model, width = max_ma - min_ma, height = 1, fill = as.factor(poly))) +
  scale_x_reverse("Time (Ma)", limits = c(541, 0)) +
  scale_y_discrete(NULL) +
  scale_fill_viridis_d(NULL, end = .9) +
  coord_geo(dat = GTS2020_periods, lwd = 1,
            bord = c("left", "right", "bottom")) +
  theme_classic(base_size = 14) +
  theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
ggsave("./figures/metric_4_raw.png", gg_met4_raw, width = 13, height = 4.5)

# Metric #4b: GAM knots -----------------------------------------
# met4b_sqs <- div_sqs %>%
#   filter(!is.na(paleolat_bin)) %>%
#   select(stage, paleolat_bin, model, qD) %>%
#   left_join(lat_bins, by = c("paleolat_bin" = "bin")) %>%
#   group_by(stage, model) %>%
#   summarise(edf = 
#     tryCatch(summary(gam(qD ~ s(mid)))$edf,
#              error = function(e) NA)
#   ) %>%
#   left_join(time_bins, by = c("stage" = "bin"))
# 
# gg_met4b_sqs <- ggplot(met4b_sqs) +
#   geom_tile(aes(x = mid_ma, y = model, width = max_ma - min_ma, height = 1, fill = edf)) +
#   scale_x_reverse("Time (Ma)", limits = c(541, 0)) +
#   scale_y_discrete(NULL) +
#   scale_fill_viridis_c(NULL, end = .9) +
#   coord_geo(dat = GTS2020_periods, lwd = 1,
#             bord = c("left", "right", "bottom")) +
#   theme_classic(base_size = 14) +
#   theme_will(legend.position = "top", legend.margin = margin(-5, -5, -5, -5))
# ggsave("./figures/metric_4b.png", gg_met4b_sqs, width = 13, height = 4.5)


## Sampling through time metric ---------------------------------
  #Load data
occdf <- readRDS("./data/processed/pbdb_data.RDS")
colldf <- occdf %>% select(collection_no, bin_assignment)
colldf <- unique(colldf)
row.names(colldf) <- 1:nrow(colldf)
  #Evaluate number of collection through time
n_col <- unlist(lapply(X = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                       FUN = function(x){
                         idx <- which(colldf$bin_assignment == x)
                         return(length(idx))
                       }))

nb_coll.df <- data.frame(time_bin = unique(sort(colldf$bin_assignment, decreasing = TRUE)),
                         number_of_collections = n_col)
nb_coll.df$mid_time <- unlist(lapply(X = nb_coll.df$time_bin,
                                     FUN = function(x){
                                       return(time_bins$mid_ma[which(time_bins$bin == x)])
                                     }))

col_plot <- ggplot(data = nb_coll.df, aes(x = mid_time, y = number_of_collections)) +
  scale_x_reverse(limits = c(542, -0.7),
                  breaks = c(0, 100, 200, 300, 400, 500),
                  labels = c(0, 100, 200, 300, 400, 500)) +
  scale_y_continuous(limits = c(0, 4100),
                     breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, 2000, 3000, 4000)) +
  geom_point(size = 2, colour = "blue") +
  geom_line(linewidth = 1, colour = "blue") +
  labs(x = "Time (Ma)",
       y = "Number of collections") +
  theme_will(axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14)) +
  coord_geo(lwd = 1)
ggsave("./figures/Number_of_collections.png", col_plot, width = 13, height = 6)
