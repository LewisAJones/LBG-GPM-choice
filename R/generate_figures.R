# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: generate-figures.R
# Aim: Generate figures
# Last updated: 2023-12-01
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Figure 1: Avg pairwise palaeolat diff and avg geodesic dist 
source("./R/subscripts/04_lat_diff_and_GDD_figure.R")
rm(list = ls())

# Figure 3 and S2: Stage-level latitudinal divrsity
source("./R/subscripts/05_Diversity_figures.R")
rm(list = ls())

# Figure 4, 5, S3, S4 and S5: Diversity and pairwise differences heatmaps
source("./R/subscripts/06_Diversity_and_Difference_heatmaps.R")
rm(list = ls())

# Figure 5, 6, S6 and S7: Max diversity bin and rank order differences
source("./R/subscripts/07_metrics_figures.R")
rm(list = ls())

# Figure S1: Stage-level number of collections
source("./R/subscripts/08_nb_localities_figure.R")
rm(list = ls())
