# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: generate-figures.R
# Aim: Generate figures
# Last updated: 2023-12-01
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Figure 1: Avg pairwise palaeolat diff and avg geodesic dist 
source("./R/subscripts/04_lat_diff_and_GDD_figure.R")
rm(list = ls())

# Figure 2 and S1: Stage-level latitudinal divrsity
source("./R/subscripts/05_Diversity_figures.R")
rm(list = ls())

# Figure 3, 4, S2, S3: Diversity and pairwise differences heatmaps
source("./R/subscripts/06_Diversity_and_Difference_heatmaps.R")
rm(list = ls())

# Figure 5, 6 and S4: Max diversity bin and rank order differences
source("./R/subscripts/07_metrics_figures.R")
rm(list = ls())

# Figure S5: Stage-level number of collections
source("./R/subscripts/08_nb_collections_figure.R")
rm(list = ls())
