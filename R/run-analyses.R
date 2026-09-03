# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: run-analyses.R
# Last updated: 2026-09-03
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Run analyses ----------------------------------------------------------
# data prep
source("R/subscripts/01_data_prep.R")
rm(list = ls())
# raw LBG
source("R/subscripts/02_raw_LBG.R")
rm(list = ls())
# SQS LBG
source("R/subscripts/03_SQS_LBG.R")
rm(list = ls())
# Latitudinal difference figures
source("R/subscripts/04_lat_diff_figure.R")
rm(list = ls())
# Diversity figures
source("R/subscripts/05_diversity_figures.R")
rm(list = ls())
# Diversity and difference heatmaps
source("R/subscripts/06_diversity_and_difference_heatmaps.R")
rm(list = ls())
# Metrics figures
source("R/subscripts/07_metrics_figures.R")
rm(list = ls())
# Number of collections through time figures
source("R/subscripts/08_nb_localities_figure.R")
rm(list = ls())
# Summary stats
source("R/subscripts/09_summary_stats.R")
rm(list = ls())
