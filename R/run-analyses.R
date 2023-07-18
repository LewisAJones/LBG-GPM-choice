# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: run-analyses.R
# Last updated: 2023-03-21
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Run analyses ----------------------------------------------------------
# data prep
source("./R/subscripts/01_data_prep.R")
rm(list = ls())
# raw LBG
source("./R/subscripts/02_raw_LBG.R")
rm(list = ls())
# SQS LBG
source("./R/subscripts/03_SQS_LBG.R")
rm(list = ls())
# diversity figures
source("./R/subscripts/04_diversity_figures.R")
rm(list = ls())
# metrics figures
source("./R/subscripts/05_metrics_figures.R")
rm(list = ls())
# average GDD
source("./R/subscripts/07_Average_GDD.R")
