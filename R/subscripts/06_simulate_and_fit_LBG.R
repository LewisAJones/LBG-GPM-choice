# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_simulation_fitting.R
# Last updated: 2023-03-27
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Libraries -------------------------------------------------------------
library(dplyr)

# Simulate biodiversity gradients ---------------------------------------
# emp_div is the empirical diversity data (at a given time) we will extract the simulation parameters from
# n: Number of data points to generate  
# div_column_key is the key of the diversity column in the emp_div dataset
# emp_lats is the equal-area lat bins dataframe used in the empirical diversiy assessment
# mid_lats_key is the key of the mid latitude in equal-area lat bins dataframe
# lat_bin_keys is the key of the latitudinal bin column in both empirical lats and diversity datasets (HAS TO BE THE SAME)
# type is the type of gradient to simulate: "unimodal", "bimodal", "flat" or "MIXED"
# if type="MIXED", the user has to specify the shape of the northern and southern gradients to simulate: "gaussian" or "flat"
simulate_biodiv_grad <- function(emp_data, n, div_column_key, emp_lats, mid_lat_key, lat_bin_keys, type, northern=NA, southern=NA) {
  # Sequence of latitudes used and hemisphere specification
  latitude <- seq(-90, 90, length.out = n)
  if(0 %in% latitude){ #easier to partition northern and southern hemispheres
    latitude <- latitude[-which(latitude == 0)]
  }  
  # Extract mean and sd from emp_data
  ## Global scale
  meanDIV <- mean(emp_data[, c(div_column_key)], na.rm = TRUE)
  sdDIV <- sd(emp_data[, c(div_column_key)], na.rm = TRUE)
  ## Northern hemisphere
  northern_bins <- emp_lats[which(emp_lats[, c(mid_lat_key)] > 0), c(lat_bin_key)]
  northern_subset <- emp_data[which(emp_data[, c(lat_bin_key)] %in% northern_bins), ]
  meanDIV_North <- mean(northern_subset[, c(div_column_key)], na.rm = TRUE)
  sdDIV_North <- sd(northern_subset[, c(div_column_key)], na.rm = TRUE)
  ## Southern hemisphere
  southern_bins <- emp_lats[which(emp_lats[, c(mid_lat_key)] < 0), c(lat_bin_key)]
  southern_subset <- emp_data[which(emp_data[, c(lat_bin_key)] %in% southern_bins), ]
  meanDIV_South <- mean(southern_subset[, c(div_column_key)], na.rm = TRUE)
  sdDIV_South <- sd(southern_subset[, c(div_column_key)], na.rm = TRUE)
  
  # Simulate a unimodal gradient using a gaussian distribution
  if (type == "unimodal") {
    species_richness <- simulate_gaussian(latitude = latitude,
                                          mean = meanDIV,
                                          sd = sdDIV)
  } 
  # Simulate a bimodal gradient using a combination of two Gaussian distributions with different means and sds extracted from the data
  if (type == "bimodal") {
    ## Northern hemisphere
    species_richness <- dnorm(latitude[which(latitude > 0)],
                              mean = meanDIV_North, 
                              sd = sdDIV_North)
    ## Southern hemisphere
    species_richness <- species_richness + dnorm(latitude[which(latitude < 0)], 
                                                 mean = meanDIV_South, 
                                                 sd = sdDIV_South)
  } 
  # Simulate a flat gradient
  if (type == "flat") {
    ## Set the species richness to a constant value, mean diversity
    species_richness <- rep(meanDIV, n)
  }
  # Dissociate both hemispheres
  if (type == "MIXED") {
    if ( (northern == "gaussian") & (southern == "flat") ){
      ## Northern hemisphere
      species_richness <- dnorm(latitude[which(latitude > 0)], 
                                mean = meanDIV_North, 
                                sd = sdDIV_North)
      ## Southern hemisphere
      species_richness <- species_richness + rep(meanDIV_South, n)
    }
    else if ( (northern == "flat") & (southern == "gaussian") ){
      ## Northern hemisphere
      species_richness <- rep(meanDIV_North, n)
      ## Southern hemisphere
      species_richness <- species_richness + dnorm(latitude[which(latitude < 0)], 
                                                   mean = meanDIV_South, 
                                                   sd = sdDIV_South)
    }
    else(stop("northern and southern arguments have to either be set to 'gaussian' or 'flat'. If set to the same string, please refer to other 'type' arguments"))
  }
  # Proportional richness
  species_richness <- species_richness / sum(species_richness)
  # Return the simulated data as a dataframe
  return(data.frame(latitude = latitude, species_richness = species_richness))
}
# Examples
# Empirical data
div_raw <- readRDS("./data/processed/genus_counts.RDS")
lats <- readRDS("./data/lat_bins.RDS")
lats <- rename(lats, paleolat_bin = "bin")
# Simulate gradients based on the 20th time stage and PALEOMAP 
unimodal_data <- simulate_biodiv_grad(emp_data = div_raw[which((div_raw$stage_bin == 20)&(div_raw$model == "PALEOMAP")), ],
                                      n = 1000,
                                      div_column_key = "n_genera",
                                      mid_lat_key = "mid",
                                      emp_lats = lats,
                                      lat_bin_keys = "paleolat_bin",
                                      type = "unimodal")

bimodal_data <- simulate_biodiv_grad(emp_data = div_raw[which((div_raw$stage_bin == 20)&(div_raw$model == "PALEOMAP")), ],
                                     n = 1000,
                                     div_column_key = "n_genera",
                                     mid_lat_key = "mid", 
                                     emp_lats = lats,
                                     lat_bin_keys = "paleolat_bin",
                                     type = "bimodal")

flat_data <- simulate_biodiv_grad(emp_data = div_raw[which((div_raw$stage_bin == 20)&(div_raw$model == "PALEOMAP")), ],
                                  n = 1000,
                                  div_column_key = "n_genera",
                                  mid_lat_key = "mid",
                                  emp_lats = lats,
                                  lat_bin_keys = "paleolat_bin",
                                  type = "flat")
# Plot
plot(unimodal_data)
plot(bimodal_data)
plot(flat_data)

# Compare curves using root mean squared error (RMSE) 
# empirical_data is our reconstructed gradient
# simulated_data is the simulated gradient
# x_col is the name of our x column
# y_col is the name of our y column
compare_curves <- function(empirical_data, simulated_data, x_col = "latitude", y_col = "species_richness") {
  # Get the x and y values from the empirical data
  x_emp <- empirical_data[[x_col]]
  y_emp <- empirical_data[[y_col]]
  # Interpolate the simulated data onto the x values of the empirical data
  y_sim <- approx(simulated_data[[x_col]], simulated_data[[y_col]], xout = x_emp)$y
  # Calculate the root mean squared error (RMSE) between the empirical and simulated data
  rmse <- sqrt(mean((y_emp - y_sim)^2))
  # Rescale RMSE to the number of data points (this is unneccesary, it is just there
  # so we aren't dealing with tiny numbers for now)
  rmse <- rmse * nrow(simulated_data)
  # Return the RMSE value
  return(rmse)
}

# Compare the simulated data to the empirical data (here we are comparing two
# simulated gradients but you get the idea)
rmse <- compare_curves(unimodal_data, bimodal_data)
rmse