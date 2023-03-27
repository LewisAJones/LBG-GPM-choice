# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_simulation_fitting.R
# Last updated: 2023-03-27
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Simulate biodiversity gradients ---------------------------------------
# emp_div is the empirical diversity data (at a given time) we will extract the simulation parameters from
# div_column_key is the key of the diversity column in the emp_div dataset
# emp_lats is the equal-area lat bins dataframe used in the empirical diversiy assessment
# mid_lats_key is the key of the mid latitude in equal-area lat bins dataframe
# lat_bin_keys is the key of the latitudinal bin column in both empirical lats and diversity datasets (HAS TO BE THE SAME)
# type is the type of gradient to simulate: "unimodal", "bimodal", "flat" or "MIXED"
# if type="MIXED", the user has to specify the shape of the northern and southern gradients to simulate: "gaussian" or "flat"
simulate_biodiv_grad <- function(emp_data, div_column_key, emp_lats, mid_lats_key, lat_bin_keys, type, northern=NA, southern=NA) {
  # Number of data points to generate  
  n <- nrow(emp_lats)
  # Sequence of latitudes used end hemisphere specification
  latitude <- emp_lats[, c(mid_lats_key)]
  northern_hemisphere <- emp_lats[which(latitude > 0), c(lat_bin_keys)]
  southern_hemisphere <- emp_lats[which(latitude < 0), c(lat_bin_keys)]
  # Initialise the species richness vector
  species_richness <- rep(0, n)
  # Simulate a unimodal gradient
  if (type == "unimodal") {
    meanDIV <- mean(emp_data[, c(div_column_key)], na.rm = TRUE)
    sdDIV <- sd(emp_data[, c(div_column_key)])
    # Use a Gaussian distribution to simulate a unimodal gradient
    species_richness <- dnorm(latitude, mean = meanDIV, sd = sdDIV)
  } 
  # Simulate a bimodal gradient using a combination of two Gaussian distributions with different means and sds extracted from the data
  if (type == "bimodal") {
    #Northern hemisphere
    northern_bins <- which(emp_data[, c(lat_bin_keys)] %in% northern_hemisphere)
    meanDIV_North <- mean(emp_data[northern_bins, c(div_column_key)], na.rm = TRUE)
    sdDIV_North <- sd(emp_data[northern_bins, c(div_column_key)])
    species_richness <- dnorm(latitude[which(latitude > 0)], mean = meanDIV_North, sd = sdDIV_North)
    #Southern hemisphere
    southern_bins <- which(emp_data[, c(lat_bin_keys)] %in% southern_hemisphere)
    meanDIV_South <- mean(emp_data[southern_bins, c(div_column_key)], na.rm = TRUE)
    sdDIV_South <- sd(emp_data[southern_bins, c(div_column_key)])
    species_richness <- species_richness + dnorm(latitude[which(latitude > 0)], mean = meanDIV_South, sd = sdDIV_South)
  } 
  # Simulate a flat gradient
  if (type == "flat") {
    # Set the species richness to a constant value, set to mean diversity
    species_richness <- rep(10, n)
  }
  # Proportional richness
  species_richness <- species_richness / sum(species_richness)
  # Return the simulated data as a dataframe
  return(data.frame(latitude = latitude, species_richness = species_richness))
}
# Examples
# Simulate gradients
unimodal_data <- simulate_biodiv_grad(n = 1000, type = "unimodal")
bimodal_data <- simulate_biodiv_grad(n = 1000, type = "bimodal")
flat_data <- simulate_biodiv_grad(n = 1000, type = "flat")
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