# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_simulation_fitting.R
# Last updated: 2023-03-27
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Simulate biodiversity gradients ---------------------------------------
# n is the number of data points to generate
# type is the type of gradient to simulate: "unimodal", "bimodal", or "flat"
simulate_biodiv_grad <- function(n, type) {
  # Generate a sequence of latitude values ranging from -90 to 90 degrees
  latitude <- seq(-90, 90, length.out = n)
  # Initialize the species richness vector
  species_richness <- rep(0, n)
  # Simulate a unimodal gradient
  if (type == "unimodal") {
    # Use a Gaussian distribution to simulate a unimodal gradient
    species_richness <- 5*dnorm(latitude, mean = 0, sd = 30)
  } 
  # Simulate a bimodal gradient
  if (type == "bimodal") {
    # Use a combination of two Gaussian distributions with different means
    species_richness <- 5*dnorm(latitude, mean = 45, sd = 15) + 
      5*dnorm(latitude, mean = -45, sd = 15)
  } 
  # Simulate a flat gradient
  if (type == "flat") {
    # Set the species richness to a constant value
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