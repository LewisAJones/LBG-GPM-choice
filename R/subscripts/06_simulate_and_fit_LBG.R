# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: 06_simulation_fitting.R
# Last updated: 2023-03-27
# Repository: https://github.com/LewisAJones/LBG-GPM-choice

# Libraries -------------------------------------------------------------
library(dplyr)

## Sigma finder ---------------------------------------------------------
  # Unimodal case
sigma_unimodal <- function(data, mu = 0){ #data is a vector of length 6, containing estimates of diversity for one time bin and one model
  #"tropical diversity", i.e, between -19.48 and 19.48° lat (bin 3 and 4)
  div_trop <- sum(data[3:4])/sum(data)
  #Find sigma
  sigma_prop <- sigma_prop <- seq(0.01, 1000, 0.01)
  trop_estimated <- function(sigma){
    capricorn <- pnorm(q = -19.48, mean = mu, sd = sigma)
    cancer <- pnorm(q = 19.48, mean = mu, sd = sigma)
    return(cancer - capricorn)
  }
  CDF_emp <- unlist(lapply(X = sigma_prop,
                           FUN = trop_estimated))
  CDF_emp <- abs(CDF_emp - div_trop)
  sigma_emp <- sigma_prop[which.min(CDF_emp)]
  return(sigma_emp)
}
  #Bimodal_case
sigma_bimodal <- function(data, loc){ #loc = "NORTH" or "SOUTH"
  mu <- 45
  if(loc == "NORTH"){
    sub_div <- data[1:3]
    div_2 <- data[1]
  }
  else if(loc == "SOUTH"){
    sub_div <- data[4:6]
    div_2 <- data[4]
  }
  else{
    stop("Argument 'loc' unspecified.")
  }
  # Diversity before 2nd bin (included)
  div_2 <- sum(div_2)/sum(sub_div)
  #Find sigma
  sigma_prop <- sigma_prop <- seq(0.01, 1000, 0.01)
  trop_estimated <- function(sigma){
    return(pnorm(q = 19.48, mean = mu, sd = sigma))
  }
  CDF_emp <- unlist(lapply(X = sigma_prop,
                           FUN = trop_estimated))
  CDF_emp <- abs(CDF_emp - div_2)
  sigma_emp <- sigma_prop[which.min(CDF_emp)]
  return(sigma_emp)
}


# Simulate biodiversity gradients ---------------------------------------
# emp_data is a vector of length 6, containing empirical estimates of diversity for one time bin and one model
# n: Number of data points to generate (default 1000)
# type is the type of gradient to simulate: "unimodal", "bimodal", "flat" or "MIXED"
# if type="MIXED", the user has to specify the shape of the northern and southern gradients to simulate: "gaussian" or "flat"
simulate_biodiv_grad <- function(emp_data, n=1000, type, northern=NA, southern=NA) {
  # Sequence of latitudes used and hemisphere specification
  latitude <- seq(from = 90, to = -90, length.out = n)
  # Simulate a unimodal gradient using a gaussian distribution
  if (type == "unimodal") {
    sigma <- sigma_unimodal(data = emp_data, mu = 0)
    species_richness <- dnorm(latitude, mean = 0, sd = sigma)
  } 
  # Simulate a bimodal gradient using a combination of two Gaussian distributions with different means and sds extracted from the data
  if (type == "bimodal") {
    ## Northern hemisphere (cf. lat_bins dataframe)
    sigma_north <- sigma_bimodal(data = emp_data, loc = "NORTH")
    species_richness_north <- dnorm(latitude[which(latitude > 0)], mean = 45, sd = sigma_north)
    ## Southern hemisphere
    sigma_south <- sigma_bimodal(data = emp_data, loc = "SOUTH")
    species_richness_south <- dnorm(latitude[which(latitude <= 0)], mean = -45, sd = sigma_south)
    ## Combine
    species_richness <- append(species_richness_north,
                               species_richness_south)
  } 
  # Simulate a flat gradient
  if (type == "flat") {
    ## Set the species richness to a constant value, mean diversity
    species_richness <- rep(mean(emp_data, na.rm = TRUE), n)
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
  # Return simulated data
  return(species_richness)
}

## Example (gradients based on the 1st time stage and PALEOMAP) ----------
  # EMPLIRICAL DATA
div_raw <- readRDS("./data/processed/genus_counts.RDS")
div_bin <- div_raw$n_genera[which((div_raw$stage_bin == 1) &
                                    (div_raw$model == "PALEOMAP") &
                                    (is.na(div_raw$paleolat_bin) == FALSE))]
  # UNIMODAL
unimodal <- simulate_biodiv_grad(emp_data = div_bin,
                                 n = 1000,
                                 type = "unimodal")
scaling_cst <- max(div_bin) / max(unimodal)
unimodal <- unimodal * scaling_cst
  # BIMODAL
bimodal <- simulate_biodiv_grad(emp_data = div_bin,
                                n = 1000,
                                type = "bimodal")
scal_north <- max(div_bin[1:3]) / max(bimodal[1:(as.integer(length(bimodal))/2)])
scal_south <- max(div_bin[4:6]) / max(bimodal[(as.integer(length(bimodal))/2+1):length(bimodal)])
bimodal[1:(as.integer(length(bimodal))/2)] <- bimodal[1:(as.integer(length(bimodal))/2)] * scal_north
bimodal[(as.integer(length(bimodal))/2+1):length(bimodal)] <- bimodal[(as.integer(length(bimodal))/2+1):length(bimodal)] * scal_south
  # FLAT
flat <- simulate_biodiv_grad(emp_data = div_bin,
                             n = 1000,
                             type = "flat")
scaling <- mean(div_bin, na.rm = TRUE) / unique(flat)
flat <- flat * scaling
  # PLOT
lat <- seq(from = 90, to = -90, length.out = 1000)
lat_bins <- readRDS("./data/lat_bins.RDS")

par(mfrow = c(2,2))
plot(x = lat_bins$mid, y = div_bin, type = 'b', xlab = "", ylab = "")
title(main = "Actual diversity estimates", xlab = "Latitude", ylab = "Nb. genera")
plot(x = lat, y = unimodal, type = 'b', xlab = "", ylab = "")
title(main = "Simulated Unimodal", xlab = "Latitude", ylab = "Nb. genera")
plot(x = lat, y = bimodal, type = 'b', xlab = "", ylab = "")
title(main = "Simulated Bimodal", xlab = "Latitude", ylab = "Nb. genera")
plot(x = lat, y = flat, type = 'b', xlab = "", ylab = "")
title(main = "Simulated Flat", xlab = "Latitude", ylab = "Nb. genera")

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