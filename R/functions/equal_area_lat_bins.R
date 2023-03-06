# Calculate surface area of polygon (latitudinally bounded)
# x: the minimum and maximum latitude of the segment
# crs: desired CRS to be used
bounded_surface_area <- function(x, crs = 4326) {
  # Define lngs
  lng <- c(0, 0, 90, 90)
  # Define lats
  lat <- c(x[1], x[2], x[1], x[2])
  # Create point dataframe
  p <- data.frame(lng = lng,
                  lat = lat)
  # Create polygon for area calc
  poly <- p |>
    sf::st_as_sf(coords = c("lng", "lat"), crs = crs) |>
    sf::st_bbox() |>
    sf::st_as_sfc()
  # Calculate area
  as.numeric(sf::st_area(x = poly) * 4)
}
# Example
# bounded_surface_area(x = c(80, 90), crs = 4326)

# Generate equal area latitudinal bins
# n: number of bins desired
# fun: function be used to calculate geographic area
# crs: desired CRS to be used
# by: the step size used for calculating latitudinal bins (smaller values)
# result in more evenness between bins.
equal_area_lat_bins <- function(n = 12, fun = bounded_surface_area,
                                crs = 4326, by = 0.1) {
  # Calculate surface area of the Earth assuming an oblate spheroid
  earth_area <- bounded_surface_area(x = c(-90, 90))
  # How many bins should be made?
  bin_area <- earth_area / n
  # Specify number of steps (higher number of steps, increased equality)
  vals <- seq(0, 90, by = by)
  # Set starting index
  indx <- 1
  # Half number of bins to avoid recalculating hemispheres
  n <- n / 2
  # Create empty df for populating
  bins <- data.frame(bin = 1:n, max = rep(NA, n), mid = rep(NA, n),
                     min = rep(NA, n), area = rep(NA, n))
  # Run across desired number of bins
  for (i in 1:n) {
    # Set area value to 0
    a <- 0
    # Record starting index
    start_indx <- indx
    while (a < bin_area) {
      # Increase index
      indx <- indx + 1
      # Extract latitudes
      lat <- c(vals[indx-1], vals[indx])
      # If NAs are present, index has gone beyond vals length (90 lat)
      if (indx > length(vals)) {
        indx <- indx - 1
        break
      }
      # Calculate area
      a <- a + fun(x = lat)
    }
    # Record final index
    end_indx <- indx
    # Assign breaks
    bins$max[i] <- vals[end_indx]
    bins$min[i] <- vals[start_indx]
    # Assign area value
    bins$area[i] <- a
  }
  # Create southern hemisphere df
  s_bins <- bins
  bins <- bins[order(bins$bin, decreasing = TRUE), ]
  s_bins$min <- s_bins$min * -1
  s_bins$max <- s_bins$max * -1
  # Reverse column values
  s_bins[, c("max", "min")] <- s_bins[, c("min", "max")]
  # Bind data
  bins <- rbind.data.frame(bins, s_bins)
  # Update bin numbers
  bins$bin <- 1:nrow(bins)
  # Drop row names
  row.names(bins) <- NULL
  # Add area proportions
  bins$area_prop <- bins$area / sum(bins$area)
  # Add mid bin
  bins$mid <- (bins$max + bins$min) / 2
  # Return bins
  return(bins)
}
# Example
# bins <- equal_area_lat_bins(by = 0.01)
