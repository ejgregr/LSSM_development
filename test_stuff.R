# Generate timestamps for 24 hours
start_time <- as.POSIXct("2024-03-01 00:00:00")  # Start time
end_time <-   as.POSIXct("2024-03-31 00:00:00")    # End time (24 hours later)
timestamps <- seq(from = start_time, to = end_time, by = "hour")

# Shift the sinusoidal pattern 6 hours forward
# We add a phase shift of pi/2 to the sinusoidal function to shift the peak 6 hours ahead
photon_concentration <- 100 + 50 * sin(2 * pi * ((as.numeric(format(timestamps, "%H")) / 24) - 6 / 24))

# Create the time series data frame
photon_time_series <- data.frame(Time = timestamps, Photon_Concentration = photon_concentration)

# View the first few rows of the shifted time series
plot(photon_time_series, type='l')
