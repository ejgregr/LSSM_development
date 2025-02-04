


#---- OLDER Stuff ... 

# Data structure will track a plant over 12 months. 
# All other results will be scaled to/from this. 

# An initial plant - baby sporophyte
iPlant <- dataframe( "frond_l" = 1.0,   # cm
                     "frond_w" = 0.05,  # g
                     "frond_n" = 50     # number of fronds
                     "stype_d" = 0.1)   # stype diameter for mass scaling


# Known parameters
B0 <- 0.01    # Initial biomass (kg)
K <- 10       # Carrying capacity (kg)
B_final <- 10 # Final biomass after t months
t <- 6        # Time (months)

# Logistic growth function
logistic_growth <- function(B0, K, t, r) {
  B_predicted <- K / (1 + ((K - B0) / B0) * exp(-r * t))
  return(B_predicted)
}

# Monthly growth rate ESTIMATE for 6 months - optimal
x <- NULL
for (t in 1:6) {
y <- logistic_growth( 0.01, 10, t, 1.565)
x <- c( x, y)
}
plot(x)




# Parameters
B0 <- 0.01    # Initial biomass (kg)
K <- 10       # Carrying capacity (kg)
t <- 1:6      # Time steps (months)
T <- c(10, 12, 15, 18, 20, 22)  # Example temperatures (°C)

# Environmental influence on growth rate
r_max <- 0.8  # Maximum growth rate
T_opt <- 20   # Optimal temperature (°C)
T_range <- 10 # Temperature range for growth (°C)


# Temperature-dependent growth rate function
r_temp <- function(T) {
  r_max * (1 - ((T - T_opt)^2) / (T_range^2))
}


# Logistic growth model with temperature influence
logistic_growth_temp <- function(t, T) {
  r_T <- r_temp(T)  # Calculate growth rate for each temperature
  B_t <- K / (1 + ((K - B0) / B0) * exp(-r_T * t))
  return(B_t)
}

# Apply model over time and temperatures
biomass <- sapply(T, function(temp) logistic_growth_temp(t, temp))


# Plot results
matplot(t, biomass, type = "l", lty = 1, col = rainbow(length(T)),
        xlab = "Time (months)", ylab = "Biomass (kg)",
        main = "Temperature-Dependent Logistic Growth")
legend("bottomright", legend = paste("T =", T, "°C"), col = rainbow(length(T)), lty = 1)


# And then we  scale back up to aPatch using our own allometric scaling ... whew. 



From Pontier et al. 2024:
Measured growth at 3 sites across 4 years
  
Mean growth rates per sampling events across all sites and years:
Blades : 0.83 ± 0.20 cm/day to 8.25 ± 0.78 cm/day
stipes : 0.23 ± 0.15 cm/day to 9.31 ± 0.75 cm/day


aPatch



