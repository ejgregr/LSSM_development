#----------------------------------------------------------------------------
# Script:  LSSM_ODE.R
# Created: January 2026
## Purpose: Run a simple Dynamic Energy Budget-like model based on coupled, 
# ordinary differential equations. Written with support from Gemini.
# Updates:
#
#============================================================================

#---- Parameterize and run the multiple blade ODE model -----

# 1. Pick number of blades for individual plant
n_blades <- round(rnorm(1, mean = 41.5, sd = 13.2))


# 2. Create blade population data frame. 
#   Draw unique metabolic parameters for EACH blade based on W&P
#   Reported SEs have been converted to SDs based on sample size (see report)
blade_population <- data.frame(
  max_fix   = rnorm(n_blades, 78.5, sd = 33.5),
  doc_day   = rnorm(n_blades, 10.8, sd = 4.51),
  doc_night = rnorm(n_blades, 3.65, sd = 3.84),
  resp_prop = rnorm(n_blades, 0.15, sd = 0.08)  # SE retained for stability
)

# 3. Add start days to the population
# TO DO: Explore this interval to ensure all blades done by day 60
split_interval <- 10  # Days between doubling events (approx 2 weeks)
# Assign the start_days to the blades
blade_population$start_day <- assign_geometric_start_days(n_blades, split_interval)
# Sort for cleaner plotting (Oldest blades first)
blade_population <- blade_population[order(blade_population$start_day), ]

# Add a tiny bit of noise (+/- 2 days) so they don't appear strictly robotically
jitter <- round(runif(n_blades, 0, 2))
blade_population$start_day <- pmax(1, blade_population$start_day + jitter)

# Check the distribution
table(blade_population$start_day)

# 3. Biological "Guard Rails" (Truncation)
# Prevent negative physiology or impossible values
blade_population$max_fix   <- pmax(blade_population$max_fix, 40)  # Min fixation 40 umol
blade_population$doc_day   <- pmax(blade_population$doc_day, 0)   # No negative leakage
blade_population$doc_night <- pmax(blade_population$doc_night, 0) # No negative leakage

# Clamp respiration: 
# Min 0.05 (super efficient) to Max 0.50 (dying tissue)
blade_population$resp_prop <- pmin(pmax(blade_population$resp_prop, 0.05), 0.50)

model_params <- list(
  r_max_th       = 0.13,   # Max daily rate of stipe growth (in cm/day)
  dark_fixation  = 3.75,   # Light indepdt nighttime fixation - Weigel and Pfister
  DOC_leakage    = 0.162,  # Ave of daily fixed carbon lost as DOC - W & P
  resp_rate      = 0.10,   # Estimated 
  wet_to_dry     = 0.13,   # portion of kelp that is not water 
  dry_to_C       = 0.15,   # portion of dry kelp that is carbon
  K_th           = 2.0,    # The biomass (g DW) at with stipe stops growing
  slough_min     = 0.002,  # Summer baseline (0.2% loss/day)
  slough_max     = 0.10,   # Senescence peak (10% loss/day)
  day_senescence = 240     # The day (approx late August) when loss ramps up
)

# 3. Initial state
init_B_fr <- rep(0.1, n_blades) # Start with 0.01g DW per blade
initial_state <- c(B_th = 0.05, init_B_fr)

# 4. Define time sequence (from day 1 to day 100)
times <- env_daily$day

# d. Run the ODE
output <- ode(
  y = initial_state, 
  times = times, 
  func = grow_kelp_multiple_blades, 
  parms = model_params,
  env_data = env_daily,
  blade_pars = blade_population
)




#---- Parameterize and run the single blade ODE model -----
# # Define parameters (initial estimates, tune these with the paper data)
# model_params <- list(
#   r_max_th       = 0.13,    # Max daily rate of stipe growth (in cm/day)
#   max_fixation   = 120.78,  # Peak umol C/g/h in June - Weigel and Pfister
#   dark_fixation  = 3.75,    # Light indepdt nighttime fixation - Weigel and Pfister
#   DOC_leakage    = 0.162,   # Ave of daily fixed carbon lost as DOC - W & P
#   resp_rate      = 0.10,    # Estimated 
#   carbon_fract   = 0.15,    # Estimated
#   K_th           = 2.0,     # The biomass (g DW) at with stipe stops growing
#   slough_rate_th = 0.001,
#   slough_rate_fr = 0.01 
# )
# 
# # 3. Define initial state (start with small kelp plant)
# initial_state <- c(
#   B_th = 0.005,  # initial stipe dry mass (g DW Carbon)
#   B_fr = 0.5     # initial frond dry mass (g DW Carbon)
# )
# 
# # 4. Define time sequence (from day 1 to day 100)
# times <- env_daily$day
# 
# # 5. Solve the ODE system using the 'ode' function
# output <- ode(
#   times = times,
#   y = initial_state,
#   func = grow_kelp_model,
#   parms = model_params,
#   env_data = env_daily,
#   method = "rk4" # Runge-Kutta 4th order method
# )
