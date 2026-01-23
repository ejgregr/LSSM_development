#----------------------------------------------------------------------------
# Script:  LSSM_ODE.R
# Created: January 2026
## Purpose: Run a simple Dynamic Energy Budget-like model based on coupled, 
# ordinary differential equations. Written with support from Gemini.
# Updates:
# 2026/01/16: After several days of creating an increasingly complex kelp ODE
#   including multiple fronds with different emergence times, decided to fall 
#   back to a single, super-individual because the 'splitting' of the fronds
#   into two blades made conserving the mass difficult in an ODE context. 
#   Apparently I rediscovered something called the mass deletion paradox. :\
#   Oh well. Little to be gained with this approach except that it was cool. :\
##
#============================================================================

#---- Parameterize and run the single blade ODE model -----
# Define parameters (initial estimates, tune these with the paper data)
model_params <- list(
  r_max_th       = 0.13,    # Max daily rate of stipe growth (assume proportional to Pontier's cm/day)
#  max_fixation   = 120.78,  # Max observed umol C/g/h in June - Weigel and Pfister
  max_fixation   =    78.5,  # Mean observed umol C/g/h in June - Weigel and Pfister
  dark_fixation  = 3.75,    # Light indepdt nighttime fixation (C/g/h) - Weigel and Pfister
  DOC_leakage    = 0.162,   # Ave of daily fixed carbon lost as DOC - W & P
  
  #resp_rate      = 0.15,    # Estimated
  resp_frac_GP    = 0.10, # Fraction of GP respired during c-fixing (W&P suggest resp ~10% of gross fixation)
  maint_umolC_gDW = 10,   # Maintenance respiration rate (umol C / g DW / day) at Tref
  T_ref          = 10,    # reference temperature for Q10
  Q10_resp       = 2.0,   # Q10 for respiration
  maint_size_exp = 1.0,   # exponent for biomass scaling (1 = linear)
  maint_seas_amp = 0.0,   # optional seasonal increase (0 disables)

  gDW_to_gC      = 0.313,    # Estimated. RSSM says 0.313. 
  wet_to_dry     = 0.13,    # From RSSM.
  K_th           = 2.0,     # The biomass (g WW) at which stipe stops growing
  slough_min     = 0.01,   # Min slough rate for blades (early)
  slough_max     = 0.10,     # Max slough rate for blades (late)
  senescence_day = 240     # The day (approx late August) when loss ramps up
)

# 3. Define initial state (start with small kelp plant)
initial_state <- c(
  B_th = 0.50,  # initial stipe mass (g WW)
  B_fr = 0.75   # initial frond mass (g WW)
)

# 4. Define time sequence (from day 1 to day 100)
times <- env_daily$day

# 5. Solve the ODE system using the 'ode' function
output <- ode(
  times = times,
  y = initial_state,
  func = grow_kelp_model,
  parms = model_params,
  env_data = env_daily,
  method = "rk4" # Runge-Kutta 4th order method
)


env_daily[1:20,]



#---- Nereo growth parameters (older) ----
# B_init      <- 0.025  # Initial mass of sporophyte (est. at 25 mg based on ChatGPT)
# B_max       <- 9.23   # Max biomass for an adult Nereo. Calculated from field results (Weigel and Pfister 2021)
# r_max       <- 0.065  # Maximum daily growth rate, assuming 6 month growth and mature plant is 9.23 kg 
# sp_start    <- 0.005  # Established sporophyte mass (5 g)
# wet_to_dry  <- 0.13   # Water content of Nereo (Bullen et al. 2024)
# dry_to_C    <- 0.25   # Carbon content of Nereo (dry) (Bullen et al. 2024)
# 
# #---- Environmental influence on growth rate
# T_opt     <- 10        # Optimal temperature (°C)
# T_max     <- 14        # Maximal temperature for growth (°C)
# DLI_opt   <- 30        # Optimal DLI in mol/m2/day
# DLI_range <- 20        # Range (+/-) that allows growth in mol/m2/day 



#---- Parameterize and run the multiple blade ODE model -----

# 1. Pick number of blades for individual plant
# n_blades <- round(rnorm(1, mean = 41.5, sd = 13.2))
# 
# 
# # 2. Create blade population data frame. 
# #   Draw unique metabolic parameters for EACH blade based on W&P
# #   Reported SEs have been converted to SDs based on sample size (see report)
# blade_population <- data.frame(
#   max_fix   = rnorm(n_blades, 78.5, sd = 33.5),
#   doc_day   = rnorm(n_blades, 10.8, sd = 4.51),
#   doc_night = rnorm(n_blades, 3.65, sd = 3.84),
#   resp_prop = rnorm(n_blades, 0.15, sd = 0.08)  # SE retained for stability
# )
# 
# # 3. Add start days to the population
# # TO DO: Explore this interval to ensure all blades done by day 60
# split_interval <- 10  # Days between doubling events (approx 2 weeks)
# # Assign the start_days to the blades
# blade_population$start_day <- assign_geometric_start_days(n_blades, split_interval)
# # Sort for cleaner plotting (Oldest blades first)
# blade_population <- blade_population[order(blade_population$start_day), ]
# 
# # Add a tiny bit of noise (+/- 2 days) so they don't appear strictly robotically
# jitter <- round(runif(n_blades, 0, 2))
# blade_population$start_day <- pmax(1, blade_population$start_day + jitter)
# 
# # Check the distribution
# #table(blade_population$start_day)
# #blade_population$start_day <- 1
# 
# # 3. Biological "Guard Rails" (Truncation)
# # Prevent negative physiology or impossible values
# blade_population$max_fix   <- pmax(blade_population$max_fix, 40)  # Min fixation 40 umol
# blade_population$doc_day   <- pmax(blade_population$doc_day, 0)   # No negative leakage
# blade_population$doc_night <- pmax(blade_population$doc_night, 0) # No negative leakage
# 
# # Clamp respiration: 
# # Min 0.05 (super efficient) to Max 0.50 (dying tissue)
# blade_population$resp_prop <- pmin(pmax(blade_population$resp_prop, 0.05), 0.50)
# 
# model_params <- list(
#   r_max_th       = 0.13,   # Max daily rate of stipe growth (in cm/day)
#   dark_fixation  = 3.75,   # Light indepdt nighttime fixation - Weigel and Pfister
#   DOC_leakage    = 0.162,  # Ave of daily fixed carbon lost as DOC - W & P
#   resp_rate      = 0.10,   # Estimated 
#   wet_to_dry     = 0.13,   # portion of kelp that is not water 
#   dry_to_C       = 0.25,   # portion of dry kelp that is carbon
#   K_th           = 2.0,    # The biomass (g DW) at with stipe stops growing
#   slough_min     = 0.002,  # Summer baseline (0.2% loss/day)
#   slough_max     = 0.10,   # Senescence peak (10% loss/day)
#   day_senescence = 240     # The day (approx late August) when loss ramps up
# )
# 
# # 3. Initial state
# init_B_fr <- rep(1.0, n_blades) # Start with 0.01g DW per blade
# initial_state <- c(B_th = 0.5, init_B_fr)
# 
# # 4. Define time sequence (from day 1 to day 100)
# times <- env_daily$day
# 
# # d. Run the ODE
# output <- ode(
#   y = initial_state, 
#   times = times, 
#   func = grow_kelp_multiple_blades, 
#   parms = model_params,
#   env_data = env_daily,
#   blade_pars = blade_population
# )




