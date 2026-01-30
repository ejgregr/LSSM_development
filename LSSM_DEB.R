#----------------------------------------------------------------------------
# Script:  LSSM_DEB.R
# Created: November 2025
## Purpose: Run a simple Dynamic Energy Budget model, written by ChatGPT
# Notes: Not a DEB, but a more sophisticated, 2-part logistic growth model
# Updated: Jan 2026. 
#   Deterministic (i.e., logistic growth) model not very satisfying. Updated
#   this with a 'DEB-like' model that uses Ordinary Differential Equations to 
#   predict growth based on light and temperature. 
#   Nutrients are assumed to be not limiting - a significant assumption.
#   Uses data from Weigel and Pfister (2021) to partition and parameterize 
#   carbon fixation, respiration, and DOC production.
#   Data from Pontier et al. (2023) used to attenuate growth as a function of 
#   temperature and light levels. 
#============================================================================


#----------------------------------
# Jan 27: EVERYTHING in DW please. 
# Define the nereo growth ODE system. We are growing wet weight using a single frond model.
# Includes measured growth parameters (Weigel and Pfister 2021) 
# and adjusted growth rates for temp and DLI responses (Pontier et al. 2024)
grow_kelp_model <- function(t, state, pars, env_data) {
  
  B_th <- state[1]  # g DW
  B_fr <- state[2]  # g DW
  
  with(as.list(pars), {
    
    # 1. Calculate ROW INDEX for looking up data.
    # Need to map simulation time (0, 1, 2) to row number (1, 2, 3)
    row_idx <- max(1, min(floor(t) + 1, nrow(env_data)))
    
    T_C      <- env_data$Temp[row_idx]
    day_h    <- env_data$day_hours[row_idx]
    night_h  <- env_data$night_hours[row_idx]
    real_DOY <- env_data$real_DOY[row_idx]
    
    # DLI - adjusted for reflectance
    DLI <- (1 - alpha_water) * env_data$DLI[row_idx]
    
    # --- Environmental Scaling (Pontier et al. 2024) ---
    fT <- fT_reparam(T_C, T_opt, sigma_warm)
    fL <- fL_reparam(DLI)
    
    # --- STIPE growth using an allometric (per Starko and Martone 2016)

    # target stipe biomass 
    B_th_target <- a_th_fr * (B_fr^b_th_fr)
    # Avoid division by zero at very small fronds
    B_th_target <- max(B_th_target, 1e-9)
    # Apply allometric relaxation to target 
    dB_th <- k_th * fT * max(0, B_th_target - B_th)
    
    # --- Mass-specific FROND fixation rates 
    # (g DW per hour) -> per g DW per day ---
    
    GP_day_umolC_gDW   <- max_fixation  * fT * fL * day_h
    GP_night_umolC_gDW <- dark_fixation * fT * night_h
    GP_umolC_gDW       <- GP_day_umolC_gDW + GP_night_umolC_gDW
    
    # And now translate to the Carbon fixed this day (for DOC leakage, below)
    GP_day_umolC   <- B_fr * GP_day_umolC_gDW
    GP_night_umolC <- B_fr * GP_night_umolC_gDW
    GP_umolC       <- GP_day_umolC + GP_night_umolC
    
    # --- DOC Leakage - variable as a function of light (W&P)
    # Define a daily fraction
    DOC_day_frac <- DOC_day_min + (DOC_day_max - DOC_day_min) /
      (1 + exp(DOC_DLI_k * (DLI - DOC_DLI50)))
    # Bound to [0,1]
    DOC_day_frac <- max(0, min(1, DOC_day_frac))
    
    # Calculate the DOC
    DOC_day_umolC   <- DOC_day_frac   * GP_day_umolC
    DOC_night_umolC <- DOC_night_frac * GP_night_umolC
    DOC_umolC       <- DOC_day_umolC + DOC_night_umolC
    
    #--- Calculate respiration proportion
    R_prod_umolC  <- respired_prod * GP_day_umolC # proportion of daily C fixing
    R_maint_umolC <- R_maint_umolC_gDW_h * fT * 24 * (B_fr + B_th) # Maintenance respiration
    
    # Net µmol fixed for this simulated day
    NetC_umolC <- GP_umolC - DOC_umolC - R_prod_umolC - R_maint_umolC  
    
    
    # --- Senescence - variable sloughing Logic (Sigmoid)
    # Create a smooth transition from low summer sloughing to high fall sloughing
    # Steepness is set to 0.1 (hardcoded) for a gradual 30-day shift
    slough_today <- slough_min + (slough_max - slough_min) / 
      (1 + exp(-senescence_k * (real_DOY - senescence_day)))
    
    # Turn off sloughing  (debugging)
    # slough_today <- 0
    
    # Convert frond sloughing (g DW/day) to µmol C/day for C-budget plotting
    Slough_gDW <- slough_today * B_fr
    Slough_umolC <- Slough_gDW * gDW_to_gC / 12.011 * 1e6
    
    # Convert total net C -> total biomass gain (DW)
    # µmol C -> mol C -> g C -> g DW
    tissue_growth <- NetC_umolC * 1e-6 * 12.011 / gDW_to_gC
    
    # --- Net Growth (Biomass Change) ---
    dB_fr <- tissue_growth - (slough_today * B_fr)
    
    # as.numeric() ensures name stays as defined. ode() oddity.
    return(list(
      c(dB_th, dB_fr),
      real_DOY = real_DOY,
      fL = fL, 
      fT = fT,
      DLI = DLI,
      GP_day_umolC    = as.numeric(GP_day_umolC),     # Carbon fixed during day
      GP_night_umolC  = as.numeric(GP_night_umolC),   # Carbon fixed during night
      GP_umolC        = as.numeric(GP_umolC),         # Total carbon fixed 
      DOC_umolC       = as.numeric(DOC_umolC),        # Total DOC loss
      DOC_day_umolC   = as.numeric(DOC_day_umolC),    # Day time DOC loss
      DOC_night_umolC = as.numeric(DOC_night_umolC),  # Night time DOC loss
      Slough_umolC    = as.numeric(Slough_umolC),     # Total slough loss
      R_prod_umolC    = as.numeric(R_prod_umolC),     # Respiration from production
      R_maint_umolC   = as.numeric(R_maint_umolC)     # Respiration from maintenance
    ))
  })
}


model_params <- list(
  # Frond dynamics
  max_fixation   = 120.78, # Max observed umol C/gDW/h in June - Weigel and Pfister
  dark_fixation  = 3.75,   # Light indepdt nighttime fixation (C/gDW/h) - Weigel and Pfister
  respired_prod  = 0.15,   # 
  R_maint_umolC_gDW_h = 0.05, # set/fit later
  slough_min     = 0.001,   # Min slough rate for blades (early)
  slough_max     = 0.01,   # Max slough rate for blades (late)
  senescence_day = 240,    # The day (approx late August) when loss ramps up
  senescence_k   = 0.025,   # Rate of senescence increase. .05 ~ a 60 day ramp up

    # DOC leakage as fraction of GP during DAY: decreases with light (DLI)
  # Based on Weigel and Pfister
  DOC_day_max = 0.24,   # low-light PER (e.g., ~0.23–0.24)
  DOC_day_min = 0.08,   # high-light PER (e.g., ~0.08)
  DOC_DLI50   = 43,     # Calc'd from DLI distributional data
  DOC_DLI_k   = 0.15,   # Calc'd from DLI distributional data
  #  DOC_DLI50   = 30,     # DLI at midpoint of decline (mol m^-2 d^-1)
  #  DOC_DLI_k   = 0.25,   # steepness (per DLI unit)
  
  # Night DOC leakage (simple assumption): fraction of night GP
  DOC_night_frac = 0.10, # tune later; keep small under non–nutrient-limited assumption  
  
  # Temperature restriction based on Pontier et al. 
  # To span Pontier’s stated range: sigma_warm ~ 0.75 (steep penalty),  1.7 (mild  penalty)
  T_opt      = 10,
  sigma_warm = 1.5,
  
  # Stipe dynamics
  a_th_fr = 0.5,   # coefficient a. 
  b_th_fr = 1.0,   # exponent b. 
  k_th    = 0.3,   # enough to prevent chronic stipe lag, but not instantaneous
  
  # Constants
  gDW_to_gC    = 0.313, # Estimated. RSSM says 0.313. 
  wet_to_dry   = 0.13,  # From RSSM.
  alpha_water  = 0.06  # Fraction of DLI reflected at water surface
)

# 3. Define initial state (start with small kelp plant)
# Starting biomass based on young sporophyte with established first blade = 1 to 4 g WW total.
initial_state <- c(
  B_th = 1 * model_params$wet_to_dry,  # Stated value for stipe mass in (g WW)
  B_fr = 2 * model_params$wet_to_dry   # Stated value for frond mass in (g WW)
)

# 4. Define time sequence grow_kelp(from day 1 to day 100)
times <- env_daily$day

# 5. Solve the ODE system using the 'ode' function
output <- ode(
  y = initial_state,
  times = times,
  func = grow_kelp_model,
  parms = model_params,
  env_data = env_daily,
  method = "rk4" # Runge-Kutta 4th order method
)

out_df <- as.data.frame( output )
out_df <- cbind( out_df, "B_th_WW" = out_df$B_th/model_params$wet_to_dry,
                         "B_fr_WW" = out_df$B_fr/model_params$wet_to_dry)



### Fin.














