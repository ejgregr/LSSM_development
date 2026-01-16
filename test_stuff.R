library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. GEOMETRIC ASSIGNMENT FUNCTION ---
assign_geometric_start_days <- function(n, interval) {
  start_days <- numeric(n)
  count_assigned <- 0
  generation <- 1
  while (count_assigned < n) {
    if (generation == 1) { wave_size <- 2 } else { wave_size <- 2^(generation - 1) }
    n_in_this_wave <- min(wave_size, n - count_assigned)
    current_day <- 1 + ((generation - 1) * interval)
    start_idx <- count_assigned + 1
    end_idx   <- count_assigned + n_in_this_wave
    start_days[start_idx:end_idx] <- current_day
    count_assigned <- count_assigned + n_in_this_wave
    generation <- generation + 1
  }
  return(start_days)
}

# --- 2. PARAMETER SETUP ---
# A. Define Blade Count FIRST
n_blades <- round(rnorm(1, mean = 41.5, sd = 13.2))

# B. Create Population DataFrame
blade_population <- data.frame(
  max_fix   = rnorm(n_blades, 78.5, sd = 33.5),
  doc_day   = rnorm(n_blades, 10.8, sd = 4.51),
  doc_night = rnorm(n_blades, 3.65, sd = 3.84),
  resp_prop = rnorm(n_blades, 0.15, sd = 0.08) 
)

# C. Apply Geometric Timing (10 Day Interval)
split_interval <- 10 
blade_population$start_day <- assign_geometric_start_days(n_blades, split_interval)

# D. Add Jitter (Spread the waves slightly)
jitter <- round(rnorm(n_blades, mean = 0, sd = 2)) 
blade_population$start_day <- pmax(1, blade_population$start_day + jitter)

# E. Biological "Guard Rails" 
# CRITICAL CHANGE: Raised Min Fixation to 40 to prevent "Dud" blades that don't grow
blade_population$max_fix   <- pmax(blade_population$max_fix, 40) 
blade_population$resp_prop <- pmin(pmax(blade_population$resp_prop, 0.05), 0.50)
blade_population$doc_day   <- pmax(blade_population$doc_day, 0)
blade_population$doc_night <- pmax(blade_population$doc_night, 0)

# Sort for cleaner plotting (Oldest blades first)
blade_population <- blade_population[order(blade_population$start_day), ]

# --- 3. ODE FUNCTION (WITH STAGGERED MASK) ---
grow_kelp_staggered <- function(t, state, pars, env_data, blade_pars) {
  B_th        <- state[1]
  B_fr_vector <- state[-1]
  
  with(as.list(pars), {
    # Time mapping
    row_idx  <- max(1, min(floor(t) + 1, nrow(env_data)))
    real_DOY <- env_data$day[row_idx]
    
    # Env Data
    T_val   <- env_data$Temp[row_idx]
    I_val   <- env_data$DLI[row_idx]
    day_h   <- env_data$day_hours[row_idx]
    night_h <- env_data$night_hours[row_idx]
    
    # Scalars
    fT <- ifelse(T_val <= 10.15, exp(-0.05*(T_val-10.15)^2), exp(-0.64*(T_val-10.15)))
    fL_blade <- ifelse(I_val < 20, I_val/20, ifelse(I_val > 40, max(0, 1-(I_val-40)/20), 1.0))
    fT_resp  <- 2.0^((T_val - 10) / 10)
    
    # Variable Sloughing
    current_slough <- slough_min + (slough_max - slough_min) / 
      (1 + exp(-0.1 * (real_DOY - day_senescence)))
    
    # --- ACTIVE MASK ---
    # Calculates 0 or 1 for every blade based on current time 't'
    active_mask <- ifelse(t >= blade_pars$start_day, 1, 0)
    
    # Stipe Growth
    stipe_mult <- max(0, (1 - (B_th / K_th)^2))
    dB_th <- (r_max_th * fT * stipe_mult * B_th) - (slough_rate_th * B_th)
    
    # Blade Growth
    GP_day   <- (blade_pars$max_fix * fT * fL_blade) * day_h
    GP_night <- (dark_fixation * fT) * night_h 
    GP_total <- GP_day + GP_night
    DOC_loss  <- (blade_pars$doc_day * fT * fL_blade * day_h) + (blade_pars$doc_night * fT * night_h)
    Resp_loss <- GP_total * blade_pars$resp_prop * fT_resp
    Net_C_Gain <- GP_total - DOC_loss - Resp_loss
    
    tissue_growth <- (Net_C_Gain * 1e-6 * 12.011) / dry_to_C
    
    # Apply Mask: Inactive blades have 0 growth and 0 sloughing
    dB_fr_vector <- (tissue_growth * B_fr_vector * active_mask) - 
      (current_slough * B_fr_vector * active_mask)
    
    return(list(c(dB_th, dB_fr_vector)))
  })
}

# --- 4. MODEL CONFIG & RUN ---
model_params <- list(
  r_max_th = 0.13, dark_fixation = 3.75, DOC_leakage = 0.162, 
  resp_rate = 0.10, wet_to_dry = 0.13, dry_to_C = 0.15, 
  K_th = 2.0, slough_rate_th = 0.001,
  slough_min = 0.002, slough_max = 0.10, day_senescence = 240
)

# Initial State
init_B_fr <- rep(0.1, n_blades)
initial_state <- c(B_th = 0.05, init_B_fr)

# Run ODE
output <- ode(
  y = initial_state, 
  times = 0:(nrow(env_daily)-1), # Ensure this aligns with your env_data rows
  func = grow_kelp_staggered, 
  parms = model_params,
  env_data = env_daily,
  blade_pars = blade_population
)

# --- 5. PLOT INDIVIDUAL TISSUES (LOG SCALE) ---
plot_individual_tissues(output, env_daily)