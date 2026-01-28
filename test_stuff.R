model_params <- list(
  r_max_th       = 0.13,  # Max daily rate of stipe growth (a proxy for proportion increase per day, based on Pontier's cm/day)
  max_fixation   = 120.78,# Max observed umol C/gDW/h in June - Weigel and Pfister
  dark_fixation  = 3.75,  # Light indepdt nighttime fixation (C/gDW/h) - Weigel and Pfister
  DOC_leakage_d  = 0.162, # Ave of daily fixed carbon lost as DOC - W & P
  respired_prod  = 0.15,   # was respired_d
  R_maint_umolC_gDW_h = 0, # set/fit later
  gDW_to_gC      = 0.313, # Estimated. RSSM says 0.313. 
  wet_to_dry     = 0.13,  # From RSSM.
  K_th           = 5,     # The biomass (g WW) at which stipe stops growing
  slough_min     = 0.005,  # Min slough rate for blades (early)
  slough_max     = 0.05,  # Max slough rate for blades (late)
  senescence_day = 240    # The day (approx late August) when loss ramps up
)

# 3. Define initial state (start with small kelp plant)
# Starting biomass based on young sporophyte with established first blade = 1 to 4 g WW total.
initial_state <- c(
  B_th = 1.0 * model_params$wet_to_dry,  # Stated value for stipe mass in (g WW)
  B_fr = 1.5 * model_params$wet_to_dry   # Stated value for frond mass in (g WW)
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

#plot_kelp_biomass_DW(out_df, log=T ) 
plot_kelp_biomass_WW(out_df, log=F ) 

head(out_df)

# TO DO: CHeck the temp curve - Starts funny. 
plot( out_df$B_fr )
plot( out_df$fL )
plot( out_df$slough_today )


range( env_daily$ DLI )


plot_fT_curve()
plot_fL_curve()

plot_temperature_scaling( env_daily )
plot_DLI_scaling( env_daily )



#----------------------------------
# Jan 27: EVERYTHING in DW please. 
grow_kelp_model <- function(t, state, pars, env_data) {
  
  B_th <- state[1]  # g DW
  B_fr <- state[2]  # g DW
  
  with(as.list(pars), {
    
    # 1. Calculate ROW INDEX for looking up data.
    # Need to map simulation time (0, 1, 2) to row number (1, 2, 3)
    row_idx <- max(1, min(floor(t) + 1, nrow(env_data)))
    
    T_C      <- env_data$Temp[row_idx]
    DLI      <- env_data$DLI[row_idx]
    day_h    <- env_data$day_hours[row_idx]
    night_h  <- env_data$night_hours[row_idx]
    real_DOY <- env_data$real_DOY[row_idx]
    
    # --- Environmental Scaling (Pontier et al. 2024) ---
    # Temperature scaling applies to both stipe and blade, Light only to blade. 
    # fT <- ifelse(T_C <= 10.15,
    #              exp(-0.05 * (T_C - 10.15)^2),
    #              exp(-0.64 * (T_C - 10.15)))
    # 
    # fL <- ifelse(DLI < 20, DLI/20,
    #              ifelse(DLI > 40, max(0, 1 - (DLI - 40)/20), 1.0))
    # 
    # # Restrict degree of environmental scaling (debugging)
    # fL <- max( fL, 0.5)
    # fT <- max( fT, 0.5)
    
    fT <- fT_reparam(T_C)
    fL <- fL_reparam(DLI)
    
    # --- STIPE growth. r_max_th is based on observed growth (Pontier et al. 2024)
    # K_th is in g ww so convert. 
    stipe_mult <- max(0, 1 - ( B_th / (K_th * wet_to_dry) )^2)
    dB_th      <- r_max_th * fT * stipe_mult * B_th  # g DW / day
    
    # --- Mass-specific FROND fixation rates 
    # (g DW per hour) -> per g DW per day ---
    
    GP_night_umolC_gDW <- dark_fixation * fT * night_h
    GP_day_umolC_gDW   <- max_fixation  * fT * fL * day_h
    GP_umolC_gDW       <- GP_day_umolC_gDW + GP_night_umolC_gDW
    
    # --- Senescence - variable sloughing Logic (Sigmoid) ---
    # Create a smooth transition from low summer sloughing to high fall sloughing
    # Steepness is set to 0.1 (hardcoded) for a gradual 30-day shift
    slough_today <- slough_min + (slough_max - slough_min) / 
      (1 + exp(-0.1 * (real_DOY - senescence_day)))
    
    # Turn off sloughing  (debugging)
    #slough_today <- 0
    
    # Convert to TOTAL fluxes for fronds (photosynthetic tissue)
    GP_umolC  <- B_fr          * GP_umolC_gDW
    DOC_umolC <- DOC_leakage_d * GP_umolC # proportion fixed lost as DOC
    R_prod_umolC  <- respired_prod * GP_umolC # proportion fixed respired
    R_maint_umolC <- R_maint_umolC_gDW_h * fT * 24 * (B_fr + B_th) # Maintenance respiration
    
    # Net µmol fixed for this simulated day
    NetC_umolC <- GP_umolC - DOC_umolC - R_prod_umolC - R_maint_umolC  

    # Convert total net C -> total biomass gain (DW)
    # µmol C -> mol C -> g C -> g DW
    tissue_growth <- NetC_umolC * 1e-6 * 12.011 / gDW_to_gC
    
    
    # --- Net Growth (Biomass Change) ---
    dB_fr <- tissue_growth - (slough_today * B_fr)
    
    
    return(list(
      c(dB_th, dB_fr),
      GP_umolC_gDW = GP_umolC_gDW,
      NetC_umolC   = NetC_umolC,
      fL = fL, 
      fT = fT,
      slough_today = slough_today,
      R_prod_umolC  = R_prod_umolC,
      R_maint_umolC = R_maint_umolC
    ))
  })
}



#---- Plotting functions ----

# Takes temp from env_daily df, calculates and plots associated scaling factor
plot_temperature_scaling <- function( env_dat ){

  # Plotting bits 
  DOY     <- env_dat$real_DOY
  T_daily <- env_dat$Temp
  fT_daily <- fT_reparam( T_daily )

  # Set up plot margins to allow second y-axis
  par(mar = c(5, 4, 4, 4) + 0.1)

  # --- Left axis: Temperature ---
  plot( 
    DOY, T_daily,
    type = "l", col = "steelblue", lwd = 2,
    xlab = "Day of Year", ylab = "Temperature (°C)", bty = "l" )
  
  # --- Right axis: fT ---
  par(new = TRUE)
  
  plot(
    DOY, fT_daily,
    type = "l", col = "firebrick", lwd = 2,
    axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1) )
  
  axis(
    side=4, at=seq(0, 1, by=0.2), col = "firebrick", col.axis = "firebrick" )
  
  mtext(
    "Temperature scaling factor (fT)",
    side = 4,
    line = 3,
    col = "firebrick"
  )
  
  # --- Legend ---
#  legend(
#    "topright",
#    legend = c("Temperature", "Temperature scaling (fT)"),
#    col = c("steelblue", "firebrick"), lwd = 2, bty = "n" )
} 

# Takes DLI from env_daily df, calculates and plots associated scaling factor
plot_DLI_scaling <- function( env_dat ){
  
  # Days to plot
  DOY <- env_dat$real_DOY
  DLI       <- env_dat$DLI
  fL_daily  <- fL_reparam( DLI )
  
  # Set up plot margins to allow second y-axis
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  # --- Left axis: DLI ---
  plot( 
    DOY, DLI,
    type = "l", col = "steelblue", lwd = 2,
    xlab = "Day of Year", ylab = "Daily Light Interval (mol m⁻² d⁻¹)", bty = "l" )
  
  # --- Right axis: fL ---
  par(new = TRUE)
  
  plot(
    DOY, fL_daily,
    type = "l", col = "firebrick", lwd = 2,
    axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1) )
  
  axis(
    side=4, at=seq(0, 1, by=0.2), col = "firebrick", col.axis = "firebrick" )
  
  mtext(
    "DLI scaling factor (fL)",
    side = 4,
    line = 3,
    col = "firebrick"
  )
  
  # --- Legend ---
  #  legend(
  #    "topright",
  #    legend = c("Temperature", "Temperature scaling (fT)"),
  #    col = c("steelblue", "firebrick"), lwd = 2, bty = "n" )
} 


# Parameterless functions to plot the environmental curves.

plot_fT_curve <- function(T_opt = 10.0, warm_loss_per_C = 0.23,
                          T_min = 8, T_max = 25, T_step = 0.1) {
  
  df <- data.frame(Temp = seq(T_min, T_max, by = T_step))
  df$fT <- fT_reparam(df$Temp, T_opt = T_opt, warm_loss_per_C = warm_loss_per_C)
  
  ggplot(df, aes(x = Temp, y = fT)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = T_opt, linetype = "dashed") +
    labs(
      x = "Temperature (°C)",
      y = "Temperature scaling factor (fT)"
    ) +
    ylim(0, 1) +
    theme_classic()
}

plot_fL_curve <- function(
    L_low = 20,
    L_high = 40,
    low_loss_per_10 = 0.23,
    high_loss_per_10 = 0.23,
    DLI_min = 0,
    DLI_max = 70,
    DLI_step = 0.5
) {
  df <- data.frame(DLI = seq(DLI_min, DLI_max, by = DLI_step))
  df$fL <- fL_reparam(
    DLI = df$DLI,
    L_low = L_low,
    L_high = L_high,
    low_gain_per_10 = low_loss_per_10,
    high_loss_per_10 = high_loss_per_10
  )
  
  ggplot(df, aes(x = DLI, y = fL)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = c(L_low, L_high), linetype = "dashed") +
    labs(
      x = expression("Daily Light Integral (mol m"^-2*" d"^-1*")"),
      y = "Light scaling factor (fL)"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic()
}



#--- ChatGPT functions to scale T and DLI effects according to Pontier et al. ----

# --- Helper: convert "p% change per unit" into exponential rate constant ---
# If multiplier per unit is m (e.g., 0.77 per +1°C), then k = -log(m) / unit
k_from_multiplier <- function(multiplier, unit = 1) {
  -log(multiplier) / unit
}

# --- Temperature scaling: flat at/below optimum; exponential decay above optimum ---
# Based on ~23% loss per +1°C above the optimum for blade growth near ~11–12°C (Table 1) :contentReference[oaicite:1]{index=1}
fT_reparam <- function(T_C,
                       T_opt = 10.0,
                       warm_loss_per_C = 0.23) {  # 23% loss per +1°C
  # multiplier per +1°C on warm side
  m_warm <- 1 - warm_loss_per_C          # e.g., 0.77
  k_warm <- k_from_multiplier(m_warm, 1) # per °C
  
  ifelse(T_C <= T_opt,
         1.0,
         exp(-k_warm * (T_C - T_opt)))
}

# --- DLI scaling for blades: hump-shaped with low- and high-light penalties ---
# Approx ~23% gain per +10 DLI near ~20 (low-light side) and ~23% loss per +10 above ~40 (high-light side) :contentReference[oaicite:2]{index=2}
fL_reparam <- function(DLI,
                       L_low = 20,
                       L_high = 40,
                       low_gain_per_10 = 0.23,   # interpreted as "being 10 below optimal costs ~23%"
                       high_loss_per_10 = 0.23) {
  
  # Convert to multipliers:
  # low side: at (L_low - 10) -> multiplier ~ (1 - 0.23) = 0.77
  m_low10  <- 1 - low_gain_per_10
  k_low    <- k_from_multiplier(m_low10, 10)   # per (mol m^-2 d^-1)
  
  # high side: at (L_high + 10) -> multiplier ~0.77
  m_high10 <- 1 - high_loss_per_10
  k_high   <- k_from_multiplier(m_high10, 10)  # per (mol m^-2 d^-1)
  
  out <- ifelse(
    DLI < L_low,
    exp(-k_low * (L_low - DLI)),           # decays below L_low
    ifelse(
      DLI <= L_high,
      1.0,                                 # optimum plateau
      exp(-k_high * (DLI - L_high))         # decays above L_high
    )
  )
  
  # Keep bounded [0,1]
  pmax(0, pmin(1, out))
}












