#----------------------------------------------------------------------------
# Script:  LSSM_DEB.R
# Created: November 2025
## Purpose: Run a simple Dynamic Energy Budget model, written by ChatGPT
# Notes: Not a DEB, but a more sophisticated, 2-part logistic growth model
#----------------------------------------------------------------------------
# Updates
# Jan 2026:
#   Deterministic (i.e., logistic growth) model not very satisfying. Updated
#   this with a 'DEB-like' model that uses Ordinary Differential Equations to
#   predict growth based on light and temperature.
#   Nutrients are assumed to be not limiting - a significant assumption.
#   Uses data from Weigel and Pfister (2021) to partition and parameterize
#   carbon fixation, respiration, and DOC production.
#   Data from Pontier et al. (2024) used to attenuate growth as a function of
#   temperature and light levels.
# Oct 27: Restart after summer hiatus. Last modifications Jun 23.
#   Re-parameterized: stipe-frond allometry (Stekoll et al. 2006), maintenance
#   respiration (Fales et al. 2023), sloughing, temperature response, and
#   initial state.
#   - Nutrients: Trialed a Michaelis-Menten limiter and removed as no site nutrient 
#     data meant the limiter was controlled by two unconstrained parameters.
#     Nutrients remain NOT limiting.
#   - Temperature: The forcing from the mooring sensor falls within Pontier's 
#     2.5-3.1 m logger range (summer 9.3-15.7 vs 9.35-16.3 C) and agrees with the
#     deeper logger on the same mooring so no depth or regional correction is warranted.
#     fT is 1 at or below T_opt with a half-Gaussian above it. The sigma_warm is
#     calibrated to Pontier et al. (2024) Table 1. See tech report.
#   - Sloughing: half-cosine ramp in real_DOY. No sloughing (slough_min) before
#     slough_start, rising smoothly to slough_max at slough_end (zero slope at
#     both ends); frond only. Parameters are assumptions with no data
#     behind them, chosen so that growth bends down by mid-September.
#
# MODEL STRUCTURE (single 'super frond' + stipe; daily time step, rk4)
#   State variables (both g DW):  B_th = stipe,  B_fr = frond
#   Frond C budget (per day):
#     Net C = GP - DOC - R_prod - R_maint   (umol C/day)
#   converted to g DW and reduced by sloughing:
#     dB_fr = Net C (g DW) - slough_today * B_fr
#   Stipe: relaxes toward an allometric target set by B_fr; one-way (cannot
#   shrink), and its growth is NOT debited from the carbon budget.
#   Unit convention: everything is DW. WW = DW / wet_to_dry (output only).
#============================================================================

#--- REQUIRES
# deSolve (for ode)
# env_daily - a df describing growing conditions, one row per day, with columns
#   day, real_DOY, Temp (C), DLI (mol/m2/d, pre-reflectance), day_hours, night_hours.
#   This is made in the main script from dli_df and par_df, which come from LSSM_Water_Analysis.Rproj
# fT_reparam(), fL_reparam() - temperature and light limiters, defined in LSSM_config2.R
#   (Pontier et al. 2024: fT plateau at 1 up to T_opt, Gaussian decline above;
#    fL plateau between DLI 20 and 40, lower and higher DLI penalised)

#----- Helper functions

# Target stipe biomass (g DW) for a given frond biomass (g DW).
# Piecewise, based on the Stekoll et al. (2006) wet-weight regression for blade vs stipe mass:
#   large plants:  B_th = th_slope * B_fr + th_int
#   small plants:  B_th = th_max_ratio * B_fr (proportional), which joins the
#                  line at B_fr = th_B_cross (50 g WW). Below this the regression
#                  is an extrapolation (it does not pass through the origin).
stipe_target <- function(B_fr, pars) {
  pmax(1e-9, pmin(pars$th_max_ratio * B_fr,
                  pars$th_slope * B_fr + pars$th_int))
}


#----- Main function and parameters
#----------------------------------
# Jan 27: EVERYTHING in DW please.
# Define the nereo growth ODE system. We are growing dry weight using a single frond model.
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
    # Both limiters are multipliers on fixation, 0-1.
    # Nutrients are not limiting (no term).
    fT <- fT_reparam(T_C, T_opt, sigma_warm)
    fL <- fL_reparam(DLI)

    # --- STIPE growth using Nereo biomass relationships from Stekoll et al. (2006)
    # target stipe biomass
    B_th_target <- stipe_target(B_fr, pars)

    # Apply allometric relaxation to target (one-way: stipe cannot shrink).
    # Temperature-scaled. No carbon cost is charged for stipe growth.
    dB_th <- k_th * fT * max(0, B_th_target - B_th)

    # --- Mass-specific FROND fixation rates
    # (umol C per g DW per hour) -> umol C per g DW per day
    GP_day_umolC_gDW   <- max_fixation  * fT * fL * day_h
    GP_night_umolC_gDW <- dark_fixation * fT * night_h
    GP_umolC_gDW       <- GP_day_umolC_gDW + GP_night_umolC_gDW

    # Now calculate the Carbon fixed this day (umol C/day, whole frond; for DOC leakage, below)
    GP_day_umolC   <- B_fr * GP_day_umolC_gDW
    GP_night_umolC <- B_fr * GP_night_umolC_gDW
    GP_umolC       <- GP_day_umolC + GP_night_umolC

    # --- DOC Leakage - variable as a function of light (W&P)
    # Define a daily fraction
    DOC_day_frac <- DOC_day_min + (DOC_day_max - DOC_day_min) /
                    (1 + exp(DOC_DLI_k * (DLI - DOC_DLI50)))
    DOC_day_frac <- max(0, min(1, DOC_day_frac)) # Bound to [0,1]

    # Calculate the DOC
    DOC_day_umolC   <- DOC_day_frac   * GP_day_umolC
    DOC_night_umolC <- DOC_night_frac * GP_night_umolC
    DOC_umolC       <- DOC_day_umolC + DOC_night_umolC

    # --- Respiration (two components)
    # Maintenance: scales with temperature via Q10, and frond biomass plus a 
    # fraction of stipe biomass. NOTE: This is independent of fT which limits growth/fixation.
    # Stipe_resp_frac is the stipe's mass-specific maintenance rate relative to the frond
    # (1 = same as frond, 0 = none). No measurements, based on structural assumptions.
    fR <- Q10_resp ^ ((T_C - T_ref_resp) / 10)
    R_maint_umolC <- R_maint_umolC_gDW_h * fR * (B_fr + (stipe_resp_frac * B_th)) * 24

    # Production (growth) respiration: a fixed proportion of daytime fixation
    R_prod_umolC  <- respired_prod * GP_day_umolC
    
    # Net µmol fixed for this simulated day
    NetC_umolC <- GP_umolC - DOC_umolC - R_prod_umolC - R_maint_umolC

    # Convert total net C -> total biomass gain (DW)
    # µmol C -> mol C -> g C -> g DW
    tissue_growth <- NetC_umolC * 12.011 * 1e-6 / gDW_to_gC

    # --- Senescence: half-cosine (sinusoidal) variable rate ramp for sloughing.
    # slough_min until slough_start, rising smoothly to slough_max at slough_end
    # (zero slope at both ends). Applies to the frond only.
    ramp <- 0.5 * (1 - cos(pi * pmin(1, pmax(0, (real_DOY - slough_start) /
                                               (slough_end - slough_start)))))
    slough_today <- slough_min + (slough_max - slough_min) * ramp
    
    # Convert frond sloughing (g DW/day) to µmol C/day for C-budget plotting
    Slough_gDW <- slough_today * B_fr
    Slough_umolC <- Slough_gDW * gDW_to_gC / 12.011 * 1e6

    # --- Net Growth (Biomass Change) ---
    dB_fr <- tissue_growth - (slough_today * B_fr)

    # as.numeric() ensures name stays as defined. ode() oddity.
    # Output order: derivatives first, then diagnostics (these become columns in ode() output)
    return(list(
      c(dB_th, dB_fr),
      real_DOY = real_DOY,
      fL = fL,
      fT = fT,
      fR = fR,                                        # respiration temperature multiplier
      DLI = DLI,                                      # DLI after reflectance loss
      net_daily_gDW   = as.numeric(dB_fr+dB_th),      # total net g DW change (frond and stipe)
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
  max_fixation   = 120.78, # umol C/gDW/h. Top of W&P range of gross 13C fixation (48-121; mean 78.5)
  dark_fixation  = 3.75,   # Light-independent nighttime fixation (umol C/gDW/h) - Weigel and Pfister
  
  # Respiration 
  respired_prod  = 0.08,   # Production (growth) respiration, as fraction of daytime GP.
  R_maint_umolC_gDW_h = 6, # umol C/gDW/h at T_ref_resp. Fales et al. 2023 dark resp,
                           # 0.2-0.4 mg O2/g/h = 6-12 (RQ = 1).
  T_ref_resp     = 13,    # reference temperature (C); mid-range of Fales 13-16C data
  Q10_resp       = 1,     # 1 = no temperature dependence of respiration. Not fitted
  stipe_resp_frac = 0.33, # The stipe's mass-specific respiration rate relative to the frond
  
  # POM loss
  slough_min   = 0,      # d^-1; no sloughing of young, vigourous plants
  slough_start = 226,    # DOY (~Aug 15): sloughing begins to rise
  slough_end   = 285,    # DOY (~Oct 15): sloughing reaches slough_max
  slough_max   = 0.05,   # d^-1 (assumption, no data)

  # DOC leakage as fraction of GP during DAY: decreases with light (DLI)
  # Based on Weigel and Pfister
  DOC_day_max = 0.24,   # low-light PER (e.g., ~0.23–0.24)
  DOC_day_min = 0.08,   # high-light PER (e.g., ~0.08)
  DOC_DLI50   = 43,     # DLI at midpoint of decline. calc'd from DLI distributional data
  DOC_DLI_k   = 0.15,   # Steepness per DLI unit, calc'd from DLI distributional data

  # Night DOC leakage (simple assumption): fraction of night GP
  DOC_night_frac = 0.10, # tune later; keep small under non–nutrient-limited assumption

  # Temperature restriction on fixation rate based on Pontier et al.
  T_opt      = 10,    # C; optimum 9-11 C in Pontier
  sigma_warm = 2.4,   # width of the Gaussian decline above T_opt, calibrated to Pontier Table 1

  # Stipe dynamics (allometry: Stekoll et al. 2006; see stipe_target() above)
  th_slope = 0.3212,  # g DW stipe per g DW frond (regression slope after removing blades)
  k_th     = 0.3,     # Relaxation rate (per day) - enough to prevent stipe lag, but not instantaneous

  # Constants
  gDW_to_gC    = 0.313, # From RSSM.
  wet_to_dry   = 0.13,  # From RSSM. Applied to both stipe and frond
  alpha_water  = 0.06   # Fraction of DLI reflected at water surface
)

# Derived stipe allometry parameters (all g DW)
model_params$th_int       <- 0.1052 * 1000 * model_params$wet_to_dry  # 13.68 g DW: regression intercept (0.1052 kg WW)
model_params$th_B_cross   <- 50 * model_params$wet_to_dry             # 6.5 g DW: frond mass (50 g WW) where small-plant branch meets the line
model_params$th_max_ratio <- model_params$th_slope + model_params$th_int / model_params$th_B_cross  # 2.43: stipe:frond ratio for small plants


#----- RUNNNING the code ...
# 1. Define initial state (start with small kelp plant)
# Starting biomass (DW) based on young sporophyte with established first frond = 1 to 4 g WW total.
# B_fr0 is the frond mass; the stipe starts at its allometric target.
B_fr0 <- 15 * model_params$wet_to_dry

initial_state <- c(
  B_th = stipe_target(B_fr0, model_params),
  B_fr = B_fr0)

# 2. Define time sequence: one step per day of env_daily
times <- env_daily$day

# 3. Solve the ODE system using the 'ode' function
output <- ode(
  y = initial_state,
  times = times,
  func = grow_kelp_model,
  parms = model_params,
  env_data = env_daily,
  method = "rk4" # Runge-Kutta 4th order method
)

kelp_df <- as.data.frame( output )
kelp_df <- cbind( kelp_df, "B_th_WW" = kelp_df$B_th/model_params$wet_to_dry,
                           "B_fr_WW" = kelp_df$B_fr/model_params$wet_to_dry)


### Fin.
