#----------------------------------------------------------------------------
# Script:  LSSM_DEB.R
# Created: November 2025
## Purpose: Run a simple Dynamic Energy Budget model, written by ChatGPT
# Notes: Not a DEB, but a more sophisticated, 2-part logistic growth model
#============================================================================

#---- FUNCTIONS ----

# Plot grams DW of fronds, stipe, and total
FullKelpPlot <- function ( kelp_out ) {
# Example: assuming your model output is called kelp_out
plot( as.Date(kelp_out$days), kelp_out$B_gDW,
      type = "l", lwd = 3, col = "darkgreen",
      xlab = "Date", ylab = "Biomass (g DW)",
      main = "Nereocystis Biomass Dynamics")

lines(as.Date(kelp_out$days), kelp_out$B_stipe_gDW,
      col = "brown", lwd = 2, lty = 2)

lines(as.Date(kelp_out$days), kelp_out$B_fronds_gDW,
      col = "forestgreen", lwd = 2, lty = 3)

legend("topleft",
       legend = c("Total", "Stipe", "Fronds"),
       col = c("darkgreen", "brown", "forestgreen"),
       lty = c(1, 2, 3), lwd = 2, bty = "n")
}

#--- Fourth model - Separates stipe and canopy growth.
grow_kelp4 <- function(df,
                       # ---- Biomass book-keeping ----
                       Bt_init = 0.005,              # g DW initial stipe
                       Bf_init = 0.5,                # g DW initial fronds (when triggered)
                       B_max_wet = 9230,             # g wet ~9.23 kg
                       wet_to_dry = 0.13,
                       dry_to_C = 0.25,
                       alpha = 1.0,                  # fraction of new dry mass that is carbon
                       
                       # ---- Growth kinetics ----
                       r_max_th = 0.2,
                       r_max_fr = 0.15,              # initial estimate = 0.065 based on earlier math
                       T_min = 4, T_opt = 8, T_max = 14,
                       k_L = 10,
                       
                       # ---- Logistic caps ----
                       stipe_max_frac = 0.20,
                       frond_start_gDW  = 1,          # stipe mass needed to begin fronds
                       
                       # ---- Losses ----
                       slough_rate_th = 0.0,
                       slough_rate_fr = 0.01
) {
  
  n <- nrow(df)
  stopifnot(all(c("days","temp","salt","DLI") %in% names(df)))
  
  # Convert caps
  B_max_dry <- B_max_wet * wet_to_dry        # total DW
  B_th_cap  <- stipe_max_frac * B_max_dry  # max stipe DW
  
  # Storage
  B_th <- numeric(n)
  B_fr <- numeric(n)
  
  # ---- Initialize correctly ----
  B_th[1] <- min(Bt_init, B_th_cap)
  B_fr[1] <- 0   # fronds start later
  
  # Tracking
  r_eff_th <- numeric(n)
  r_eff_fr <- numeric(n)
  fT_vec   <- numeric(n)
  fL_vec   <- numeric(n)
  
  # Limiter functions
  fT <- function(T){
    if (T <= T_min || T >= T_max) return(0)
    ((T - T_min)/(T_opt - T_min)) * ((T_max - T)/(T_max - T_opt))
  }
  fL <- function(I){ I / (I + k_L) }
  
  # ---- Time loop ----
  for (i in 2:n) {
    T  <- df$temp[i]
    I  <- df$DLI[i]
    FT <- fT(T)
    FL <- fL(I)
    
    # --- stipe growth ---
    r_th <- r_max_th * FT
    growth_th <- B_th[i-1] * r_th * (1 - B_th[i-1] / B_th_cap)   # logistic
    loss_th   <- slough_rate_th * B_th[i-1]
    B_th[i]   <- max(0, min(B_th_cap, B_th[i-1] + growth_th - loss_th))
    
    # --- Frond growth ---
    frond_on <- B_th[i] >= frond_start_gDW
    
    # Trigger fronds EXACTLY once
    if (frond_on && B_fr[i-1] == 0) {
      B_fr[i-1] <- Bf_init
    }
    
    B_total_prev <- B_th[i-1] + B_fr[i-1]
    
    logistic_term <- (1 - B_total_prev / B_max_dry)
    logistic_term <- max(0, logistic_term)   # prevent negative logistic
    
    r_fr <- if (frond_on) r_max_fr * FT * FL else 0
    growth_fr <- B_fr[i-1] * r_fr * logistic_term
    loss_fr   <- slough_rate_fr * B_fr[i-1]
    
    # enforce total biomass cap
    B_fr_raw <- B_fr[i-1] + growth_fr - loss_fr
    B_fr[i] <- max(0, min(B_max_dry - B_th[i], B_fr_raw))
    
    # Store growth modifiers
    r_eff_th[i] <- r_th
    r_eff_fr[i] <- r_fr
    fT_vec[i]   <- FT
    fL_vec[i]   <- FL
  }
  
  # ---- Build output ----
  df$B_stipe_gDW <- B_th
  df$B_fronds_gDW  <- B_fr
  df$B_gDW         <- B_th + B_fr
  df$B_kgWW        <- (df$B_gDW / wet_to_dry) / 1000
  
  df$r_eff_th <- r_eff_th
  df$r_eff_fr <- r_eff_fr
  df$fT       <- fT_vec
  df$fL       <- fL_vec
  
  df
}


