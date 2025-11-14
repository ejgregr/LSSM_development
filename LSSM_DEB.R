#----------------------------------------------------------------------------
# Script:  LSSM_DEB.R
# Created: November 2025
## Purpose: Run a simple Dynamic Energy Budget model, written by ChatGPT
# Notes: Not a DEB, but a more sophisticated, 2-part logistic growth model
#============================================================================


kelp_grow <- grow_kelp4(daily_ocean, 
                            r_max_th = 0.2,
                            r_max_fr = 0.15,
                            thallus_max_frac = 0.2,
                            frond_start_gDW = 1)

plot(kelp_grow$days, kelp_grow$B_kgWW, type="l", ylab="Biomass (kg wet)", xlab="Date")

FullKelpPlot( kelp_grow )

# Examine growth inhibitors
matplot(kelp_grow$days[-1], kelp_grow[-1,c("fT","fL")], type="l", lwd=2,
        col=c("red","blue"), ylab="Modifier (0–1)", xlab="Date")
legend("topright", legend=c("Temp limit","Light limit"),
       col=c("red","blue"), lty=1, lwd=2)


names( kelp_grow )
range( kelp_grow$days )

#---- FUNCTIONS ----

# Plotting function to see fronds and thallus
FullKelpPlot <- function ( kelp_out ) {
# Example: assuming your model output is called kelp_out
plot( as.Date(kelp_out$days), kelp_out$B_kgWW,
      type = "l", lwd = 3, col = "darkgreen",
      xlab = "Date", ylab = "Biomass (g DW)",
      main = "Nereocystis Biomass Dynamics")

plot( as.Date(kelp_out$days), kelp_out$B_gDW,
      type = "l", lwd = 3, col = "darkgreen",
      xlab = "Date", ylab = "Biomass (g DW)",
      main = "Nereocystis Biomass Dynamics")

lines(as.Date(kelp_out$days), kelp_out$B_thallus_gDW,
      col = "brown", lwd = 2, lty = 2)

lines(as.Date(kelp_out$days), kelp_out$B_fronds_gDW,
      col = "forestgreen", lwd = 2, lty = 3)

legend("topleft",
       legend = c("Total", "Thallus", "Fronds"),
       col = c("darkgreen", "brown", "forestgreen"),
       lty = c(1, 2, 3), lwd = 2, bty = "n")
}

#--- Fourth model - Separates thallus and canopy growth.
grow_kelp4 <- function(df,
                   # ---- Biomass book-keeping (your originals) ----
                   Bt_init = 0.025,              # g DW for initial sporophite 
                   Bf_init = 0.5,                  # g DW to initialize frond growth 
                   B_max_wet = 9230,             # g wet adult ~9.23 kg (total)
                   wet_to_dry = 0.13,
                   dry_to_C = 0.25,              # g C / g DW
                   alpha = 1.0, #0.7,                  # fraction of new dry mass that is newly fixed C
                   # ---- Growth kinetics (split) ----
                   r_max_th = 0.06,              # d^-1 thallus max specific growth
                   r_max_fr = 0.08,              # d^-1 frond  max specific growth
                   T_min = 4, T_opt = 8, T_max = 14,  # temp triangle (slightly cooler opt for earlier growth)
                   k_L = 10,                     # DLI half-saturation (mol m^-2 d^-1) for fronds
                   # ---- Logistic caps & trigger ----
                   thallus_max_frac = 0.15,      # thallus maximum fraction of B_max_dry
                   frond_start_gDW = 3,         # fronds start once thallus ≥ this (g DW)
                   # ---- Losses ----
                   slough_rate_th = 0.003,       # d^-1 thallus background loss
                   slough_rate_fr = 0.012 )       # d^-1 frond loss (higher, blades turn over)
{
  n <- nrow(df)
  if (!all(c("days","temp","salt","DLI") %in% names(df))) {
    stop("df must contain columns: days, temp, salt, and DLI")
  }
  
  # Convert caps to dry weight
  B_max_dry <- B_max_wet * wet_to_dry                  # g DW (total)
  B_th_cap  <- thallus_max_frac * B_max_dry            # g DW (thallus-only cap)
  
  # Initialise state
  B_th <- numeric(n)   # thallus DW
  B_fr <- numeric(n)   # frond   DW
  B_th[1] <- min(B_init, B_th_cap)   # start all in thallus
  B_fr[1] <- max(0, B_init - B_th[1])
  
  r_eff_th <- numeric(n)
  r_eff_fr <- numeric(n)
  fT_vec   <- numeric(n)
  fL_vec   <- numeric(n)
  
  # Helper limiters
  fT <- function(T){
    if (T <= T_min || T >= T_max) return(0)
    ((T - T_min)/(T_opt - T_min)) * ((T_max - T)/(T_max - T_opt))
  }
  fL <- function(I){ I / (I + k_L) }
  
  # Time loop
  for (i in 2:n) {
    T <- df$temp[i]
    I <- df$DLI[i]
    FT <- fT(T)               # temperature limitation 0–1
    FL <- fL(I)               # light limitation 0–1 (for fronds)
    
    # --- Thallus growth (primarily temperature-driven; weakly tied to light implicitly through total cap) ---
    r_th <- r_max_th * FT
    growth_th <- B_th[i-1] * r_th * (1 - B_th[i-1] / B_th_cap)
    loss_th   <- slough_rate_th * B_th[i-1]
    dB_th     <- growth_th - loss_th
    B_th[i]   <- max(0, min(B_th_cap, B_th[i-1] + dB_th))
    
    # --- Frond growth (light + temperature + total logistic cap) ---
    # Trigger frond growth when thallus is established
    frond_on <- B_th[i] >= frond_start_gDW
    
    if (frond_on && B_fr[i-1] == 0) {
      B_fr[i-1] <- Bf_init    # Initial frond mass in g DW
    }
    
    # Set frond growth rate for this time step
    r_fr <- if (frond_on) r_max_fr * FT * FL else 0
    
    B_tot_prev <- B_th[i-1] + B_fr[i-1]
    growth_fr  <- B_fr[i-1] * r_fr * (1 - B_tot_prev / B_max_dry)
    loss_fr    <- slough_rate_fr * B_fr[i-1]
    dB_fr      <- growth_fr - loss_fr
    B_fr[i]    <- max(0, min(B_max_dry - B_th[i], B_fr[i-1] + dB_fr))
    
    # --- Carbon uptake (only from net positive ΔB_total) ---
    dB_total <- (B_th[i] + B_fr[i]) - (B_th[i-1] + B_fr[i-1])
    C_umol <- 0
    if (dB_total > 0) {
      C_g   <- dB_total * dry_to_C * alpha   # g C fixed
      C_umol <- C_g * 1e6 / 12               # µmol C
    }
    
    r_eff_th[i] <- r_th
    r_eff_fr[i] <- r_fr
    fT_vec[i]   <- FT
    fL_vec[i]   <- FL
  }
  
  # Pack outputs (totals + parts)
  df$B_thallus_gDW <- B_th
  df$B_fronds_gDW  <- B_fr
  df$B_gDW         <- B_th + B_fr
  df$B_kgWW        <- (df$B_gDW / wet_to_dry) / 1000
  
  # growth modifiers (optional but helpful)
  df$r_eff_th      <- r_eff_th
  df$r_eff_fr      <- r_eff_fr
  df$fT            <- fT_vec
  df$fL            <- fL_vec
  
  return(df)
}


