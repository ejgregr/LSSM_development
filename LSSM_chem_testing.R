library(seacarb)

# ---------------------------------------------------------
# FUNCTION: simulate_kelp_pH()
# ---------------------------------------------------------
# Inputs:
#   daily_ocean  = df with columns: days, temp, salt, pCO2, totAlk
#   kelp_grow    = df with column B_gDW (total dry biomass)
#   volume_kg    = mass of seawater affected (kg)
#   C_frac       = carbon fraction of kelp dry weight (default = 0.33)
#
# Output:
#   dataframe with:
#      day, T, S, TA, DIC_with_kelp, pH_with_kelp, pH_no_kelp
# ---------------------------------------------------------


flag_CA <- 24 # pCO2 and Alkalinity 
flag_CD <- 25 # pCO2 and DIC  
flag_AD <- 15 # Alkalinity and DIC 

result <- simulate_kelp_pH(
  daily_ocean = daily_ocean,
  kelp_grow   = kelp_grow,
  volume_kg   = 50   # choose any water volume
)

tail(result)



simulate_kelp_pH <- function(daily_ocean, kelp_grow,
                             volume_kg,
                             C_frac = 0.33) {
  
  stopifnot(all(daily_ocean$days == kelp_grow$days))
  n <- nrow(daily_ocean)
  
  # compute net daily biomass production (gDW)
  dB <- c(0, diff(kelp_grow$B_gDW))
  
  # convert to grams C
  dC_g <- dB * C_frac
  
  # convert to moles C ( = moles DIC )
  dC_mol <- dC_g / 12.01
  
  # convert to moles DIC removed per kg of water
  drawdown_mol_kg <- (dC_mol) / volume_kg
  
  # output storage
  out <- data.frame(
    day = daily_ocean$days,
    T = daily_ocean$temp,
    S = daily_ocean$salt,
    TA = daily_ocean$totAlk,
    DIC_no_kelp = NA,
    DIC_with_kelp = NA,
    pH_no_kelp = NA,
    pH_with_kelp = NA,
    delta_pH = NA
  )
  
  for (i in 1:n) {
    
    S <- daily_ocean$salt[i]
    T <- daily_ocean$temp[i]
    TA <- daily_ocean$totAlk[i]
    pCO2 <- daily_ocean$pCO2[i]
    
    # baseline carbonate system
    base <- carb(flag = 24, var1 = TA, var2 = pCO2, S = S, T = T)
    
    out$DIC_no_kelp[i] <- base$DIC
    out$pH_no_kelp[i]  <- base$pH
    
    # DIC removed today (in mol/kg)
    dd_mol_kg <- drawdown_mol_kg[i]
    
    # new DIC
    DIC_new <- base$DIC - dd_mol_kg
    if (DIC_new < 0) DIC_new <- 0   # safety
    
    out$DIC_with_kelp[i] <- DIC_new
    
    # pH after DIC removal
    post <- carb(flag = 15, var1 = TA, var2 = DIC_new, S = S, T = T)
    
    out$pH_with_kelp[i] <- post$pH
  }
  
  out$delta_pH <- out$pH_with_kelp - out$pH_no_kelp
  
  out
}
