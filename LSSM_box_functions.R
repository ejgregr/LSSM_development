#################################################################################
# Script:  LSSM_box_functions.R
# Created: June 3, 2025. EJG
# functions to support the box model
# ## Updates: 
# 2025/xx/xx: 
#################################################################################



applyVerticalMixing <- function(box, conc_surf, conc_bot, vol_surf, vol_bot, alpha_mix) {
  # Skip if either layer has zero volume
  if (vol_surf[box] == 0 || vol_bot[box] == 0) return(list(
    conc_surf = conc_surf[box],
    conc_bot = conc_bot[box]
  ))
  
  # Mixing strength for this box
  alpha <- alpha_mix[box]
  
  # Volume to exchange (based on smaller layer)
  v_ex <- alpha * min(vol_surf[box], vol_bot[box])
  
  # Mass exchanged (µmol)
  mass_s_to_b <- v_ex * conc_surf[box]
  mass_b_to_s <- v_ex * conc_bot[box]
  
  # Updated total mass in each layer
  mass_s_new <- conc_surf[box] * vol_surf[box] - mass_s_to_b + mass_b_to_s
  mass_b_new <- conc_bot[box]  * vol_bot[box]  - mass_b_to_s + mass_s_to_b
  
  # Updated concentrations
  conc_s <- mass_s_new / vol_surf[box]
  conc_b <- mass_b_new / vol_bot[box]
  
  return(list(conc_surf = conc_s, conc_bot = conc_b))
}



populateFlowMatrices <- function(){
  # Populate flow matrices
  for (i in seq_len(nrow(flow_df))) {
    from <- flow_df$From[i]
    to   <- flow_df$To[i]
    tide <- flow_df$Tide[i]
    rf_s <- flow_df$RF_surface[i]
    rf_b <- flow_df$RF_bottom[i]
    exch <- flow_df$exch_fraction[i]
    
    val_surf <- rf_s * exch
    val_bot  <- rf_b * exch
#    print( val_surf)
#    print( val_bot)
    
    if (tide == "ebb") {
      ebb_surf[from, to] <- val_surf
      ebb_bot[from, to]  <- val_bot
    } else if (tide == "flood") {
      flood_surf[from, to] <- val_surf
      flood_bot[from, to]  <- val_bot
    }
  }
  return(list(
    ebb_surf = ebb_surf,
    ebb_bot = ebb_bot,
    flood_surf = flood_surf,
    flood_bot = flood_bot
  ))
}