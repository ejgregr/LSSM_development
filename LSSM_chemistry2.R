#################################################################################
# Script: LSSM_Chemistry2.R
# Created: During second pass thru the chemistry code.
# Purpose: ----- Version 2 - Examine chemistry based on output of plant growth model -----
# Notes:
# - Early 2026 version expanded on the simplistic prototype and supported the first team model. 
#----------------------------------------------------------------------------
# Updates
# 26Nov: Fully reviewed and revised for final contract report.  

#################################################################################

#---- Required inputs:
#str( kelp_df )    # simulated Kelp growth, from LSSM_DEB2.R
#str( disc_data )  # created as part of output from Water_Analysis.proj
#str( env_daily )  # Required for the tea simulation

# Get rid of vestigial stuff to support debugging ... 
#rm(list = c("x", "sims", "inout_dic", "pairs", "sim_out" ))

#----------------------------------------------------------------------------
# Workflow:
# 1. Add some chemistry numbers
# 2. Combine inside and outside samples for the kelp stations
#    These are for station by station comparisons on the 'big plot'
# 3. Prepare the tea-making landscape 
#----------------------------------------------------------------------------

# Remove early samples prior to mooring deployment
disc_data <- subset(disc_data, as.Date(DateTime) >= as.Date( start_date ) & 
                              as.Date(DateTime) <= as.Date( end_date ))

#-----  1. Complete the chemistry for the discrete water samples ----
chem_data <- add_DIC_Alk( disc_data )
str(chem_data)


#---- Temporarly patch TA and DIC while awaiting info on params and 15 Batch from Wiley
chem_data$TA  <- chem_data$TA_lab       # lab TA (umol/kg)
chem_data$DIC <- chem_data$TCO2_lab     # lab Adj. TCO2 (umol/kg)
chem_data$pH  <- ph_from_ta_dic(chem_data$TA, chem_data$DIC, chem_data$S_lab, chem_data$T_lab)
#---- End of discrete data patch

summary(chem_data$pH - chem_data$pH_lab)

# For downstream compatability. I'm afraid I just 'overbuilt' the patch. :\
chem_data$Temp     <- chem_data$T_lab
chem_data$Salinity <- chem_data$S_lab
chem_data$pCO2     <- chem_data$pCO2_lab
chem_data$OArag    <- chem_data$OArag_lab

# Show the TA/Salinity analysis. m_TA is the regression fit.
m_TA <- fit_TA_salinity(chem_data)

#-----  2. Organise the data by in/out pairs  ----
# Pair In/Out by Site + Date
# NB: CP_In in May is a replicate and an odd duck in the pairs DF.
pairs <- make_pairs_out_in( chem_data )

#---- 3. Create the simulation landscape ----
# This is a grid of kelp density vs. exposure, used to make 'tea' by steeping
# a variable number of kelp plants in a fixed volume of water (the density)
# Exposure is the time the steep time - no mixing or flushing at this point.


#*** Latest work is in test_stuff(). 
# GONE BACK to water analysis to check for TOC2 ( Current

# Requires pairs2, kelp_df, env_daily)
sim_out <- simulate_dens_landscape(
  pairs_df = pairs,
  kelp_df = kelp_df,
  env_daily = env_daily,
  kelp_density_plants_m2_vec = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
  exposure_hours_vec         = c(1, 2, 3)
#  exposure_hours_vec         = c(1, 2, 4, 6, 8, 10, 12, 16, 20, 24)
)

unique( sim_out$sims$Site )

# Available dates for specified site
date_by_site <- sim_out$sims %>%
  dplyr::filter(Site == 'CP') %>%
  dplyr::distinct(Date) %>%
  dplyr::arrange(Date)

str( sim_out$sims )


# -------------------- Early density plotting routines -------------------

#----- Density-based plot for selected station and date
plot_density_exposure_dpH( sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22") )

# Interpolate the desired surface ...
x <- make_interp_surface(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22") )
y <- make_obs_dpH(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22") )
plot_interp_density_exposure_dpH( x, y ) 

# The precursor plot to the BigP with all stations over time. 
plot_interp_density_exposure_by_date(
  sims = sim_out$sims, site_pick = "CP" )



#---- FIN. ----
