#----- Version 2 - Examine chemistry based on output of plant growth model -----
# 2026/02/09. Close to done. 

#--- Get rid of vestigial stuff to support debugging
#rm(list = c("x", "sims", "inout_dic", "pairs", "sim_out" ))

#---- Required inputs:
str( kelp_df )    # simulated Kelp growth
str( inout_data ) # loaded in LSSM_config2.R



# =============================================================================
# Workflow:
# =============================================================================

# (A) Add DIC to discrete water samples
inout_dic <- add_DIC_no_grouping(inout_data)

# (B) Pair In/Out by Site + Date
# NB: CP IN on May is a replicate and an odd duck in the pairs DF.
pairs <- make_pairs_out_in(inout_dic)

#----------------- Create the landscape ------------------------------
# (C) Run response landscape (requires kelp_df + env_daily)

# Either with water and kelp mass separately ...  
# sim_out <- simulate_landscape(
#    pairs_df = pairs,
#    kelp_df = kelp_df,
#    env_daily = env_daily,
#    water_mass_kg_vec = c(10, 25, 50),
#    kelp_biomass_gdw_vec = c(5, 10, 15),
#    exposure_hours_vec = c(6, 12, 24)
#  )

# Or (simpler) with just density  ...  
# BUT, needs reasonable density. 
sim_out <- simulate_dens_landscape(
  pairs_df = pairs,
  kelp_df = kelp_df,
  env_daily = env_daily,
  kelp_density_gdwkg = c(0.1, 1, 5, 10, 20, 50, 100, 500),
  exposure_hours_vec = c(1, 6, 9, 12, 18, 24),
  water_mass_kg = 1
)

unique( sim_out$sims$Site )

# Available dates for specified site
date_by_site <- sim_out$sims %>%
  dplyr::filter(Site == 'CP') %>%
  dplyr::distinct(Date) %>%
  dplyr::arrange(Date)

# Calculate match error between observed inside and modelled dpH

plot_landscape_dpH_facets( sim_out$sims, date_pick = as.Date("2025-07-11"), exposure_pick = 24 )



#----- Density-based plots ... ----
plot_density_exposure_dpH( sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22") )



# ---- Initial density plotting routines -----

# Interpolate the desired surface ...
x <- make_interp_surface(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22") )
y <- make_obs_dpH(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22") )
plot_interp_density_exposure_dpH( x, y ) 


x <- "CP"
plot_interp_density_exposure_by_date(
  sims = sim_out$sims,
  site_pick = x
)



#---- FIN. ----
