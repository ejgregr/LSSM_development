#################################################################################
# Script:  LSSM_main2.R
# Created: January 13 2026. EJG
# 
# Main revision of the earlier LSSM, focused on using 2025 data. 
#
## Updates: 
# Newly available data include
#   1. currents (speed and direction) for two tidal model sites
#   2. DLI from satellite to better reflect growing conditions
#       NOTE: this fixed the DLI which now matches the values reported by Pontier.
# Jan29: Nice stable, parameterized version of LSSM 'DEB'. Easy to run and 
#   ready for chemistry.
#################################################################################

rm(list=ls(all=T)) # Erase environment.
set.seed(42)       # Setting seed for reproducibility

# Load packages and functions ... 
setwd( "c:/Data/Git/LSSM_development")
source( "LSSM_config2.R" )

# Set temporal window to match moorings ... 
start_date <- "2023-05-04 00:00:00 PDT"
end_date   <- "2023-10-11 00:00:00 PDT"

#---- Part 1: Load and prep environmental drivers ----
# Loads temperature, salinity, CO2, and oxygen from moorings, 
# as well as light (as DLI and PAR).
source( "LSSM_config2.R" )

#---- Environmental data prep --------
DLI_env  <- data.frame( "Datestamp" = DLI_df$Date, "DLI" = DLI_df$mean )
temp_env <- data.frame( "Datestamp" = date( DST_focal1$DateTime ), "Temp" = DST_focal1$Temp )
PAR_env  <- data.frame( "Datestamp" = par_df$Timestamp, "PAR" = par_df$mean, "SD" = par_df$stdDev )
env_daily <- make_env_daily( DLI_env, temp_env )
head(env_daily)

#Add real DOY for plotting etc.
env_daily <- cbind(env_daily, "real_DOY" = as.integer(format(as.Date(env_daily$Datestamp), "%j")) )

# Daylight hours for the days in env_daily 
day_hours <- get_daylight_hours( PAR_env )
#plot(env_daily$DLI)

# Join daylight hours to available daily environmental data 
env_daily <- env_daily %>%
  left_join(day_hours, by = "Datestamp")

head(env_daily)

# Test temp effect on growth
env_daily$Temp <- env_daily$Temp - 2
plot( env_daily$Temp )

#----------------------- Load and run the ODE model  ---------------------------
# uses dataframe env_daily
source( "LSSM_DEB.R" )
# creates dataframe kelp_df
#---------------------- Visualization and diagnostics --------------------------
names( kelp_df )
plot_kelp_biomass_WW(kelp_df, log=F ) 
plot( kelp_df$net_daily_gDW)

head( kelp_df )

# Jan29: Reasonably well parameterized, except DOC leakage and respiration
# fall of the cliff with higher temps (beyond 12C). This is new with the Gaussian
# updated fT curve. But it is also reasonable, and logical give the T constraint.
# Check into this later, with respect to field temps and moorings. 
# Chat says to Increase sigma_warm until: fT ≈ 0.3–0.5 at late-season temperatures

plot_fT_curve( model_params$T_opt, model_params$sigma_warm )
plot_fL_curve()

plot_C_losses(kelp_df, log_y=T )
plot_temperature_scaling( env_daily )
plot_DLI_scaling( env_daily )
summary(fT_reparam(env_daily$Temp))


#---- Part 3: Calculate chemistry changes during  growing season ----
source( "LSSM_chemistry2.R" )





#------------ Plot the BIG plot of tea landscape with paired samples -------

# Interpolate all the surfaces from the simulation ...
surf_all <- build_interp_grid_data(sim_out$sims, nx = 200, ny = 200, z_col = "dpH_end")

# Find all available contours ... 
# This assumes a global variable 
contour_ok <- summarise_contour_availability_all_sites(sim_out$sims, nx = 200, ny = 200, z_col = "dpH_end")

# Filter out the surfaces with no contour ... 
# NOTE: 
surf_all_ok <- surf_all %>%
  dplyr::inner_join(
    contour_ok %>%
      dplyr::filter(has_contour) %>%
      dplyr::select(Site, Date),   # <- keep only the keys from contour_ok
    by = c("Site", "Date")
  )

# Create label positions for dpH values ...
label_df <- surf_all_ok %>%
  dplyr::group_by(Site, Date) %>%
  dplyr::summarise(
    dpH_obs = dpH_obs[1],
    x_lab = min(exposure_hours, na.rm = TRUE),
    y_lab = max(kelp_density_gdwkg, na.rm = TRUE),
    .groups = "drop"
  )


# Plot only those with contours ... 
bigP <- plot_interp_grid_sites_dates(surf_all_ok, z_col = "dpH_end", log10_density = TRUE)

#--- PLOT all incubations ... 
#bigP <- plot_interp_grid_sites_dates(surf_all, z_col = "dpH_end", log10_density = TRUE)

ggsave(
  "deltapH_and_density_exposure_grid.png",
  plot = bigP,
  width = 14,          # keep width similar
  height = 25,         # ~3× taller than default
  units = "in",
  dpi = 300
)











#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 

Sys.setenv(TEXLIVE_AUTO_UPDATE = "0")

today <- format(Sys.Date(), "%Y-%m-%d")

rmarkdown::render( "LSSM_documentation.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = PDF_dir,
                   output_file = paste0( "LSSM_testing_", today ))

# FIN.
