#################################################################################
# Script:  LSSM_main2.R
# Created: January 13 2026. EJG
# 
# Main revision of the earlier LSSM, focused on using 2025 data. 
#
## Updates: 
# Newly available data include
#   1. currents (speed and direction) for two tidal model sites
#   2. DLI from sattelite to better reflect growing conditions
#       NOTE: this fixed the DLI which now matches the values reported by Pontier.
#
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
# Loads temperature, salinity, CO2, and oxygen from moorings
source( "LSSM_config2.R" )

#---- Environmental data prep --------
DLI_env  <- data.frame( "Datestamp" = DLI_df$Date, "DLI" = DLI_df$mean )
temp_env <- data.frame( "Datestamp" = date( DST_focal1$DateTime ), "Temp" = DST_focal1$Temp )
PAR_env  <- data.frame( "Datestamp" = par_df$Timestamp, "PAR" = par_df$mean, "SD" = par_df$stdDev )
env_daily <- make_env_daily( DLI_env, temp_env )
head(env_daily)

# Daylight hours for the days in env_daily 
day_hours <- get_daylight_hours( PAR_env )

# Join daylight hours to available daily environmental data 
env_daily <- env_daily %>%
  left_join(day_hours, by = "Datestamp")


#---- Nereo growth parameters (older) ----
# B_init      <- 0.025  # Initial mass of sporophyte (est. at 25 mg based on ChatGPT)
# B_max       <- 9.23   # Max biomass for an adult Nereo. Calculated from field results (Weigel and Pfister 2021)
# r_max       <- 0.065  # Maximum daily growth rate, assuming 6 month growth and mature plant is 9.23 kg 
# sp_start    <- 0.005  # Established sporophyte mass (5 g)
# wet_to_dry  <- 0.13   # Water content of Nereo (Bullen et al. 2024)
# dry_to_C    <- 0.25   # Carbon content of Nereo (dry) (Bullen et al. 2024)
# 
# #---- Environmental influence on growth rate
# T_opt     <- 10        # Optimal temperature (°C)
# T_max     <- 14        # Maximal temperature for growth (°C)
# DLI_opt   <- 30        # Optimal DLI in mol/m2/day
# DLI_range <- 20        # Range (+/-) that allows growth in mol/m2/day 




#-------------- Run the ODE model  ---------------------------------
source( "LSSM_ODE.R" )

#----- First view of output  ----- 
out_df <- as.data.frame( output )
bio_df <- get_biomass_timeseries( out_df )
kelp_plot <- plot_kelp_biomass(bio_df)
print(kelp_plot)



#-------------- Visualization and diagnostics ---------------------------------

#Ensure output has multiple blades 
ncol( output )


plot_individual_tissues(output, env_daily)


first_blade_day1 <- output[2, 3]
last_blade_day1  <- output[2, ncol(output)]

final_row <- output[nrow(output), ]
manual_frond_sum <- sum(final_row[3:ncol(output)])
stipe_val <- final_row[2]

#----- Visualization using output as a df ----- 
out_df <- as.data.frame( output )

#--- 
bio_df <- get_biomass_timeseries( out_df )
kelp_plot <- plot_kelp_biomass(bio_df)
print(kelp_plot)


#--- Show blade distribution  ---
p_dist <- plot_blade_distribution( out_df )
print(p_dist)


# --- Plot kelp with sloughing ---
# Add real DOY to bio_df
start_day <- env_daily$day[1] # e.g., 135
bio_df$Real_DOY <- bio_df$time + start_day

plot_kelp_dynamics(bio_df, env_daily, sen_day = 240)






#---- Part 3: Calculate chemistry changes during  growing season ----
source( "LSSM_chemistry.R" )





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
