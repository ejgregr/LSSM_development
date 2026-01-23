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

#Add real DOY for plotting etc.
env_daily <- cbind(env_daily, "real_DOY" = as.integer(format(as.Date(env_daily$Datestamp), "%j")) )
#  real_DOY <- env_data$day[row_idx]

# Daylight hours for the days in env_daily 
day_hours <- get_daylight_hours( PAR_env )
#plot(env_daily$DLI)

# Join daylight hours to available daily environmental data 
env_daily <- env_daily %>%
  left_join(day_hours, by = "Datestamp")



head(env_daily)
#-------------- Load and run the ODE model  ---------------------
source( "LSSM_ODE.R" )
out_df <- as.data.frame( output )
out_df$time <- env_daily$real_DOY

names( out_df )
out_df[1:20,]
#-------------- Visualization and diagnostics ---------------------------------
#----- First view of output  ----- 
#plot_kelp_biomass(out_df ) 
plot_kelp_biomass(out_df, log=T ) 



#--- Diagnostics

plot(out_df$fL)

dt <- c(NA, diff(out_df$time))

plot(out_df$net_specific)

cum_log_growth <- cumsum(ifelse(is.na(dt), 0, out_df$net_specific * dt))

plot(cum_log_growth)

head(out_df)
bio_df <- out_df
# --- Plot kelp with sloughing ---
# Add real DOY to bio_df
start_day <- env_daily$day[1] # e.g., 135
bio_df$Real_DOY <- bio_df$time + start_day

plot_kelp_dynamics(bio_df, env_daily, sen_day = 240)


#--- Check slough rate vs. blade light limitation 
plot(out_df$time, out_df$fL, type="l", xlab="day", ylab="fL_blade")
plot(out_df$time, out_df$slough, type="l", xlab="day", ylab="current_slough (1/day)")
plot(out_df$time, out_df$Growth_Rate_Fr - out_df$slough, type="l",
     xlab="day", ylab="net specific rate (1/day)")
abline(h=0, lty=2)



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
