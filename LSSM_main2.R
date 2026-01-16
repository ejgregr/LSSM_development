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


#---- Part 1: Load environmental drivers ----
# Loads temperature, salinity, CO2, and oxygen from moorings
source( "LSSM_drivers.R" )





#---- Show Hourly temp and salinity for BATI moorings 5 and 7 ----
BATI5$Mooring <- 'BATI5'
BATI7$Mooring <- 'BATI7'
plotme <- rbind( BATI5, BATI7)

# Plot full (hourly) temperature and salinity time series
PlotBATI( plotme, var='temp', head='Mooring time series of temperature (C)')
PlotBATI( plotme, var='salt', head='Mooring time series of salinity (psu)') 

# Add mooring ID for plotting daily mooring time series and merge
dBATI5$Mooring <- 'BATI5'
dBATI7$Mooring <- 'BATI7'
plotme <- rbind( dBATI5, dBATI7)

PlotBATI(plotme, var='temp_mean', head='Mooring daily time series of temperature (C)')
PlotBATI(plotme, var='salt_mean', head='Mooring daily time series of salinity (psu)')

#---- Create and show alkalinity from salinity for BATI moorings 5 and 7----
# Using equation from Evans et al. 2015:  TA = 48.7709*S + 606.23 (μmol kg-1)

daily_TA <- data.frame( date = dBATI5$datestamp )
daily_TA$B5 <- CalcAlk( dBATI5$salt_mean )
daily_TA$B7 <- CalcAlk( dBATI7$salt_mean )

plot( daily_TA$date, daily_TA$B5, type='l', main= "Average, Daily Total Alkalinity (BATI5)",
      xlab="", ylab="mol / kg" )

plot( daily_TA$date, daily_TA$B7, type='l', main= "Average, Daily Total Alkalinity (BATI7)",
      xlab="", ylab="mol / kg" )

#---- Part 1b: Create and Add solar input ----
# Temporal resolution and extents depend on mooring data
# Light data available from sensors but needs cleaning up. Simulate for now.

BATI5 <- AddSurfaceLight( df = BATI5, lat =latitude, lon = longitude,
                      datetime_col = 'datestamp', tz = bc_time )

# Now update the daily data frame with daily light index
# Function returns df of [date, DLI]
x <- ParToDLI(df = BATI5, datetime_col = "datestamp", 
              par_col = "PAR_umol_m2_s", tz = bc_time)
dBATI5 <- cbind( dBATI5, "DLI" = x[[2]] )
dBATI7 <- cbind( dBATI7, "DLI" = x[[2]] )



#---- Show hourly PAR and DLI ----
par(mfrow=c(2,1), mar=c(2,4,2,1) )
plot(BATI5$datestamp, BATI5$PAR_umol_m2_s,
     type = "l", xlab = "Time", ylab = "Surface PAR (µmol m⁻² s⁻¹)",
     main= "Hourly Photosynthetically Active Radiation (simulated)")

plot( dBATI5$date, dBATI5$DLI, type='l', main= "Daily Light Intensity (simulated)", 
      xlab="", ylab="mol photons / m² d" )


#---- Part 1c: Show carbon data from the MV Columbia before and after interpolation
plot(ak_dat$date, ak_dat$CO2mn, type = "o", main = "Full Ak Ferry daily seawater pCO2 climatology",
     xlab = "Date", ylab = "uatm", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = ak_dat$date, format = "%m-%d")

plot(daily_swCO2$Date, daily_swCO2$dCO2, type = "o", main = "Summer Ak Ferry daily pCO2 - interpolated",
     xlab = "Date", ylab = "uatm", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = daily_swCO2$Date, format = "%m-%d")


#---- Build an Ocean df to house all the data ----

# All carb() related values should be in carb() units (mol/kg) and day length:
daily_ocean <- data.frame( 'days'   = dBATI5$date,
                           'temp'   = dBATI5$temp_mean, 
                           'salt'   = dBATI5$salt_mean,
                           'DLI'    = dBATI5$DLI,
                           'pCO2'   = daily_swCO2$dCO2,
                           'totAlk' = daily_TA$B5 )

head(daily_ocean, 10)


#---- Part 2 Grow a kelp plant during main growing season (MAY to SEPT) ----
# Go with BATI5 mooring as that's closer to the carbon data from the Columbia

# Load growth function, with light and temp inhibitions, and a plotting function.
source( "LSSM_DEB.R" )

#Notes: Fronds have to grow faster than stipe for stipe to reach cap (i.e., max fraction)
kelp_grow <- grow_kelp4(daily_ocean )

#write.csv(kelp_grow, paste0( DEB_dir, "/kelp_grow_Nov15.csv"), row.names = FALSE)

kelp_grow$B_fronds_gDW

plot( c( 0, diff( kelp_grow$B_fronds_gDW )) )



par(mfrow=c(1,1), cex=1.2)
FullKelpPlot( kelp_grow )

# Examine growth inhibitors
matplot(kelp_grow$days[-1], kelp_grow[-1,c("fT","fL")], type="l", lwd=2,
        col=c("red","blue"), ylab="Modifier (0–1)", xlab="Date")
legend("bottomleft", legend=c("Temp limit","Light limit"),
       col=c("red","blue"), lty=1, lwd=2)


daily_kelp <- c( 0, diff( kelp_grow$B_fronds_gDW ))

matplot(
  kelp_grow$days[-1],
  cbind( kelp_grow$fT[-1], kelp_grow$fL[-1], daily_kelp[-1] ),
  type = "l", lwd = 2, col = c("red", "blue", "darkgreen"),
  xlab = "Date", ylab = "Modifier (0–1)"
)

legend(
  "bottomleft",
  legend = c("Temp limit", "Light limit", "Daily kelp growth"),
  col = c("red", "blue", "darkgreen"), lty = 1, lwd = 2
)



# --- First plot: fT and fL on left axis ---
par(mar = c(5, 5, 4, 6), cex=1.3)
matplot(
  kelp_grow$days[-1],
  kelp_grow[-1, c("fT", "fL")],
  type = "l", lwd = 2, lty = 1,
  col = c("red", "blue"),
  xlab = "Date",
  ylab = "Modifier (0–1)"
)

legend( x=19480, y=0.95,
       legend = c("Temp limit", "Light limit"),
       col = c("red", "blue"), lty = 1, lwd = 2)

# --- Add kelp growth on the right axis ---
par(new = TRUE)   # allow overplotting

plot(
  kelp_grow$days[-1],
  daily_kelp[-1],
  type = "l", lwd = 2, col = "darkgreen",
  axes = FALSE, xlab = "", ylab = ""
)

axis(side = 4, las=1)  # right y-axis
mtext("Daily kelp growth (g DW/day)", side = 4, line = 3, cex=1.3)

legend( x=19480, y=4,
       legend = "Kelp growth",
       col = "darkgreen", lty = 1, lwd = 2)





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
