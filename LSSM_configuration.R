#----------------------------------------------------------------------------
# Script:  LSSM_configuration.R
# Created: September 2024
#
# Purpose: Create and load initial data structures for the LSSM. 
#
# Notes:
#  - 

#================================== Load require packages =================================
# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "lubridate", "dplyr",
                        "rmarkdown","knitr", "tinytex", "kableExtra",
                        "seacarb", "gsw", "truncnorm")

# Other packages that might be useful. 
# "tidyr", "raster", "stringr", "rasterVis",
# "RColorBrewer", "factoextra", "ggpubr", "cluster", 
# "diffeR", "vegan", "ranger", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
# "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr", "binr", 'gwxtab'

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)
version$version.string


#==== Configuration ====

moText <- c( "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Elemental stuff in moles, complex/diverse molecules in g. 
ounits <- list( sizeu = 'ha',
                tempu = 'C',
                saltu = 'ppt',
                DO2u  = 'mol/m3',
                DICu  = 'mol/m3',
                DOCu  = 'g/m3',
                POCu  = 'g/m3',
                NOXu  = 'g/m3',
                ph    = 'pH'    )

# Constants
latitude  <- 49.2827 # Vancouver's latitude
longitude <- -123.1207 # Vancouver's longitude

data_dir <- "C:/Data/Git/LSSM/Data"
DEB_dir  <- "C:/Data/Git/LSSM/DEB"

DEB <- 



#==== Functions ====

#==== Load and prep data ====
Load2023MooringData <- function(){
  ctd<- read.csv(paste(data_dir, "ctd_surface_cond_moorings2023.csv", sep="/"))
  
  # set data types
  ctd$date_time<- ymd_hms(ctd$date_time,  tz = "America/Vancouver")
  ctd$temperature_C<- as.numeric(ctd$temperature_C)
  ctd$salinity_psu<- as.numeric(ctd$salinity_psu)
  ctd$depth_m<- as.numeric(ctd$depth_m)
  
  # Add year month and day in separated columns
  ctd$year<- year(ctd$date_time)
  ctd$month<- month(ctd$date_time)
  ctd$day<- day(ctd$date_time)
  ctd$dmy = dmy(paste(ctd$day, ctd$month, ctd$year, sep="-"))
  colnames(ctd)[13]<- "ymd"
  
  # classify moorings based on environmental clusters to merge and plot 
  # ctd$cluster<- NA
  # 
  # ctd[ctd$site == "mooring1", "cluster"]<- "5"
  # ctd[ctd$site == "mooring2", "cluster"]<- "5"
  # ctd[ctd$site == "mooring3", "cluster"]<- "5"
  # ctd[ctd$site == "mooring4", "cluster"]<- "2"
  # ctd[ctd$site == "mooring5", "cluster"]<- "2"
  # ctd[ctd$site == "mooring6", "cluster"]<- "4"
  # ctd[ctd$site == "mooring7", "cluster"]<- "4"
  # ctd[ctd$site == "mooring8", "cluster"]<- "5"
  
  return(ctd)
}

# Create a DF that can be passed to the growth model.
# Resolution and extents depend on the temperature data.
# Currently dependent on BATI sensor measures processed by Barbosa.
PrepSensorData <- function( moor_dat, smo, emo){
  
  # First pull hourly temp data for specified start/end months
  m_idx <- month(moor_dat$date_time) >= smo & month(moor_dat$date_time) <= emo
  ts_out <- data.frame( 
    cbind( "month" = month( moor_dat[ m_idx, "date_time" ]),
           "day"   = day( moor_dat[ m_idx, "date_time" ]),
           "hour"  = hour( moor_dat[ m_idx, "date_time" ]),
           "temp"  = moor_dat[ m_idx, "temperature_C" ],
           "salt"  = moor_dat[ m_idx, "salinity_psu" ])
  )
  
  # Now ensure unique tuples of month, day, hour, average if necessary. 
  ts_out <- ts_out %>%
    group_by(month, day, hour) %>%
    summarise(
      temp = mean(temp, na.rm = TRUE),   # Mean of temp
      salt = mean(salt, na.rm = TRUE)    # Mean of salt
    )
  
  return(ts_out)
}  


#==== Simulating insolation ====
# Function to calculate the solar declination angle (in degrees) based on the day of the year
solar_declination <- function(day_of_year) {
  23.44 * sin((360 / 365) * (day_of_year - 81) * pi / 180)
}

# Function to calculate the hour angle (in degrees) based on the local solar time
hour_angle <- function(time) {
  # Local Solar Time (LST) is approximated here as the hour in UTC + longitude / 15
  local_time <- hour(time) + (minute(time) / 60)
  LST <- local_time + longitude / 15
  15 * (LST - 12) # Convert to degrees
}

# Function to calculate the solar elevation angle (in degrees)
solar_elevation_angle <- function(time, lat, lon) {
  day_of_year <- yday(time) # Day of the year
  declination <- solar_declination(day_of_year)
  hour_ang <- hour_angle(time)
  
# Convert to radians for trigonometric calculations
  declination_rad <- declination * pi / 180
  latitude_rad <- lat * pi / 180
  hour_ang_rad <- hour_ang * pi / 180
  
# Calculate the solar elevation angle in radians
  sin_alpha <- sin(latitude_rad) * sin(declination_rad) +
    cos(latitude_rad) * cos(declination_rad) * cos(hour_ang_rad)
  
# Convert back to degrees
  solar_angle <- asin(sin_alpha) * 180 / pi
  
# And replace negative values with 0 (i.e., just dark)
  solar_angle[ solar_angle < 0] <- 0

  return( solar_angle ) # Return the angle, ensuring it's non-negative
#return(max(solar_angle, 0)) # Return the angle, ensuring it's non-negative
}

# Function to calculate PAR (Photon Flux Density in mol photons/m²/s) based on the solar elevation angle
CalculatePhotons <- function(solar_angle) {
  # Assume a clear-sky model for simplicity:
#  if (solar_angle <= 0) return(0) # No sunlight when the sun is below the horizon
  
  # Empirical estimate of solar irradiance based on solar angle (in W/m²)
  I <- 1361 * sin(solar_angle * pi / 180) # 1361 W/m² is the solar constant
  
  # Approximate PAR as 45% of total solar irradiance
  I_PAR <- 0.45 * I
  
  # Convert to Photon Flux Density (mol photons/m²/s)
  # Using a conversion factor (4.57e-6 mol/W/s) based on average photon energy
  photon_flux_density <- I_PAR * 4.57e-6 * 3600
  return(photon_flux_density)
}


#==== Plotting and data display ====

PlotInputs <- function( indat, ptitle ){
  colnames( indat ) <- c( "temp", "nutrients", "DIC", "light")
  see_in <- data.frame(
              cbind( "hour" = 1:dim(indat)[[1]], indat))
  
  melt_in <- melt( see_in, id.vars ="hour" )

  ggplot(melt_in) +
    geom_point( aes(x=hour, y=value), size=0.4, color="darkgreen") +
    facet_wrap(~variable, ncol = 2, scales = "free_y") +
    labs(x = "Cumulative hours", title = ptitle ) +
 #   theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}

PlotAllDEBResults <- function( deb_out ){
  
  par(mfrow=c(3,3), mar=c(1,4,3,2))
  for (i in colnames(deb_out)[3:8] ){
    plot( deb_out[, i]~deb_out$days, type='l', xlab="", xaxt = "n", ylab=i ) }
  par(mar=c(4,4,0,2))
  for (i in colnames(deb_out)[9:11] ){
    plot( deb_out[, i]~deb_out$days, type='l', xlab="Days", ylab=i ) }
  
  par(mar=c(1,4,3,2))
  for (i in colnames(deb_out)[12:17] ){
    plot( deb_out[, i]~deb_out$days, type='l', xlab="", xaxt = "n", ylab=i ) }
  par(mar=c(4,4,0,2))
  for (i in colnames(deb_out)[18:20] ){
    plot( deb_out[, i]~deb_out$days, type='l', xlab="Days", ylab=i ) }
  
  par(mar=c(1,4,3,2))
  for (i in colnames(deb_out)[21:23] ){
    plot( deb_out[, i]~deb_out$days, type='l', xlab="", xaxt = "n", ylab=i ) }
  par(mar=c(4,4,0,2))
  for (i in colnames(deb_out)[24:26] ){
    plot( deb_out[, i]~deb_out$days, type='l', xlab="Days", ylab=i ) }
}


# A function (from 'gsw') to convert conductivity to psu, note dependence on T and P.
#   salinity = gsw.SP_from_C(conductivity, temperature, pressure)
# BUT conductivity values seem to be 2 orders of magnitude off?


#--------------  STUB FUNCTIONS ---------------

# Input either a shape file name or the shape file itself

InitializeClusters <- function( mapdat ){
  
  aCluster <- list( cname = "", 
                    size = 10.0,
                    temp = 15.0,
                    salt = 30.0,
                    NOX  = 99.0,
                    DO2  = 99.0,
                    DIC  = 99.0,
                    DOC  = 99.0,
                    POC  = 99.0,
                    DCO2 = 1.0,
                    carbA  = 1.0, 
                    carb   = 1.0,
                    bicarb = 1.0,
                    ph     = 7.0
  )  
  
  z <- vector( "list", 6)
  
  for (i in 1:length(mapdat)) {
    j <- aCluster
    j$cname <- paste0( "cluster ", i )
    z[[i]] <- j
  }
  
  return(z)
}


# Create necessary GLOBAL data strUctures and populate with initial state.
# These structures will be lists (1 to n) of a list of attributes.
InitializeSimulation <- function( nMonths, firstMo, clusts ){
  
  aState <- list( month = "", 
                  clusters = clusts
  )
  
  oStates <- vector("list", nMonths)
  
  moIdx <- monthIndex( firstMo, moText )
  
  for (i in 1:nMonths) {
    oStates[[i]] <- aState
    oStates[[i]]$month <- moText[moIdx]
    moIdx <- moIdx+1
  }
  
  return( oStates )
}



monthIndex <- function(aMonth, moString) {
  # Find the indices where the target string matches elements in the string vector
  index <- which(moString == aMonth)
  # If the string is not found, return NA
  if (length(index) == 0) {
    return(NA)
  } else {
    return(index)
  }
}


# FIN

# FIN


# Fin.