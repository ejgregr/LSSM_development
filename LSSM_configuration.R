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
required.packages <- c( "ggplot2", "reshape2", "lubridate", "dplyr", "stringr",
                        "rmarkdown","knitr", "tinytex", "kableExtra", #currently causing trouble.
                        "seacarb", "gsw", "suncalc")

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
#Broughton lat and lon for light simulation
latitude  <- 50.67   #  50 deg, 40'
longitude <- -127.50 # 127 deg, 30'
bc_time <- 'America/Vancouver'

data_dir <- "C:/Data/Git/LSSM_development/Data"
DEB_dir  <- "C:/Data/Git/LSSM_development/DEB"
PDF_dir  <- "C:/Data/Git/LSSM_development/PDF" # for Markdown script

#---- Nereo growth parameters
B_init      <- 0.025  # Initial mass of sporophyte (est. at 25 mg based on ChatGPT)
B_max       <- 9.23   # Max biomass for an adult Nereo. Calculated from field results (Weigel and Pfister 2021)
r_max       <- 0.065  # Maximum daily growth rate, assuming 6 month growth and mature plant is 9.23 kg 
sp_start    <- 0.005  # Established sporophyte mass (5 g)
wet_to_dry  <- 0.13   # Water content of Nereo (Bullen et al. 2024)
dry_to_C    <- 0.25   # Carbon content of Nereo (dry) (Bullen et al. 2024)

#---- Environmental influence on growth rate
T_opt     <- 10        # Optimal temperature (°C)
T_max     <- 14        # Maximal temperature for growth (°C)
DLI_opt   <- 30        # Optimal DLI in mol/m2/day
DLI_range <- 20        # Range (+/-) that allows growth in mol/m2/day 


#---- Anticipated parameters ----
# Elemental stuff in moles, complex/diverse molecules in g. 
ounits <- list( sizeu = 'ha',
                tempu = 'C',
                saltu = 'ppt',
                DO2u  = 'mol/m3',
                DICu  = 'mol/m3',
                DOCu  = 'g/m3',
                POCu  = 'g/m3',
                NOXu  = 'g/m3',
                ph    = 'pH')


#-------------------------- Functions ----------------------------------------

#---- Utility growth functions ----
# Logistic growth, with option for temperature and light factors.
logistic_growth <- function(B0, K, r, step, temp=1, lite=1) {
  B_predicted <- K / (1 + ((K - B0) / B0) * exp(-r * step * temp * lite))
  return(B_predicted)
}

# Calculate temperature effect on growth rate
# Pontier shows growth rate relatively stable below 10C, declines to about 1/2 by 14C
t_scale <- function(T) {
  scaled <- ifelse(T <= 10, 1, 0.5 / (T / T_max * 0.7))
  return( scaled )
}

# Light-dependent growth rate function
# Pontier shows growth rate peaks ~30 DLI. Simplify their decline to be 1/2 on both sides  
DLI_scale <- function(lite) {
  #(0.75*((lite - DLI_opt)^.5) / (DLI_range^2))
  1 - (((lite - DLI_opt) / DLI_range)^2) * 0.5
}


#==== Load and prep data =====

# Uses hard-coded name for source text file. 
# Assumes seawater CO2 measurements are what we want. 
# sw CO2 measured in ppm which is euqivalent to uatm requred by carb().
# Summarizes several years of data into a 'climatology', and interpolates missing days.
LoadAkFerryCO2Data <- function(){
#  x<- read.csv(paste(data_dir, "HakaiColumbiaFerryResearch.txt", sep="/"))
  x<- read.csv(paste(data_dir, "HakaiColumbiaFerryResearch_V2.txt", sep="/"))
  
  # Keep just date, T, S, and SW_pCO2_wet_SST ... 
  xx <- x[, c(1,4,5,9)]
  names(xx)
  xx <- xx[ complete.cases( xx ), ]
  aa <- as.Date( xx$time )
  
  # Build a dataframe
  #foo       <- data.frame( "date" = aa, "SW_CO2" = b)
  #CO2_daily <- data.frame( "month"=month(foo$date), "day"=day(foo$date), "CO2"=foo$SW_CO2  )
  CO2_daily <- data.frame("month"=month(aa), "day"=day(aa),
                          "CO2"  = xx$SW_pCO2_wet_SST,
                          "Temp" = xx$TSG_T,
                          "Salt" = xx$TSG_Salinity)
  
  # Get the daily mean for available dates from Oct 2017 to Oct 2019
  CO2_daily <- CO2_daily %>%
    group_by(month, day) %>%
    summarise(
      CO2mn = mean(CO2,  na.rm = TRUE),   # Mean of SW CO2
      Tmn   = mean(Temp, na.rm = TRUE),
      Smn   = mean(Salt, na.rm = TRUE),
    )
  
  # Need to rebuild the date after summarising
  aa <- as.Date(paste("2023", CO2_daily$month, CO2_daily$day, sep='-'), format = "%Y-%m-%d")
  CO2_daily$date <- aa
  
  return(CO2_daily)
}

# Take the daily AK data make it match the BATI mooring days (i.e., extrapolate)
PrepDailyAKFerryCO2Data <- function(dCO2, all_days){
  
  # Start with some dates on the CO2 data. Needs a faux year
  # full_dates <- as.Date(paste("2023", dCO2$month, dCO2$day, sep='-'), format = "%Y-%m-%d")
  # Expand the frame to include all the days in all_days
  expanded_data <- data.frame(
    Date = all_days,
    #    dCO2 = ifelse(all_days %in% full_dates, dCO2$CO2mn[match(all_days, full_dates)], NA)
    dCO2 = ifelse(all_days %in% dCO2$date, dCO2$CO2mn[match(all_days, dCO2$date)], NA)
  )
  
  # Interpolate the NAs created above
  expanded_data$dCO2 <- interpolateC02NAs(expanded_data$dCO2)
  return(expanded_data)
}

# Interpolate missing daily values for ambient seawater pCO2 values
#   Interpolated values are the mean of the last 3, next valid values. 
interpolateC02NAs <- function(series) {
  for (i in seq_along(series)) {
    if (is.na(series[i])) {
      # Get the three previous valid values
      prev_values <- tail(series[1:(i-1)][!is.na(series[1:(i-1)])], 3)
      
      # Get the next valid value
      next_value <- head(series[(i+1):length(series)][!is.na(series[(i+1):length(series)])], 1)
      
      # Compute the average if we have enough valid values
      if (length(prev_values) == 3 && length(next_value) == 1) {
        series[i] <- mean(c(prev_values, next_value), na.rm = TRUE)
      } else {
        series[i] <- mean(c(prev_values), na.rm = TRUE)
      }
    }
  }
  return(series)
}

Load2023MooringData <- function(){
  ctd<- read.csv(paste(data_dir, "ctd_surface_cond_moorings2023.csv", sep="/"))
  
  # set data types
  ctd$date_time<- ymd_hms(ctd$date_time,  tz = "America/Vancouver")
  ctd$temperature_C<- as.numeric(ctd$temperature_C)
  ctd$salinity_psu<- as.numeric(ctd$salinity_psu)
  ctd$depth_m<- as.numeric(ctd$depth_m)
  
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

# Extract Temp and Salinty from specified BATI mooring. 
PrepBATIMooringData <- function( moor_dat ){
  ts_out <- data.frame( 
    cbind( "unixdate" = moor_dat$date_time,
           "temp"    = moor_dat$temperature_C,
           "salt"    = moor_dat$salinity_psu,
          "conduct"  = moor_dat$conduct_ms.cm )
  )
  # Now ensure unique tuples of datetime, average if necessary. 
  ts_out <- ts_out %>%
    group_by(unixdate) %>%
    summarise(
      temp = mean(temp, na.rm = TRUE),   # Mean of temp
      salt = mean(salt, na.rm = TRUE),    # Mean of salt
      conduct = mean(conduct, na.rm = TRUE)    # Mean of conductivity
    )
  # And add a readable datestamp
  ts_out$datestamp <- as.POSIXct(ts_out$unixdate, origin = "1970-01-01", tz = "America/Vancouver")
  
  return(ts_out)
}  

# Summarise mooring data to daily by ChatGPT
DailyStats <- function(df, tz = "America/Vancouver") {
  df %>%
    # convert timestamp to calendar day
    mutate(date = as.Date(format(datestamp, tz = tz, format = "%Y-%m-%d"))) %>%
    # separate daily values by dataset
    group_by(date) %>%                    
    summarise(
      across(
        where(is.numeric),
        list(mean = ~mean(.x, na.rm = TRUE),
             sd   = ~sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) %>%
    # Convert Date back to POSIXct (datestamp) at noon
    mutate(datestamp = as.POSIXct(paste0(date, " 12:00:00"), tz = tz)) %>%
    select(datestamp, everything()) %>%        # put datestamp first
    # Remove unixdate stats if present
    select(-matches("^unixdate_"))
}

# Calculate total alkalinity from salinity
CalcAlk <- function( salt ) {
  TA <- (48.7709 * salt + 606.23) # μmol kg-1 
  TA <- TA / 1e6 # mol kg-1
  return( TA )  
}

# Add surface light (PAR) to the provided data frame
AddSurfaceLight <- function(df, lat, lon, datetime_col = "datestamp", tz = "UTC") {
  
  # Ensure datetime is POSIXct
  df <- df %>%
    mutate(datetime_local = as.POSIXct(.data[[datetime_col]], tz = tz))
  
  # Get sun position
  sun <- getSunlightPosition(date = df$datetime_local, lat = lat, lon = lon)
  
  # Convert altitude to degrees
  alt_deg <- sun$altitude * 180/pi
  
  # Very simple surface PAR approximation (µmol photons m⁻² s⁻¹)
  PAR <- ifelse(sun$altitude > 0, 2.1 * 1000 * sin(sun$altitude), 0)
  
  # Add back into dataframe
  df$PAR_umol_m2_s <- PAR
  return(df)
}

# PAR → DLI (mol m^-2 day^-1), per calendar day
# df needs a POSIXct time column and a PAR column (µmol m^-2 s^-1)
# Assumes regularly spaced timestamps
ParToDLI <- function(df, datetime_col = "datestamp", par_col = "PAR_umol_m2_s", tz = "UTC") {
  # extract & order
  ts  <- as.POSIXct(df[[datetime_col]], tz = tz)
  par <- as.numeric(df[[par_col]])
  o <- order(ts); ts <- ts[o]; par <- par[o]
  
  # infer sampling step (seconds) from the median interval
  if (length(ts) < 2) stop("Need at least two timestamps")
  dt_sec <- as.numeric(stats::median(diff(ts), na.rm = TRUE), units = "secs")
  
  # sum PAR within each day, then convert: (µmol s^-1) * s → µmol; /1e6 → mol
  day   <- as.Date(ts, tz = tz)
  par_sum <- tapply(par, day, function(x) sum(x, na.rm = TRUE))
  dli <- (par_sum * dt_sec) / 1e6
  
  data.frame(Date = as.Date(names(dli)), DLI_mol_m2_day = as.numeric(dli), row.names = NULL)
}


#==== Plotting and data display ====

PlotBATI <- function(df_all,
                      var = c("temp", "salt", "conduct", "temp_mean", "salt_mean"),
                      head = "",
                      start = NULL, end = NULL) {
  
  var <- match.arg(var)
  
  # Optional date filters
  if (!is.null(start)) {
    df_all <- df_all[df_all$datestamp >= as.POSIXct(start), ]
  }
  if (!is.null(end)) {
    df_all <- df_all[df_all$datestamp <= as.POSIXct(end), ]
  }
  
  ggplot(df_all, aes(x = datestamp, y = .data[[var]], color = Mooring)) +
    geom_line(linewidth = 0.8) +
    labs(
      x = "Time",
      y = var,
      title = head,
      color = "Dataset"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}


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
# BUT conductivity values seem to be 2 orders of magnitude off?
#gsw.SP_from_C(conductivity, temperature, pressure)




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
                    ph     = 7.0 )  
  
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


# Find the indices where the target string matches elements in the string vector
monthIndex <- function(aMonth, moString) {
  index <- which(moString == aMonth)
  # If the string is not found, return NA
  if (length(index) == 0) {
    return(NA)
  } else {
    return(index)
  }
}


# Fin.


