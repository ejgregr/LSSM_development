#----------------------------------------------------------------------------
# Script:  LSSM_config2.R
# Created: January 13 2026. 
# Purpose: Support second iteration of the nereo growth and chemistry work
#
# Notes:
#  - 
#================================== Load require packages =================================
# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "lubridate", "dplyr", "stringr",
                        "rmarkdown","knitr", "tinytex", "kableExtra", #currently causing trouble.
                        "seacarb", "gsw")

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
bc_time <- 'America/Vancouver'

data_dir <- "C:/Data/Git/LSSM_development/Data"
DEB_dir  <- "C:/Data/Git/LSSM_water_analysis/Results"
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

