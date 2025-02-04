#################################################################################
# Script:  LSSM_main.R
# Created: June 10 2024. EJG
# 
# Control script for the Local Seaweed Services Model. 
# Will eventually draw extensively on the regional model developed by Cam Bullen. 
#
## Updates: 
# 2025/01/06: Happy New Year. Revisiting after a few month hiatus.
# 2025/01/24: After a few weeks of development, have made some progress on 
# the parameter-based model. Reviewed Pontier, and Weigel and Pfister. And have
# conceptual design of what to do with water chemistry and dilution. 
# --> Deferring the life-cycle model until we sort out daily chemistry.
#
# TO DO: 
# - Fix day length so can remove kludges
# - Fix the q_mult fudge factor in PrepLightSim() used to match
#   calc'd DLI to reported by Pontier
# - Daily DIC drawdown;
# - Daily change in water parcel chemistry
# - 
#################################################################################

rm(list=ls(all=T)) # Erase environment.
set.seed(42)       # Setting seed for reproducibility

options(warn = -1) # CAREFUL. But the group messages are really annoying. :\

# Load packages and functions ... 
setwd( "c:/Data/Git/LSSM_development")
source( "LSSM_configuration.R" )


#---- Part 1: Create a DF of predictors to be used by the various growth models. ---- 
# Begin with peak growth phase of life cycle. 
# Time series (May to Sept) tied to data from BATI moorings, so start with those.

#---- Part 1a) Mooring data preparation ----
ctd_BATI <- Load2023MooringData()
t_stn5 <- ctd_BATI[ ctd_BATI$site == "mooring5", ]
t_stn6 <- ctd_BATI[ ctd_BATI$site == "mooring6", ]

BATI5 <- PrepBATIMooringData( t_stn5, start_date, end_date )
BATI6 <- PrepBATIMooringData( t_stn6, start_date, end_date )

length(hour_stamps)
length(BATI5$temp)

#---- Show Hourly temp and salinity for BATI moorings 5 and 6 ----

# KLUDGE: Shortening hour_stamps by 1h. 
x <- hour_stamps[-length(hour_stamps)]

par(mfrow=c(2,1), mar=c(4,4,2,1) )

plot( x, BATI5$temp, type='l', main= "BATI Mooring 5", xlab="", ylab="Temperature (C)" )
plot( x, BATI6$temp, type='l', main= "BATI Mooring 6", xlab="", ylab="Temperature (C)" )

plot( x, BATI5$salt, type='l', main= "BATI Mooring 5", xlab="", ylab="Salinity (psu)" )
plot( x, BATI6$salt, type='l', main= "BATI Mooring 6", xlab="", ylab="Salinity (psu)" )

# Calculate and show Daily mooring data ... needed for parametric model.
day_BATI5 <- MooringtoDays(BATI5)
x <- day_stamps
plot( x, day_BATI5$temp, type='l', main= "Daily Temperature (BATI 5)", xlab="", ylab="(C)" )
plot( x, day_BATI5$salt, type='l', main= "Daily Salinity (BATI 5)", xlab="", ylab="(psu)" )


#---- Part 1b) Organize light data ---- 
# Temporal resolution and extents depend on the mooring temperature data
# Light data available from sensors but needs cleaning up. Simulate for now.

hourly_PAR <- HourlyLightSim( hour_stamps, latitude, longitude )
daily_DLI <- DailyLight( day_stamps, hourly_PAR )

# (Kludges) 
daily_DLI <- daily_DLI[-c(151,152)] # Drop last 2 elements to match lengths 
daily_DLI[ length(daily_DLI) ] <- 17 # Fist last day cuz incomplete so biased.

#---- Show hourly PAR and  DLI ----
length(hourly_PAR)
plot( hour_stamps, hourly_PAR, type='l', main= "Hourly Photosynthetically Active Radiation (simulated)", 
      xlab="", ylab="µmol photons m⁻² s⁻²") 

# KLUDGE: Shortening timestamps by 1. 
plot( day_stamps, daily_DLI, type='l', main= "Daily Light Intensity (simulated)", 
      xlab="", ylab="mol photons / m² d" )

#---- Part 1c) Get alkalinity using BATI Salt time series ----
# Using equation from Evans et al. 2015:  TA = 48.7709*S + 606.23 (μmol kg-1) 

daily_TA <- (48.7709 * day_BATI5$salt + 606.23) # μmol kg-1 
daily_TA <- daily_TA / 1000 # mol kg-1

plot( day_stamps, daily_TA, type='l', main= "Average, Daily Total Alkalinity", 
      xlab="", ylab="mol / kg" )

# Empty plot to fill the space 
plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE, main = "Empty Plot")


#---- Part 1d) Estimate ambient pCO2 from Ak Ferry data ---- 
# Not strictly related to kelp growth, but necessary to estimate pH effect

ak_dat      <- LoadAkFerryCO2Data()
daily_swCO2 <- PrepAKFerryCO2Data(ak_dat, day_stamps)
# convert ppm to uatm (microatmospheres) at std atm pressure for carb()
daily_swCO2 <- daily_swCO2 * 1.01325 

# Show the Ak ferry data, before and after interpolation  
plot(ak_dat$date, ak_dat$CO2mn, type = "o", main = "Ak Ferry daily Seawater pCO2 climatology",
     xlab = "Date", ylab = "uatm", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = ak_dat$date, format = "%m-%d")

plot(daily_swCO2$Date, daily_swCO2$dCO2, type = "o", main = "Ak Ferry daily seawater pCO2 climatology (interpolated)",
     xlab = "Date", ylab = "uatm", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = daily_swCO2$Date, format = "%m-%d")


#----- Part 2 Grow a kelp plant during main growing season (MAY to SEPT) -----
# Parametric model, currently using simple logistic growth, with 
# optional temperature and light inhibition factors

# Calculate temperature and light variability (inhibition factors)
temp_fact <- t_scale( day_BATI5$temp )
DLI_fact  <- DLI_scale( daily_DLI )

# Show the inhibition factors 
plot( temp_fact, type = 'l' )
plot( DLI_fact, type = 'l' )


# Straight up logistic growth, no inhibition factors
y <- NULL
for (t in 1:length( temp_fact )) {
  one <- logistic_growth( B_init, B_max, r_max, t)
  y <- c(y, one)
}
log_simp <- y

# Logistic growth with temperature and light inhibition factors
y <- NULL
for (t in 1:length( temp_fact ) ) {
  one <- logistic_growth( B_init, B_max, r_max, t, temp_fact[t], DLI_fact[t])
  y <- c(y, one)
}
log_env <- y

# Comparison of simple vs. inhibited plant growth
plot( day_stamps,log_simp*1000, type='l',xlab = 'Days', ylab = 'grams', main='A plant - simple logistic growth')
plot( day_stamps, log_env*1000, type='l',xlab = 'Days', ylab = 'grams', main='A plant - inhibited logistic growth')

# NOTE: log_env implies 'shrinkage' rather that slowed growth. This is at least 
# partially because a population would shrink if r is reduced, but doesn't work
# exactly the same with growth. Needs a think. 

# --> Work with un-inhibited growth for now.

#----- Part 3 - Chemistry of the plant growth model -----
# NOTE seacarb() requires everything in mol/kg, and uses flags to identify the 
# two inputs from which to calculate carbonate chemistry. Of interest here:

flag_CA <- 4  # CO2 and Alkalinity 
flag_CD <- 25 # CO2 and DIC  
flag_AD <- 15 # Alkalinity and DIC 

# MOLES OF DIC FIXED BY KELP ... 
# Starting with the daily kelp biomass, get daily grams DIC fixed.   
gDIC_fixed <- log_simp * 1000 * wet_to_dry * dry_to_C  

#Using mol weights (g/mol) of structure and C reserve from DEB, estimate daily mols C
#Mol wt of structure = 27.51; mol wt of C reserve = 30. Use average

molDIC_fixed <- gDIC_fixed / (27.51+30) / 2

#Now the cumulative DIC fixed over 150 days (is the same shape as growth)
plot( day_stamps, molDIC_fixed )

#Use diff() to return difference btwn consecutive elements. So now its DIC fixed/day
delkDIC <- diff( molDIC_fixed ) 
plot( day_stamps[-length(day_stamps)], delkDIC, ylab = "mol DIC fixed / day" )
#This is equivalent to the rate of kelp growth. 

# AMBIENT DIC  ...
# Illustrate with 3 plots
par(mfrow=c(3,1) )

# Convert  ambient pCO2 to DIC ... Needs Alkalinity!
par(mar=c(0,4,3,1) )
plot( daily_TA, type='l', xlab='', xaxt = "n", main="Ambient pCO2 to DIC" )
par(mar=c(0,4,1,1) )
plot( daily_swCO2, type='l', xlab='',xaxt = "n" )
# Note the dropping of the date column, and the need to drop 2 dates
x <- daily_swCO2[-dim(daily_swCO2)[[1]],2 ]
x <- x[-length(x)]
y <- daily_TA[-length( daily_TA )]

# carb() wants pCO2 in uatm. Ferry data is in 

amb_DIC <- carb( flag_CA, x, y )$DIC
par(mar=c(3,4,1,1) )
plot( day_stamps[-length(day_stamps)], amb_DIC, type='l', xlab='' )

# We can look at daily changes in pH from ambient attributable to kelp. 
# This effect is proportional to the size of the plant and its rate of growth. 
# (I don't think this is all captured with the logistic growth model)
length(amb_ph)
amb_ph  <- carb( flag_AD, y, amb_DIC )$pH
kmod_ph <- carb( flag_AD, y, amb_DIC-delkDIC )$pH
plot( day_stamps[-length(day_stamps)], xlab="",
      amb_ph, type='l', col='red', main='Daily change in pH in 1 m2 of water')
lines( day_stamps[-length(day_stamps)],
       kmod_ph, type='l', col='green' )
legend("bottomleft", legend = c("Ambient", "Kelp-affected"), col = c("red", "green"), lwd = 2)




# Above is all based on 1 kg of water. 
# At what mass/volume is the reduction in DIC not significant?


# Example usage:
DIC_removed <- 156.1458  # Total DIC removed by kelp (moles)
DIC_ambient <- 2.572458  # Mean ambient DIC (moles/kg)
dilution_mass(DIC_removed, DIC_ambient, epsilon = 0.05)  # 5% threshold







x <- seq(10, round( dilution_mass( sum(molDIC_fixed), mean(amb_DIC), 0.01 )), by=100)
yK <- DIC_removed / x
yA <- 
z <- mean( daily_TA)
csystem <- carb( flag_AD, z, y )


# Function to estimate the 
dilution_mass <- function(DIC_removed, DIC_ambient, epsilon = 0.01) {
  # Ensure epsilon is positive and reasonable
  if (epsilon <= 0 || epsilon > 1) {
    stop("Epsilon should be between 0 and 1 (e.g., 0.01 for 1%, 0.05 for 5% significance)")
  }
  # Calculate required water mass (kg)
  M <- DIC_removed / (epsilon * DIC_ambient)
  
  return(M)
}



#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 

today <- format(Sys.Date(), "%Y-%m-%d")

rmarkdown::render( "LSSM_documentation.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = DEB_dir,
                   output_file = paste0( "LSSM_DEB_testing_", today ))

#FIN
