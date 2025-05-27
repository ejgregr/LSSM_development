#################################################################################
# Script:  LSSM_drivers.R
# Created: January 21, 2024. EJG
# Preparation and wrangling of predictors (Temperature, light, nutrients) at 
# different temporal resolutions and annual extents to drive Nereo growth over time. 
# GOAL: Create a DF of predictors to be used by the various growth models
# Descriptions of source data and model design things in the RMD file.
# Assumptions:
# Focus on sporophyte to adult growth phase of life cycle. 
#
# Updates: 
# 2025/05/27: Add work to date on loading and displaying available data.
#################################################################################

# NB: hour_stamps and day_stamps created in Configuration script.

#---- Part 1) Mooring data preparation ----
# Time series (May to Sept = 150 days) tied to data from BATI moorings, so start with those.
ctd_BATI <- Load2023MooringData()
t_stn5 <- ctd_BATI[ ctd_BATI$site == "mooring5", ]
t_stn6 <- ctd_BATI[ ctd_BATI$site == "mooring6", ]
t_stn7 <- ctd_BATI[ ctd_BATI$site == "mooring7", ]

BATI5 <- PrepBATIMooringData( t_stn5, start_date, end_date )
BATI6 <- PrepBATIMooringData( t_stn6, start_date, end_date )
BATI7 <- PrepBATIMooringData( t_stn7, start_date, end_date )

length(hour_stamps)
length(BATI5$temp)

#---- Show Hourly temp and salinity for BATI moorings 5 and 7 ----

par(mfrow=c(2,1), mar=c(4,4,2,1) )

plot( hour_stamps, BATI5$temp, type='l', main= "BATI Mooring 5", xlab="", ylab="Temperature (C)" )
plot( hour_stamps, BATI7$temp, type='l', main= "BATI Mooring 7", xlab="", ylab="Temperature (C)" )

plot( hour_stamps, BATI5$salt, type='l', main= "BATI Mooring 5", xlab="", ylab="Salinity (psu)" )
plot( hour_stamps, BATI7$salt, type='l', main= "BATI Mooring 7", xlab="", ylab="Salinity (psu)" )

# Calculate and show Daily mooring data ... needed for parametric model.
day_BATI5 <- MooringtoDays(BATI5)
day_BATI7 <- MooringtoDays(BATI7)

par(mfrow=c(2,2), mar=c(1,4,2,1) )
plot( day_stamps, day_BATI5$temp, type='l', main= "Daily T/S (BATI 5)", xaxt = "n", xlab="", ylab="Temperature (C)" )
plot( day_stamps, day_BATI7$temp, type='l', main= "Daily T/S (BATI 7)", xaxt = "n", xlab="", ylab="" )
par(mar=c(4,4,0,1) )
plot( day_stamps, day_BATI5$salt, type='l', main= "", xlab="", ylab="Salinity (PSU)" )
plot( day_stamps, day_BATI7$salt, type='l', main= "", xlab="", ylab="" )

#---- Part 2) Simulate and Organise light data ---- 
# Temporal resolution and extents depend on the mooring temperature data
# Light data available from sensors but needs cleaning up. Simulate for now.

hourly_PAR <- HourlyLightSim( hour_stamps, latitude, longitude )
daily_DLI <- DailyLight( day_stamps, hourly_PAR )
length(hourly_PAR)

# (Kludges) 
daily_DLI <- daily_DLI[-c(151,152)] # Drop last 2 elements to match lengths 
daily_DLI[ length(daily_DLI) ] <- 17 # Fist last day cuz incomplete so biased.
length(daily_DLI)

#---- Show hourly PAR and DLI ----
par(mfrow=c(2,1), mar=c(2,4,2,1) )

plot( hour_stamps, hourly_PAR, type='l', main= "Hourly Photosynthetically Active Radiation (simulated)", 
      xlab="", ylab="µmol photons m⁻² s⁻²") 

plot( day_stamps, daily_DLI, type='l', main= "Daily Light Intensity (simulated)", 
      xlab="", ylab="mol photons / m² d" )

#---- Part 3) Get alkalinity using BATI Salt time series ----
#Using equation from Evans et al. 2015:  TA = 48.7709*S + 606.23 (μmol kg-1)

daily_TA <- data.frame( "B5" = CalcAlk( day_BATI5$salt ),
                        "B7" = CalcAlk( day_BATI7$salt ))

plot( day_stamps, daily_TA$B5, type='l', main= "Average, Daily Total Alkalinity (BATI5)",
      xlab="", ylab="mol / kg" )

plot( day_stamps, daily_TA$B7, type='l', main= "Average, Daily Total Alkalinity (BATI7)",
      xlab="", ylab="mol / kg" )

#----- Part 4) Estimate ambient pCO2 from Ak Ferry data ----
# Not strictly related to kelp growth, but necessary to estimate pH effect

ak_dat      <- LoadAkFerryCO2Data()
daily_swCO2 <- PrepDailyAKFerryCO2Data(ak_dat, day_stamps)
length(daily_swCO2$Date)

# Show the Ak ferry data, before and after interpolation
plot(ak_dat$date, ak_dat$CO2mn, type = "o", main = "Full Ak Ferry daily seawater pCO2 climatology",
     xlab = "Date", ylab = "uatm", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = ak_dat$date, format = "%m-%d")

plot(daily_swCO2$Date, daily_swCO2$dCO2, type = "o", main = "Summer Ak Ferry daily pCO2 - interpolated",
     xlab = "Date", ylab = "uatm", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = daily_swCO2$Date, format = "%m-%d")

#----- Part 5) Data preparation for kelp growth and ocean chemistry ----


# All carb() related values should be in carb() units (mol/kg) and day length:
daily_ocean <- data.frame( 'days' = day_stamps,
                           'temp' = day_BATI5$temp, 
                           'salt' = day_BATI5$salt,
                           'pCO2' = daily_swCO2$dCO2,
                           'totA' = daily_TA$B5 )

# Fin.