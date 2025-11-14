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
# 2025/11/06: Considerably shortened as now does only data loading. Rest moved to main. 
#################################################################################

#---- Part 1) Mooring data preparation ----
# Prototype model will be tied to temporal extents of mooring data ((May 4 to Oct 11 = 160 days)
# NOTE: Correct application of tz is critical.

start_date <- "2023-05-04 00:00:00 PDT"
end_date   <- "2023-10-11 00:00:00 PDT"

ctd_BATI <- Load2023MooringData()
t_stn5 <- ctd_BATI[ ctd_BATI$site == "mooring5", ]
t_stn6 <- ctd_BATI[ ctd_BATI$site == "mooring6", ]
t_stn7 <- ctd_BATI[ ctd_BATI$site == "mooring7", ]

# Pull and rename columns, ensure hourly data 
BATI5 <- PrepBATIMooringData( t_stn5 )
BATI6 <- PrepBATIMooringData( t_stn6 )
BATI7 <- PrepBATIMooringData( t_stn7 )

# Trim mooring data to date range 
idx <- ( BATI5$datestamp  >= start_date ) & 
  ( BATI5$datestamp  <  end_date )

BATI5 <- BATI5[ idx, ]
BATI6 <- BATI6[ idx, ]
BATI7 <- BATI7[ idx, ]

# Calculate Daily mooring data ... needed for parametric model.
# Create dfs for daily conditions with mean and sd
dBATI5 <- DailyStats( BATI5 )
dBATI7 <- DailyStats( BATI7 )

#----- Part 2 Estimate ambient pCO2 from Ak Ferry data ----
# Not strictly related to kelp growth, but necessary to estimate pH effect
ak_dat      <- LoadAkFerryCO2Data()
daily_swCO2 <- PrepDailyAKFerryCO2Data(ak_dat, dBATI5$date)


# Fin.