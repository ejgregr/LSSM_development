#################################################################################
# Script:  LSSM_drivers.R
# Created: January 21, 2024. EJG
# Preparation and wrangling of predictors (Temperature, light, nutrients) at 
# different temporal resolutions and annual extents to drive Nereo growth over time. 
# Descriptions of source data and model design things in the RMD file.
# Assumptions:
# 
# ## Updates: 
# 2025/xx/xx: 
#################################################################################










#------ Extra stuff ----

# Aggregate days for all months of the year to calculate monthly daylight hours
# Relies on daily day light hours (from getAnnualDLI)
# NEEDED to groundtruth Weigel and Pfister's growth.
getMonthlyHours <- function( day_lite ){
  # Get ready ... 
  days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  month_idx <- rep(1:12, times = days_in_month)
  
  # Aggregate daily values by month
  month_lite <- tapply(day_lite, month_idx, sum)  # Sum for each month
  
  return( month_lite)
}
