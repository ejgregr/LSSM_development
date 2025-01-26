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

rm('CO2')

x<- read.csv(paste(data_dir, "HakaiColumbiaFerryResearch.txt", sep="/"))

a <- ifelse(nchar(x$s.PC_Date) < 6, paste0("0", x$s.PC_Date), x$s.PC_Date)
aa <- as.Date( a, format = "%d%m%y" )
b <- x$s.calibrated_SW_xCO2_dry

foo  <- data.frame( "date" = aa, "SW_CO2" = b)

CO2_daily <- data.frame( "month" = month(foo$date), "day"= day(foo$date), "CO2" = foo$SW_CO2  )

# Get the daily mean for available dates from Oct 2017 to Oct 2019
CO2_daily <- CO2_daily %>%
  group_by(month, day) %>%
  summarise(
    CO2mn = mean(CO2, na.rm = TRUE),   # Mean of SW CO2
  )

CO2_daily <-cbind(CO2_daily, "date" = paste(CO2_daily$month, CO2_daily$day, sep='-') )

# Create s string with faux year to plot dates on x-axis
full_dates <- as.Date(paste0("2023-", CO2_daily$date), format = "%Y-%m-%d")
# Interpolate missing days 
CO2_int <- interpolateC02NAs( CO2_daily$CO2mn )

plot(full_dates, CO2_int, type = "o", main = "Ak Ferry daily CO2 time series (interpolated)",
     xlab = "Date", ylab = "Measure", col = "blue", lwd = 2, xaxt = "n")
axis.Date(1, at = full_dates, format = "%m-%d")






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
