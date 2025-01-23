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

# Create a DF that can be passed to the growth model.
# Resolution and extents depend on the temperature data.
# Currently dependent on BATI sensor measures processed by Barbosa.


#---- Organize temperature data ---- 
x <- PrepBATIMooringData( t_stn5, start_date, end_date)
plot(x$temp, type='l' )
plot(x$temp, type='l', xlab = 'Hours', ylab='Temperature')

y <- x %>%
  group_by(month, day) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),   # Mean of temp
    salt = mean(salt, na.rm = TRUE)    # Mean of salt
  )
plot(y$temp, type='l', xlab = 'Days', ylab='Temperature')

dim(y)
t_daily <- y$temp
#---- Show t_daily ----
plot( t_daily, type='l' )

#---- Organize light data ---- 

timestamps <- seq(from = start_date, to = end_date, by = "hour")

# Hourly light levels for analytic period
light_PAR <- CalculatePhotons( 
  solar_elevation_angle( timestamps, latitude, longitude )
)

plot( timestamps, light_PAR, type='l' )  
par( new=T )
length(light_PAR)
365*24
#May 1 = 121, Sep 30 = 273
lines( timestamps[(121:273)*24 ],light_PAR[(121:273)*24 ], type='l', col = 'Red' )  


# Aggregate hours to days
paste( "Working with", length(light_PAR) /24, "days of light data ...")

par_daily <- split(light_PAR, ceiling(seq_along(light_PAR) / 24))

# initialize target ... 
light_DLI <- numeric(length(par_daily))

# Calculate average PAR for each day
for (i in seq_along(par_daily)) {
  ave_par <- mean(par_daily[[i]])
  q_mult <- 1.75
  light_DLI[i] <- (ave_par * q_mult * 3600 *24) / 1e6
}

#---- Show light_DLI ----
plot( light_DLI, type='l' )  
# Do some date math to bring the date thru for presentation.



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
