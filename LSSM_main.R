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
# - Fix day length so can remove kludge
# - Fix the q_mult fudge factor in PrepLightSim() used to match
#   calc'd DLI to reported by Pontier
# - Daily DIC drawdown;
# - Daily change in water parcel chemistry
# - 
#################################################################################

rm(list=ls(all=T)) # Erase environment.
set.seed(42)       # Setting seed for reproducibility

# Load packages and functions ... 
setwd( "c:/Data/Git/LSSM_development")
source( "LSSM_configuration.R" )


#---- Mooring data preparation ----
ctd_BATI <- Load2023MooringData()
t_stn5 <- ctd_BATI[ ctd_BATI$site == "mooring5", ]
t_stn6 <- ctd_BATI[ ctd_BATI$site == "mooring6", ]

BATI5 <- PrepBATIMooringData( t_stn5, start_date, end_date )
BATI6 <- PrepBATIMooringData( t_stn6, start_date, end_date )

length(timestamps)
length(BATI5$temp)

par(mfrow=c(2,1), mar=c(4,4,2,1) )
plot( timestamps[-length(timestamps)], BATI5$temp, type='l', main= "BATI Mooring 5", xlab="", ylab="Temperature (C)" )
plot( timestamps[-length(timestamps)], BATI6$temp, type='l', main= "BATI Mooring 6", xlab="", ylab="Temperature (C)" )

plot( timestamps[-length(timestamps)], BATI5$salt, type='l', main= "BATI Mooring 5", xlab="", ylab="Salinity (psu)" )
plot( timestamps[-length(timestamps)], BATI6$salt, type='l', main= "BATI Mooring 6", xlab="", ylab="Salinity (psu)" )


#---- Create a DF that can be used by the various growth models. ---- 
# There are 3 phases of growth for a life cycle.
# CURRENTLY, focus on peak growth. 
#   temporal resolution and extents depend on the mooring temperature data
#   Light data also available but needs cleaning up. 

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
t_daily <- y$temp
s_daily <- y$salt

length(day_stamps)
length(t_daily)
#---- Show t_daily ----
plot(day_stamps[-length(day_stamps)], t_daily, type='l', xlab = 'Days', ylab='Temperature')

#---- Organize light data ---- 
light_DLI <- PrepLightSim( timestamps, latitude, longitude )

# (Kludges) 
light_DLI <- light_DLI[-c(151,152)] # Drop last 2 elements to match lengths 
light_DLI[ length(light_DLI) ] <- 17 # Fist last day cuz incomplete so biased.

#---- Show light_DLI ----
plot( day_stamps[-length(day_stamps)], light_DLI, type='l' )  


#----- Grow a kelp plant during PRIMARY growing season (MAY to SEPT) -----

#----- Part 1 - Using simple logistic growth -----
# logistic growth function has optional temperature and light inhibition factors

# Calculate temperature and light variability (inhibition factors)
temp_fact <- t_scale( t_daily )
DLI_fact  <- DLI_scale( light_DLI )

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
par(mfcol=c(1,2))
plot( day_stamps[-length(day_stamps)],log_simp*1000, type='l',xlab = 'Days', ylab = 'grams', main='A plant - simple logistic growth')
plot( day_stamps[-length(day_stamps)], log_env*1000, type='l',xlab = 'Days', ylab = 'grams', main='A plant - inhibited logistic growth')


#----- Part 2 - Chemistry of the plant growth model -----
# NOTE seacarb() requires everything in mol/kg
# Starting with daily plant biomass, get daily grams DIC fixed.   
gDIC_fixed <- log_env * 1000 * wet_to_dry * dry_to_C  

# Using mol weights (g/mol) of structure and C reserve from DEB, estimate daily mols C
# Mol wt of structure = 27.51; mol wt of C reserve = 30. Use average

molDIC_fixed <- gDIC_fixed / (27.51+30) / 2

# Get alkalinity for BATI time series
# Using equation from Evans et al. 2015:  TA = 48.7709*S + 606.23 (μmol kg-1) 

TA_daily <- (48.7709*s_daily + 606.23) # μmol kg-1 
TA_daily <- TA_daily / 1000 # mol kg-1



carb( )


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
