#################################################################################
# Script:  LSSM_main.R
# Created: June 10 2024. EJG
# 
# Control script for the Local Seaweed Services Model. 
# Draws extensively on the regional developed by Cam Bullen. 
#
## Updates: 
# 2025/01/06: Happy New Year. Revisiting after a few month hiatus.

# TO DO: 
#  


usethis::create_from_github(
  "https://github.com/ejgregr/LSSM_development.git",
  destdir = "c:/data/git/LSSM_development/"
)

https://github.com/ejgregr/LSSM_development.git


#################################################################################

rm(list=ls(all=T))  # Erase environment.
set.seed(42)  # Setting seed for reproducibility

# Load packages and functions ... 
setwd( "c:/Data/Git/LSSM")
source( "LSSM_configuration.R" )

#==== Evolving high level code for main model ====

#==== Initialization
# mapData <- vector("integer", 6)
# 
# clusters <- InitializeClusters( mapData )
# InitializeSimulation( runLength, startMonth, clusters )

# Process interactions at monthly time step
# for (m in 1:runLength) {
#   
#   for (c in 1:clusters) {
#     
#     kelpStates( n+1 ) <- growKelp( kelpStates(n),  oceanStates(n), ecoStates(n) )
#     ecoStates( n+1 )  <- ecoMagic( kelpStates(n+1),  )
#     
#     UpdateOceanStates( kelpStates, ecoStates )
#   }
# }
# 
# 
# oceanStates <- InitializeSimulation( runLength, startMonth, clusters)
# 
# for (i in 1:length(oceanStates)) {
#   cat(oceanStates[[ i ]]$month, " ")
# }





#==== WORKING ====
# Current DEB version (rundeb_KELP.exe - 2024-09-21) uses 1 hour time steps, and runs for 1 year.
# Also has input and output files hard-coded, and expected in the same directory. 

#==== Example runs of the DEB growth model ====

# First example input from Laura ...  
#deb_input <- read.table( paste0(DEB_dir, "/original_env_conditions.txt"), header = FALSE)
#PlotInputs( deb_input, "Original test data from Laura")

#==== Create input data. ====
# Time period: May thru September (5 months) - determined by the 2023 BATI temp data. 
# Temperature taken from 2023 BATI mooring data compiled by Romina. 
#   Mooring 5 is 'reference' and mooring 6 is 'inside' the archipelago
# Nutrients and DIC are simulated as constant with variability
# Light is simulated as non-varying daily sinusoidal

ctd_BATI <- Load2023MooringData()
t_stn5 <- ctd_BATI[ ctd_BATI$site == "mooring5", ]
t_stn6 <- ctd_BATI[ ctd_BATI$site == "mooring6", ]

BATI5 <- PrepSensorData( t_stn5, 5, 9 )
head(BATI5)
par(mfrow=c(2,1), mar=c(4,4,1,1) )
plot( BATI5$temp, type='l', xlab="",  xaxt = "n", ylab="Temperature (C)" )
plot( BATI5$salt, type='l', xlab="Cumulative hours", ylab="Salinity (psu)") 
#    axis( 1, at = seq(-2, 2, by = 0.5), labels = BATI5$month) )

#==== Estimate light input for timeseries.

# Time series is based on BATI data loaded above. 
# First full day is May 4th. Last full day is Sep 30.
start_date <- as.Date("2024-05-04", format = "%Y-%m-%d")
#end_date   <- "2024-09-30"
end_date   <- as.Date("2024-06-04", format = "%Y-%m-%d")

# Generate timestamps for 24 hours

ts_start <- as.POSIXct( paste0( start_date, " 00:00:00") )  # Start timeseries
ts_end   <- as.POSIXct( paste0( end_date,   " 23:00:00") )  # End timeseries
timestamps <- seq(from = ts_start, to = ts_end, by = "hour")

light_data <- CalculatePhotons( 
                  solar_elevation_angle( timestamps, latitude, longitude )
                  )
head(BATI5)

# Restrict BATI data to the start and end dates ... 
start_month_day <- month(start_date) * 100 + day(start_date)
end_month_day   <- month(end_date) * 100 + day(end_date)

subBATI <- BATI5[ (BATI5$month * 100 + BATI5$day) >= start_month_day & 
                 (BATI5$month * 100 + BATI5$day) <= end_month_day, ]

length(light_data)
dim(subBATI)


bdat <- subBATI
# Complete the input data frame by adding light
x <- dim(bdat)[[1]]
# Order of columns is temp, nutrients, DIC, light
grow_dat <- cbind( 
            "temp"  <- 282.15 + bdat$temp, # T is in Kelvin
            # Typical [NO3-] in coastal waters:  1 to 10 × 10⁻⁶ mol/
            "nutr"  <- rnorm(x, mean = 5e-6, sd = 5e-6),
            # Typical [DIC] in coastal waters:  1.8 to 2.5 × 10⁻³ mol/L
            "DIC"   <- rnorm(x, mean = 2.15e-3, sd = 2.15e-3),
            # Estimated photon concentration
#            "light" <- rnorm(dim(BATI5)[[1]], mean = 282.15, sd = 3.54)
# NOTE: light_data longer than BATI data by a few hundred. Maybe BATI missing a few days?
#       just take what you need for now ... shouldn't be off by much. 
            "light" <- light_data[1:x]
            )

# Ceiling any negative nutrient or DIC values 
grow_dat[, 2] <- ifelse( grow_dat[, 2] >=0, grow_dat[, 2], 0 )
grow_dat[, 3] <- ifelse( grow_dat[, 3] >=0, grow_dat[, 3], 0 )

PlotInputs( grow_dat, "May to September Environmental Conditions")
dim(grow_dat)

# Write results to a space-delimited text file
write.table(grow_dat, file = paste0(DEB_dir, "/env_conditions.txt"), sep = "  ", 
            col.names= FALSE, row.names = FALSE, quote = FALSE)

# Run the DEB model
setwd( DEB_dir )
system2( "rundeb_KELP.exe" )

# This output filename is currently fixed in the MatLab script
deb_output <- read.csv( paste0(DEB_dir, "/DEB_results.csv") )

#==== DEB Outputs ====
str(deb_output)
dim(deb_output) 
deb_output <- cbind( "days" = deb_output$Timet/24, deb_output ) 
PlotAllDEBResults( deb_output )

# Plot the state variables
#   (these continuously increasing include:
names( deb_output )

dev.off()
par( mfrow=c(2,3) )
plot( deb_output$M_Vt, type='l', main = "Structural mass", xlab="Hours", ylab = "Moles" )
plot( deb_output$L_allomt, type='l', main = "Total blade length", xlab="Hours", ylab = "cm" )
plot( deb_output$Wt, type='l', main = "Total weight", xlab="Hours", ylab = "grams" )
# Calculated internally as:   W = (w_V + m_EN * w_EN + m_EC * w_EC) * M_V

plot( deb_output$m_ENt, type='l', main = "Nitrogen reserve density", xlab="Hours", ylab = "M_EN / M_v" )
plot( deb_output$m_ECt, type='l', main = "Carbon reserve density", xlab="Hours", ylab = "M_EC / M_v" )
plot( deb_output$j_EC_G, type='l', main = "Flux of Carbon to growth", xlab="Hours" )


#==== Describe the fates of the carbon assimilated ====
# define index to final t 
l_idx <- dim(deb_output)[[1]]

# Compare initial and final mass of structure, already in g
deb_output$Wt[1]
deb_output$Wt[l_idx] 

# Proportional mass increase:
(deb_output$Wt[l_idx] - deb_output$Wt[1]) / deb_output$Wt[1]

# Now weight should be equal to Mv + m_EC + m_EN, all multiplied by their molar masses ... 
From Venolia - 
carbon = 30
nitro = 54
struct = 29.89


(deb_output$M_Vt[l_idx]  *27.51 ) + 
  (deb_output$M_Vt[l_idx] * deb_output$m_ENt[l_idx] * 17.0 ) + 
    (deb_output$M_Vt[l_idx] * deb_output$m_ECt[l_idx] *30.0)



# Understanding DIC flux:
# Point is that Ct vs Ct_j. Feel like one should be CI. 

hist( deb_output$M_Vt)
plot( deb_output$j_Ct )
plot( deb_output$j_Ct * deb_output$M_Vt )



# Total carbon assimilated [molC = sum( j_CAt [molC/molV/h] * M_Vt [molV]
mC_ass <- sum( deb_output$j_CAt * deb_output$M_Vt )
mN_ass <- sum( deb_output$j_NAt * deb_output$M_Vt )
gC_ass <- mC_ass *30 # grams / mole
gN_ass <- mN_ass *17

c_res <- deb_output$m_ECt[l_idx] * deb_output$M_Vt[l_idx] # C in reserve at end
c_str <-  sum( deb_output$j_EC_G )  # total C flux to  growth, i.e., stored in the kelp structure
(c_res + c_str) * 30


# Total uptake of DIC  (j_CI = molDIC/molV/h - this is )

deb_output$j_C

head(grow_dat)


Total C assimilated = sum of carbon assimilation per time step 
= sum of (assimilation rate * volume of kelp)

This total partitions to growth and maintenance. 

Total C fixed for growth = sum of carbon used for growth per time step
= sum( j_EC_G [molC/molV/h] * M_Vt [molV] )



# Cumulative weight of kelp frond ... 
c_wt <- cumsum( deb_output$dWt )
  

PlotGrowResults( deb_output, "May to September Growth Model Results")

  
ggplot(deb_output)+
#  geom_point(aes(x= Timet, y=j_EC_Gt),size=0.4, color="blue")
  geom_point(aes(x= Timet, y=L_allomt),size=0.4, color="blue")

#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 

today <- format(Sys.Date(), "%Y-%m-%d")

rmarkdown::render( "LSSM_documentation.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = DEB_dir,
                   output_file = paste0( "LSSM_DEB_testing_", today ))





#---- Romina's Plots to visualize the data ----

ggplot(ctd)+
  geom_point(aes(x= date_time, y=temperature_C*1.2),size=0.4, color="blue")+
  # lims(y=c(5, 22))+
  geom_point( aes(x= date_time, y=salinity_psu), size=0.4, color="darkgreen")+
  facet_wrap(~site, ncol = 2)+
  scale_y_continuous("Salinity (psu)", sec.axis = sec_axis(~ (.)/1.2, name = "Temperature (C)")) +
  # scale_x_continuous("Month", breaks = 1:7) +
  labs(x= "Date", y= "temperature (C)", title = "CTDs surface (0.5m depth) - Moorings")+
  theme_classic()


ggplot(ctd)+
  geom_point(aes(x= date_time, y=temperature_C, color=site),size=0.4, alpha=0.5)+
  geom_line(aes(x= date_time, y=temperature_C, color=site),size=0.5, alpha=0.5)+
  # lims(y=c(5, 22))+
  # geom_point( aes(x= date_time, y=salinity_psu), size=0.4, color="darkgreen")+
  facet_wrap(~cluster, ncol = 3)+
  # scale_y_continuous("Salinity (psu)", sec.axis = sec_axis(~ (.)/1.2, name = "Temperature (C)")) +
  # scale_x_continuous("Month", breaks = 1:7) +
  scale_color_manual(values= c("yellow1", "yellow3", "gold", "lightblue", "blue", "lightgrey", "darkgrey","gold3"))+
  labs(x= "Date", y= "temperature (C)", title = "CTDs surface (0.5m depth) - Moorings")+
  theme_classic()


ggplot(ctd)+
  geom_point(aes(x= date_time, y=salinity_psu, color=site),size=0.4, alpha=0.5)+
  geom_line(aes(x= date_time, y=salinity_psu, color=site),size=0.5, alpha=0.5)+
  facet_wrap(~cluster, ncol = 3)+
  # scale_y_continuous("Salinity (psu)", sec.axis = sec_axis(~ (.)/1.2, name = "Temperature (C)")) +
  # scale_x_continuous("Month", breaks = 1:7) +
  scale_color_manual(values= c("yellow2", "yellow3", "gold", "lightblue", "blue", "lightgrey", "darkgrey","gold3"))+
  labs(x= "Date", y= "Salinity (psu)", title = "CTDs surface (0.5m depth) salinity - Moorings")+
  theme_classic()




#FIN
