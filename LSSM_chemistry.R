#################################################################################
# Script:  LSSM_chemistry.R
# Created: May 27, 2025. EJG
# 
# Exploring the chemistry consequences of kelp growth.
# REQUIRES TWO TIME SERIES :
#   kelp_grow: a measure of daily kelp growth
#   daily_ocean: daily oceanic conditions relevant to carb() calculations. 
#     Combines mooring B5 with pCO2 from Ak ferry.
#
# NOTE: carb() requires vars in mol/kg, except pCO2 in uatm (micro atmospheres)
#   ALL conversions done above so oceanographic descriptors are carb() friendly.
# Ferry Data uses is "calibrated_SW_xCO2_dry"
# ## Updates: 
# 2025/11/13: Updated after restructing the data loading and parameter prep. 
#################################################################################

#----- Examine chemistry based on output of plant growth model -----

str(kelp_grow)
head(daily_ocean)
# carb() uses flags to ID inputs. Of interest here:

flag_CA <- 24 # pCO2 and Alkalinity 
flag_CD <- 25 # pCO2 and DIC  
flag_AD <- 15 # Alkalinity and DIC 

#----- Part 3a) Moles of DIC fixed by kelp ... 
# Starting with the daily kelp biomass, get daily grams DIC fixed.   
gDIC_fixed <- kelp_grow$B_kgWW / 1000 * wet_to_dry * dry_to_C  

#Using mol weights (g/mol) of structure and C reserve from DEB, estimate daily mols C
#Mol wt of structure = 27.51; mol wt of C reserve = 30. Use average
molDIC_fixed <- gDIC_fixed / (27.51+30) / 2

#Now the cumulative DIC fixed over 150 days (is the same shape as growth)
plot( kelp_grow$days, molDIC_fixed )

#Use diff() to return difference btwn consecutive elements. So now its DIC fixed/day
#This is equivalent to the rate of kelp growth. 
delkDIC <- c( 'NA', diff( molDIC_fixed )) 
#NOTE: Length is now day_stamps-1
plot( kelp_grow$days, delkDIC, ylab = "mol DIC fixed / day" )

#----- Part 3b) Ambient DIC and baseline pH

#----- Convert ambient pCO2 to DIC using pCO2 and totA. -----
oceanCarb <- carb( flag_CA, var1=daily_ocean$pCO2, var2=daily_ocean$totA,
                  S=daily_ocean$salt, T=daily_ocean$temp )

# Ambient DIC in 4 plots ... 
x <- "Ambient pCO2 to pH" 
par(mfrow=c(4,1))
par(mar=c(0,4,3,1) )
plot( daily_ocean$totA, type='l', xlab='', ylab='Alkalinity (mol/kg)', xaxt = "n", main = x)
par(mar=c(0,4,1,1) )
plot( daily_ocean$pCO2, type='l', xlab='', ylab='pCO2 (uatm)', xaxt = "n" )
par(mar=c(0,4,1,1) )
plot( x=daily_ocean$days, y=oceanCarb$DIC, type='l', xlab='', ylab='DIC (mol/kg)', xaxt = "n" )
par(mar=c(3,4,1,1) ) 
plot( x=daily_ocean$days, y=oceanCarb$pH, type='l', xlab='', ylab='pH')

#reset graphic layout 
par(mfcol = c(1,1))

# Compare pH using salinity (alkalinity), and temperature directly from the Columbia
x <- "Comparing pH from assembled ocean data (line) to MV Columbia sensor data (points)" 
csalt <- CalcAlk( ak_dat$Smn )
columbia_carb <- carb( flag_CA, var1=ak_dat$CO2mn , var2=csalt,
                   S=ak_dat$Smn, T=ak_dat$Tmn )

# Plot Columbia-based pH vs. daily_ocean data 
par( cex = 1.3) # suitable for cut/pste to ppt
plot( x=daily_ocean$days, y=oceanCarb$pH, type='l', xlab='', ylab='pH', main = x)
points( ak_dat$date, y=columbia_carb$pH, pch = 21, bg = "black")


#---- Now thinking about volume effects ----
# At what mass/volume is the reduction in DIC not significant?

# Need to start with the moles of DIC fixed per day by the kelp growth. 
# B_gDW needs to be adjusted for carbon portion of growth (0.25),
# and then for grams to moles (12.01) 
amb_ph  <- oceanCarb$pH
kmod_ph <- carb( flag_AD, var1=oceanCarb$ALK, var2=(oceanCarb$DIC-molDIC_fixed),
                 S=oceanCarb$S, T=oceanCarb$T)$pH

plot( daily_ocean$days, xlab="", ylab="pH",
      kmod_ph, type='l', col='green', main='Daily change in pH in 1 kg of water', lwd=2)
lines( daily_ocean$days, amb_ph, type='l', col='red', lwd=2 )
legend("topleft", legend = c("Ambient", "Kelp-affected"), col = c("red", "green"), lwd = 2)

# Above used to have more of an effect over time. Likely wrong given this:
plot( amb_ph - kmod_ph )
# POINT is that T and S values have a big effect on carb() i.e., kmod_ph above 
# calc'd with daily_ocean() values has massive effect.  


### NOTES:
# Buffering Capacity: The ratio of TA to DIC determines the buffering of pH changes. 
#   High TA regions will show smaller pH shifts.
# Seasonality: Your May–September timeline aligns with peak kelp growth but 
#   also with upwelling events that introduce new DIC-rich waters.


#---- Calculate effect of kelp growth on different water volumes -----

# Volumes in litres / kg of water to iterate over
volumes <- c(1, 10, 50, 100)

# This is the 'ambient' ocean conditions that we will use based on BATI5 and Ak ferry
ambient_data <- oceanCarb
oceanCarb$pH

# Create a list to store the outputs
adjusted_chem <- list()

# iterate over each volume
for(i in 1:length(volumes)) {
  
  # calculate DIC fixed per litre
  molDIC_fixed_per_kg <- molDIC_fixed / volumes[i]
  # Calculate carbonate chemistry - using Alkalinity and DIC, with DIC adjusted based on kelp carbon fixation
  chem_out <- carb(flag_AD, var1=ambient_data$ALK, var2=(ambient_data$DIC - molDIC_fixed_per_kg),
                   S=ambient_data$S, T=ambient_data$T )
  # label with the volume
  chem_out$Volume_kg <- volumes[i]
  # Calculate delta pH, the difference from ambient
  chem_out$delta_pH <-  chem_out$pH - ambient_data$pH
  # Calculate delta pCO2, the difference from ambient
  chem_out$delta_pCO2 <-  chem_out$pCO2 - ambient_data$pCO2
  # save output as list element
  adjusted_chem[[i]] <- chem_out
}

# combine the results in a single dataframe, add a day of year column too
# NB: This dataframe is essentially 160 days * # of water mass levels long - for ggplot()

# add date column to each df in list ... 
adjusted_chem <- lapply(
  adjusted_chem,
  function(df) {
    df$DoY <- kelp_grow$days
    df
  }
)

# merge the list of dfs ...
final_df <- dplyr::bind_rows(adjusted_chem)

# add a numeric DOY ... 
final_df$DoY_num <- as.numeric(format(final_df$DoY, "%j"))


ggplot(final_df) +
  geom_line(aes(x = DoY_num,
                y = delta_pH,
                colour = factor(Volume_kg),
                group = Volume_kg),
            size = 1.05) +
  theme_light() +
  labs(
    x = "Day of Year",
    colour = "Volume (kg)"
  )






# Surface plots -----------------------------------------------------------------

# pH over time, by volume
x <- ggplot(adjusted_chem) +
  geom_line(aes(x = DoY, y = delta_pH, colour = Volume_kg, group = Volume_kg)) +
  theme_light()
plot(x)

# Difference from Ambient
x <-ggplot(adjusted_chem) +
  geom_line(aes(x = DoY, y = delta_pH, colour = Volume_kg), size=1.1) + 
  labs(
    title = "Daily pH Change from Ambient for Different Water Masses",  # Add your title here
    x = "Day of Year",  # Change x-axis label
    y = "Daily Change (pH)", # Change y-axis label
    colour = "Water Mass (kg)" # Change legend title
  ) +
  theme_light() +
  theme(
    text = element_text(size = 16),  # Increases all text size
    axis.title = element_text(size = 18),  # Axis labels size
    axis.text = element_text(size = 14),   # Axis tick labels size
    axis.title.x = element_text(margin = margin(t = 10)),  # More space above x-axis title
    axis.title.y = element_text(margin = margin(r = 20)),   # More space to the right of y-axis title
    legend.text = element_text(size = 14), # Legend text size
    legend.title = element_text(size = 16) # Legend title size
  )  
plot(x)

x <- ggplot(adjusted_chem) +
  geom_point(aes(x = DIC, y = pH, colour = ALK)) +
  facet_wrap(~Volume_kg)
plot(x)

# library(plotly) # for 3D plots
# 
# glimpse(x1)
# 
# pairs(x1 %>% select(S, `T`, pH, pCO2, DIC, ALK), main = "Pairwise plot of daily chemistry parameters (x1; BATI5).")
# pairs(x2 %>% select(S, `T`, pH, pCO2, DIC, ALK), main = "Pairwise plot of daily chemistry parameters (x2; BATI6).")
# 
# 
# p <- plot_ly(data = adjusted_chem %>%
#                mutate(Volume_kg = as.numeric(as.character(Volume_kg))),
#              x = ~log10(Volume_kg), y = ~DIC, z = ~`pH`, color = ~`pH`, type = "scatter3d") %>%
#   layout(title = "Volume adjusted chemistry with kelp growth")
# 
# 
# htmlwidgets::saveWidget(p, file = "pH_DIC_by_volume.html")


# surface plot requires a matrix with evenly spaced values. example below
# m = matrix(rnorm(25),nrow = 5, ncol = 12)
# 
# plot_ly(z = ~m) %>% 
#   add_surface()


# Fin.
