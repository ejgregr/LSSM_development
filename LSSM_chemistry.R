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

# ## Updates: 
# 2025/xx/xx: 
#################################################################################

#----- Part 3 - Chemistry of the plant growth model -----
# carb() uses flags to ID inputs. Of interest here:

flag_CA <- 24 # pCO2 and Alkalinity 
flag_CD <- 25 # pCO2 and DIC  
flag_AD <- 15 # Alkalinity and DIC 

#----- Part 3a) Moles of DIC fixed by kelp ... 
# Starting with the daily kelp biomass, get daily grams DIC fixed.   
gDIC_fixed <- kelp_grow * 1000 * wet_to_dry * dry_to_C  

#Using mol weights (g/mol) of structure and C reserve from DEB, estimate daily mols C
#Mol wt of structure = 27.51; mol wt of C reserve = 30. Use average
molDIC_fixed <- gDIC_fixed / (27.51+30) / 2

#Now the cumulative DIC fixed over 150 days (is the same shape as growth)
#plot( day_stamps, molDIC_fixed )

#Use diff() to return difference btwn consecutive elements. So now its DIC fixed/day
#This is equivalent to the rate of kelp growth. 
delkDIC <- diff( molDIC_fixed ) 
#NOTE: Length is now day_stamps-1
#plot( day_stamps[-length(day_stamps)], delkDIC, ylab = "mol DIC fixed / day" )

#----- Part 3b) Ambient DIC and baseline pH

# For some reason, first day of PCO2 is NAN ... do need to fix the hour <--> day scaling. 
dim(daily_ocean)
daily_ocean$pCO2[1] <-mean(daily_ocean$pCO2[-1])


#----- Convert ambient pCO2 to DIC using pCO2 and totA. -----
oceanCarb <- carb( flag_CA, var1=daily_ocean$pCO2, var2=daily_ocean$totA,
                  S=daily_ocean$salt, T=daily_ocean$temp )

#Ambient DIC in 4 plots ... 
par(mfrow=c(4,1))
par(mar=c(0,4,3,1) )
plot( daily_ocean$totA, type='l', xlab='', ylab='Alkalinity (mol/kg)', xaxt = "n", main="Ambient pCO2 to pH" )
par(mar=c(0,4,1,1) )
plot( daily_ocean$pCO2, type='l', xlab='', ylab='pCO2 (uatm)', xaxt = "n" )
par(mar=c(0,4,1,1) )
plot( x=daily_ocean$days, y=oceanCarb$DIC, type='l', xlab='', ylab='DIC (mol/kg)', xaxt = "n" )
par(mar=c(3,4,1,1) ) 
plot( x=daily_ocean$days, y=oceanCarb$pH, type='l', xlab='', ylab='pH')


#----- some checks on DIC and pH inside and out  -----
# Compare carb() results for DIC from pCO2 for B5 (QCS) and B7 (KI). 

x1 <- carb( flag_CA, var1=daily_ocean$pCO2, var2=daily_TA$B5, S=day_BATI5$salt, T=day_BATI5$temp, 
             k1k2="m06", pHscale="SWS" )

# carb() fails (or did ...?) with salinity < 5 so make some days sufficiently salty ...  
#x1b <- carb( flag_CA, var1=daily_ocean$pCO2, var2=daily_TA$B5, S=pmax( day_BATI5$salt, 5 ), T=day_BATI5$temp, 
#             k1k2="m06", pHscale="SWS" )

x2 <- carb( flag_CA, var1=daily_ocean$pCO2, var2=daily_TA$B7, S=day_BATI7$salt, T=day_BATI7$temp )

# First, DIC at B5 vs B7 little different. So appears pCO2 is a main driver of DIC,
# regardless of Alkalinity. pCO2, Temp, and Salt constant
par(mfrow=c(3,2))
par(mar=c(0,4,3,1) )
plot(daily_TA$B5, type='l', xaxt = "n", main="pCO2 to pH - B5")
plot(daily_TA$B7, type='l', ylab="", xaxt = "n", main="pCO2 to pH - B7")

par(mar=c(0,4,1,1) )
plot(x1$DIC, type='l', xaxt = "n")
plot(x2$DIC, type='l', xaxt = "n", ylab="")

# Second, resulting pH seems really low, and why lower in KI?
par(mar=c(3,4,1,1) ) 
plot( day_stamps, x1$pH, type='l' )
plot( day_stamps, x2$pH, type='l', ylab="" )


#---- Thinking about volume effects ----
# Above is all based on 1 kg of water. 
# At what mass/volume is the reduction in DIC not significant?

# We can look at daily changes in pH from ambient attributable to kelp. 
# This effect is proportional to the size of the plant and its rate of growth. 
# (I don't think this is all captured with the logistic growth model)

amb_ph  <- carb( flag_AD, y, oceanCarb$DIC )$pH
kmod_ph <- carb( flag_AD, y, oceanCarb$DIC-delkDIC )$pH

par(mfrow=c(1,1), mar=c(4,4,2,2))
plot( day_stamps, xlab="", ylab="pH",
      amb_ph, type='l', col='red', main='Daily change in pH in 1 m² of water')
lines( day_stamps,
       kmod_ph, type='l', col='green' )
legend("bottomleft", legend = c("Ambient", "Kelp-affected"), col = c("red", "green"), lwd = 2)


### NOTES:
# Buffering Capacity: The ratio of TA to DIC determines the buffering of pH changes. 
#   High TA regions will show smaller pH shifts.
# Seasonality: Your May–September timeline aligns with peak kelp growth but 
#   also with upwelling events that introduce new DIC-rich waters.


#---- Calculate effect of kelp growth on different water volumes -----

# Volumes in litres / kg of water to iterate over
volumes <- c(10, 50, 100, 1000)

# This is the moles of DIC fixed per day by the kelp growth
molDIC_fixed

# This is the 'ambient' ocean conditions that we will use based on BATI5 and Ak ferry
# EJG May27: There was a reason you thought carb() needed to be recalculate something?
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

x <- adjusted_chem
# combine the results in a single dataframe, add a day of year column too
# NB: This dataframe is essentially 150 days * # of water mass levels long - for ggplot()
adjusted_chem <- lapply(adjusted_chem, function(x) x %>% mutate(DoY = day_stamps) )
adjusted_chem <- bind_rows(adjusted_chem)

# Make volume a factor
adjusted_chem$Volume_kg <- factor(adjusted_chem$Volume_kg)

#head(adjusted_chem)
#dim(adjusted_chem)

# Surface plots --------------------------------------------------------------------------------------------------------

# pH over time, by volume
x <- ggplot(adjusted_chem) +
  geom_line(aes(x = DoY, y = pH, colour = Volume_kg))+
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
