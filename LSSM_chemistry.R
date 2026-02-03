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
# 2026/01/30: Revisiting process to more closely link to empirical work 
#################################################################################

#----- Examine chemistry based on output of plant growth model -----

#---- Required inputs:

# Kelp growth
str(kelp_df)

# Discrete samples to support inside/outside pH comparison
# Copy loading data from LSSM_Water_Quality.prj. discrete_sample_loading.R

# SOURCE modified during development to deal with inconsistent data:
#   RPtoRK : 2 codes for the same site
#   _B     : Dropped one CP replicate from May 13
x <- read_excel( paste0( source_dir, '/Discrete_sample_data/discrete_data_for_model_development_RPtoRK_B.xls' ))

# Pull and simplify names of a relevant subset of data
disc_data <- x[, c("Station_ID", "Collection_Date", "Collection_Time_PST",
                           "NIST_Temp", "YSI_S", "TA (umol/kg)", "pCO2@insituT (uatm)", "pH (total)")]
colnames(disc_data) <- c("Station", "Date", "Time", "Temp", "Salt", "TA", "pCO2", "pH")
head(disc_data)

#Fix Temp
disc_data$Temp <- as.numeric(disc_data$Temp)
# Adjust and format Date and Time columns
# Convert Excel fractional day value -> seconds -> 24 hr time
disc_data$Time <- hms::hms(round(as.numeric(disc_data$Time) * 86400))
# Format date
disc_data$Date <- as.Date(disc_data$Date)
# Add DateTime column
disc_data$DateTime <- as.POSIXct(disc_data$Date) + disc_data$Time

# Standardize station names 
#   Replace 'in' with '_In' and 'out' with '_Out' at the end of station names
disc_data$Station <- gsub("\\s*[iI][nN]$", "_In", disc_data$Station)
disc_data$Station <- gsub("\\s*[oO][uU][tT]$", "_Out", disc_data$Station)

# Check dataframe
head( disc_data )
class( disc_data )

#--- Working with kelp bed pairs, first pull all In and Outs.
#  Parse site and in/out flag from Station
inout_data <- disc_data %>%
  mutate(
    Site = str_remove(Station, "_(In|Out)$"),
    IO   = case_when(
      str_detect(Station, "_In$")  ~ "In",
      str_detect(Station, "_Out$") ~ "Out",
      TRUE ~ NA_character_
    )
  ) %>%
  # Keep only rows with valid suffix
  filter(!is.na(IO))

head( inout_data )
str( inout_data )  



# Next concrete step:
# 1) compute DIC and TA from (pH, pCO2, T, S) using seacarb::carb()
# 2) summarize In vs Out per site (mean, sd, n)
# 3) compute delta (Out - In) in DIC and pH for each site
# 4) quick plot of delta DIC by site

# df must contain: Site, IO (In/Out), Temp, Salt, pH, pCO2
# pCO2 assumed in µatm, Temp in °C, Salt in PSU

add_dic_from_ta_ph <- function(df) {
  res <- seacarb::carb(
    flag = 8,          # TA + pH given
    var1 = df$pH,      # pH (default seacarb scale, as you requested)
    var2 = df$TA,      # TA (µmol/kg)
    S    = df$Salt,
    T    = df$Temp,
    Patm = 1,
    P    = 1,    # Set for internal consistency
    Pt   = 0,    # Set for internal consistency
    Sit  = 0     # Set for internal consistency
  )
  df$DIC <- res$DIC    # µmol/kg
  df
}

summarize_in_out <- function(df) {
  
  req_vars <- c("Site", "IO", "Temp", "Salt", "TA", "pH", "pCO2")
  
  df_clean <- df %>%
    dplyr::select(all_of(req_vars)) %>%
    tidyr::drop_na()
  
  if (nrow(df_clean) < nrow(df)) {
    message("Dropped ", nrow(df) - nrow(df_clean),
            " rows with missing carbonate inputs.")
  }
  
  df2 <- add_dic_from_ta_ph(df_clean)
  
  # Summary by site x IO (this stays as-is)
  summ <- df2 %>%
    dplyr::group_by(Site, IO) %>%
    dplyr::summarise(
      n = dplyr::n(),
      Temp_mean = mean(Temp),
      Salt_mean = mean(Salt),
      
      TA_mean = mean(TA),
      TA_sd   = sd(TA),
      
      pH_mean = mean(pH),
      pH_sd   = sd(pH),
      
      pCO2_mean = mean(pCO2),
      pCO2_sd   = sd(pCO2),
      
      DIC_mean = mean(DIC),
      DIC_sd   = sd(DIC),
      .groups = "drop"
    )
  
  # ---- Deltas table (one row per Site) ----
  deltas <- summ %>%
    dplyr::select(Site, IO, n, Temp_mean, Salt_mean, TA_mean, pH_mean, pCO2_mean, DIC_mean) %>%
    tidyr::pivot_wider(
      id_cols = c(Site),                 # force one row per Site
      names_from = IO,
      values_from = c(n, Temp_mean, Salt_mean, TA_mean, pH_mean, pCO2_mean, DIC_mean),
      names_sep = "_"
    ) %>%
    dplyr::mutate(
      d_TA   = TA_mean_Out   - TA_mean_In,
      d_pH   = pH_mean_Out   - pH_mean_In,
      d_pCO2 = pCO2_mean_Out - pCO2_mean_In,
      d_DIC  = DIC_mean_Out  - DIC_mean_In
    )
  
  list(df_with_dic = df2, summary = summ, deltas = deltas)
}

plot_delta_DIC <- function(deltas_df) {
  ggplot(deltas_df, aes(x = Site, y = d_DIC)) +
    geom_col() +
    labs(
      x = "Site",
      y = expression(Delta*"DIC (Out - In, "*mu*"mol kg"^-1*")")
    ) +
    theme_classic()
}

plot_dic_boxplots <- function(df_with_dic) {
  
  ggplot(df_with_dic, aes(x = Site, y = DIC, fill = IO)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.8) +
    labs(
      x = "Site",
      y = expression("DIC ("*mu*"mol kg"^-1*")"),
      fill = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = "top"
    )
}

add_DIC_no_grouping <- function(df) {
  
  # Required variables for this workflow (TA-based closure)
  req_vars <- c("Site", "IO", "Temp", "Salt", "TA", "pH", "pCO2")
  
  df_clean <- df %>%
    dplyr::select(all_of(req_vars)) %>%
    tidyr::drop_na()
  
  if (nrow(df_clean) < nrow(df)) {
    message("Dropped ", nrow(df) - nrow(df_clean),
            " rows with missing carbonate inputs.")
  }
  
  # Compute DIC from TA + pH
  df2 <- add_dic_from_ta_ph(df_clean)
  df2
}

plot_dic_delta_boxplots <- function(df_with_dic) {
  
  # Expect columns: Site, Date, IO (In/Out), DIC
  req <- c("Site", "Date", "IO", "DIC")
  missing <- setdiff(req, names(df_with_dic))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Pair In/Out by Site + Date
  paired <- df_with_dic %>%
    select(Site, Date, IO, DIC) %>%
    pivot_wider(
      names_from = IO,
      values_from = DIC
    ) %>%
    drop_na(In, Out) %>%
    mutate(
      d_DIC = Out - In
    )
  
  if (nrow(paired) == 0) {
    stop("No valid In/Out pairs found after pairing.")
  }
  
  ggplot(paired, aes(x = Site, y = d_DIC, fill = Site)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21) +
    labs(
      x = "Site",
      y = expression(Delta*"DIC (Out - In, "*mu*"mol kg"^-1*")")
    ) +
    theme_classic() +
    theme(
      legend.position = "none"
    )
}




# ---- Run ----
results <- summarize_in_out(inout_data)

results$summary
results$deltas
names( results$deltas )
plot_delta_DIC(results$deltas)


#--- looking a bit closer ... 

x <- add_DIC_no_grouping(inout_data)
head(x)
plot_dic_boxplots( x )
plot_dic_delta_boxplots( x )















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
delkDIC <- c( 0, diff( molDIC_fixed )) 
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
x <- "pH from 'ambient' ocean data (line) and MV Columbia sensor data (points)" 
csalt <- CalcAlk( ak_dat$Smn )
columbia_carb <- carb( flag_CA, var1=ak_dat$CO2mn , var2=csalt,
                   S=ak_dat$Smn, T=ak_dat$Tmn )

# Plot Columbia-based pH vs. daily_ocean data 
par( cex = 1.5) # suitable for cut/pste to ppt
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
  theme_light(base_size = 18) +
  labs(
    x = "Day of Year",
    colour = "Water mass (kg)",
    title = "Daily pH chance from ambient for different sea water masses"
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
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
