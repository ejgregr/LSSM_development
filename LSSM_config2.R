#----------------------------------------------------------------------------
# Script:  LSSM_config2.R
# Created: January 13 2026. 
# Purpose: Support second iteration of the nereo growth and chemistry work
#
# Notes:
#  - 
#================================== Load require packages =================================
# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "lubridate", "dplyr", "stringr",
                        "rmarkdown","knitr", "tinytex", "kableExtra", #currently causing trouble.
                        "seacarb", "gsw", "deSolve", "tidyr")

# Other packages that might be useful. 
# "tidyr", "raster", "stringr", "rasterVis",
# "RColorBrewer", "factoextra", "ggpubr", "cluster", 
# "diffeR", "vegan", "ranger", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
# "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr", "binr", 'gwxtab'

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)
version$version.string

#==== Configuration ====
bc_time <- 'America/Vancouver'

data_dir <- "C:/Data/Git/LSSM_development/Data"
DEB_dir  <- "C:/Data/Git/LSSM_water_analysis/Results"
PDF_dir  <- "C:/Data/Git/LSSM_development/PDF" # for Markdown script

# Load 2025 observational data. 
input_path <- file.path(DEB_dir, "kelp_project_data.RData" )
load(input_path)


#------------------------- Functions -----------------------------

# Updated ODE model with multiple blades and blade variability based on 
# measured standard errors (Weigel and Pfister)
grow_kelp_multiple_blades <- function(t, state, pars, env_data, blade_pars) {
  # state[1] is B_th (Stipe)
  # state[2:(n+1)] are the individual B_fr values (Blades)
  B_th        <- state[1]
  B_fr_vector <- state[-1]
  
  with(as.list(pars), {
    
    # 1. Calculate the ROW INDEX (For looking up data)
    # This maps simulation time (0, 1, 2) to row number (1, 2, 3)
    row_idx <- max(1, min(floor(t) + 1, nrow(env_data)))
    
    # 2. Extract the REAL Day of Year (DOY)
    # Requires env_data to have a column named 'day' or 'DOY'
    real_DOY <- env_data$day[row_idx]
    
    T_val   <- env_data$Temp[row_idx]
    I_val   <- env_data$DLI[row_idx]
    day_h   <- env_data$day_hours[row_idx]
    night_h <- env_data$night_hours[row_idx]
    
    # --- Environmental Scaling (Pontier et al. 2024) ---
    # Temperature scaling applies to both stipe and blade ... 
    fT <- ifelse(T_val <= 10.15, exp(-0.05*(T_val-10.15)^2), exp(-0.64*(T_val-10.15)))
    # DLI scaling applies to blade only ... 
    fL_blade <- ifelse(I_val < 20, I_val/20, ifelse(I_val > 40, max(0, 1-(I_val-40)/20), 1.0))
    # Respiration Scaling (Q10) ... 
    # R scales thermodynamically (doubles every 10C). T_ref = 10C
    fT_resp  <- 2.0^((T_val - 10) / 10)
    
    # --- Variable sloughing Logic (Sigmoid) ---
    # Create a smooth transition from low summer sloughing to high fall sloughing
    # Steepness is set to 0.1 (hardcoded) for a gradual 30-day shift
    current_slough <- slough_min + (slough_max - slough_min) / 
      (1 + exp(-0.1 * (real_DOY - day_senescence)))
    
    # --- Stipe Growth (Still limited by structural capacity K_th) ---
    stipe_mult <- max(0, (1 - (B_th / K_th)^2))
    dB_th      <- r_max_th * fT * stipe_mult * B_th 
    
    
    # --- Stochastic, staggered blade growth logic  ---
    # Create a vector of 0s and 1s. 
    # If time 't' < start_day, the blade is "dormant" (0).
    active_mask <- ifelse(t >= blade_pars$start_day, 1, 0)
    
    # GP = Gross Production (umol C / g DW / day) - GP_total is a vector of 1+nblades values
    # Growth applied to all blades, but dormant ones zeroed out at the end
    GP_day   <- (blade_pars$max_fix * fT * fL_blade) * day_h
    GP_night <- (dark_fixation * fT) * night_h 
    GP_total <- GP_day + GP_night
    
    # Exudation and Respiration losses based on measured Std Errors
    DOC_loss  <- (blade_pars$doc_day * fT * fL_blade * day_h) + 
      (blade_pars$doc_night * fT * night_h)
    Resp_loss <- GP_total * blade_pars$resp_prop * fT_resp
    Net_C_Gain <- GP_total - DOC_loss - Resp_loss # umol C / gDW / day
    
    # Conversion to growth rate: (umol C -> mol C -> g C) / (g C / g DW)
    # 12.011 is the molar mass of Carbon
    tissue_growth <- (Net_C_Gain * 1e-6 * 12.011) / dry_to_C
    
    # --- Final Derivative: Growth - Sloughing ---
    # APPLY THE MASK: Dormant blades get 0 growth and 0 sloughing
    dB_fr_vector <- (tissue_growth * B_fr_vector * active_mask) - 
      (current_slough * B_fr_vector * active_mask)
    
    return(list(c(dB_th, dB_fr_vector)))
  })
}

# Generate the geometric schedule for blade emergence. 
# Calculated for a set of 'blade generations'.
# Assigns a start_day to each blade. 
# Blades are added in each step using 2^(Gen-1).
assign_geometric_start_days <- function(n, interval) {
  start_days <- numeric(n)
  count_assigned <- 0
  generation <- 1
  
  while (count_assigned < n) {
    # How many blades appear in this new wave?
    # Wave 1 (Gen 1) adds 2 blades (The first split)
    # Wave 2 (Gen 2) adds 2 more (To reach 4)
    # Wave 3 (Gen 3) adds 4 more (To reach 8)
    # Wave 4 (Gen 4) adds 8 more (To reach 16)
    
    if (generation == 1) {
      wave_size <- 2
    } else {
      wave_size <- 2^(generation - 1)
    }
    
    # Don't exceed the total n
    n_actual <- min(wave_size, n - count_assigned)
    
    # Calculate the day for this wave
    # (generation - 1) * interval gives days: 0, 14, 28, 42...
    current_day <- 1 + ((generation - 1) * interval)
    
    # Fill the vector
    start_index <- count_assigned + 1
    end_index   <- count_assigned + n_actual
    start_days[start_index:end_index] <- current_day
    
    # Increment
    count_assigned <- count_assigned + n_actual
    generation <- generation + 1
  }
  return(start_days)
}


#---- Support Functions -----

make_env_hourly <- function(PAR_df, Temp_df, tz = "UTC",
                            fill_method = c("locf", "linear")) {
  fill_method <- match.arg(fill_method)
  
  parse_utc <- function(x) {
    if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(with_tz(x, tz))
    x_chr <- as.character(x)
    suppressWarnings(parse_date_time(
      x_chr,
      orders = c("Ymd HMS", "Y-m-d H:M:S", "Y-m-d H:M:S z", "Y-m-d H:M z"),
      tz = tz
    ))
  }
  
  # Parse + clean
  PAR_clean <- PAR_df %>%
    mutate(DateTime = parse_utc(Timestamp)) %>%
    filter(!is.na(DateTime)) %>%
    transmute(DateTime, PAR = as.numeric(PAR)) %>%
    arrange(DateTime)
  
  Temp_clean <- Temp_df %>%
    mutate(DateTime = parse_utc(Timestamp)) %>%
    filter(!is.na(DateTime)) %>%
    mutate(DateTime = floor_date(DateTime, "hour")) %>%
    group_by(DateTime) %>%
    summarise(temp = mean(Temp, na.rm = TRUE), .groups = "drop") %>%
    arrange(DateTime)
  
  if (nrow(PAR_clean) == 0) stop("No PAR rows after parsing Timestamp.")
  if (nrow(Temp_clean) == 0) stop("No Temp rows after parsing Timestamp.")
  
  # Overlap window (robust)
  par_rng  <- range(PAR_clean$DateTime,  na.rm = TRUE)
  tmp_rng  <- range(Temp_clean$DateTime, na.rm = TRUE)
  t0 <- max(par_rng[1], tmp_rng[1])
  t1 <- min(par_rng[2], tmp_rng[2])
  if (t0 > t1) {
    stop(sprintf("No overlap between PAR (%s to %s) and Temp (%s to %s).",
                 par_rng[1], par_rng[2], tmp_rng[1], tmp_rng[2]))
  }
  
  # Build an hourly grid over the overlap
  grid <- tibble(DateTime = seq(from = ceiling_date(t0, "hour"),
                                to   = floor_date(t1, "hour"),
                                by   = "hour"))
  
  # Join both series to the grid
  env_hourly <- grid %>%
    left_join(PAR_clean,  by = "DateTime") %>%
    left_join(Temp_clean, by = "DateTime") %>%
    arrange(DateTime)
  
  # Fill gaps
  if (fill_method == "locf") {
    env_hourly <- env_hourly %>%
      fill(PAR, temp, .direction = "downup")
  } else {
    # linear interpolation; requires at least 2 non-NA points
    env_hourly <- env_hourly %>%
      mutate(
        PAR  = as.numeric(PAR),
        temp = as.numeric(temp)
      )
    
    env_hourly$PAR  <- approx(x = env_hourly$DateTime[!is.na(env_hourly$PAR)],
                              y = env_hourly$PAR[!is.na(env_hourly$PAR)],
                              xout = env_hourly$DateTime,
                              rule = 2)$y
    env_hourly$temp <- approx(x = env_hourly$DateTime[!is.na(env_hourly$temp)],
                              y = env_hourly$temp[!is.na(env_hourly$temp)],
                              xout = env_hourly$DateTime,
                              rule = 2)$y
  }
  
  # Hour index for deSolve forcing
  env_hourly %>%
    mutate(time_h = as.numeric(difftime(DateTime, min(DateTime), units = "hours")))
}

# NOTE DLI is already daily. Temperature is averaged.
make_env_daily <- function(dli_df, temp_df) {
  # 1. Process temperature: Convert to Date and calculate daily mean
  temp_daily <- temp_df %>%
    mutate(Datestamp = as.Date(Datestamp)) %>%
    group_by(Datestamp) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE), .groups = "drop")
  
  # 2. Process DLI: Ensure Datestamp is Date type
  dli_daily <- dli_df %>%
    mutate(Datestamp = as.Date(Datestamp))
  
  # 3. Combine: Keep only dates where both measures exist
  env_combined <- inner_join(dli_daily, temp_daily, by = "Datestamp") %>%
    arrange(Datestamp) %>%
    # Add a relative 'day' column for the ODE solver (0, 1, 2...)
    mutate(day = as.numeric(Datestamp - min(Datestamp)))
  
  return(env_combined)
}

# Calc number of daylight hours per day based on hourly PAR
get_daylight_hours <- function( par ) {
  par %>%
    # Ensure Timestamp is POSIXct and extract the Date
    mutate(Datestamp = as.Date(Datestamp)) %>%
    # Group by each day
    group_by(Datestamp) %>%
    # Count hours where PAR is greater than 0
    summarize(
      day_hours   = sum(PAR > 0, na.rm = TRUE),
      night_hours = sum(PAR <= 0, na.rm = TRUE),
      .groups = "drop"
    )
}

# Post-ODE function to summarize biomass data for plotting
get_biomass_timeseries <- function(output_df) {
  # 1. Identify blade columns 
  # (everything except time and the stipe B_th)
  blade_cols <- setdiff(names(output_df), c("time", "B_th"))
  
  # 2. Create the clean wide-format dataframe
  biomass_data <- output_df %>%
    mutate(
      # Sum across the individual stochastic blades for each row
      Frond_biomass = rowSums(across(all_of(blade_cols))),
      # Use B_th for stipe
      Stipe_biomass = B_th,
      # Total sum
      Total_biomass = Frond_biomass + Stipe_biomass
    ) %>%
    # Select only the specific columns requested
    select(time, Frond_biomass, Stipe_biomass, Total_biomass)
  
  return(biomass_data)
}

# Define the kelp growth ODE system
# Adjusted ODE function based on Table 1 (Pontier et al. 2024)
simple_grow_kelp_model <- function(t, state, pars, env_data) {
  with(as.list(c(state, pars)), {
    
    # 1. Map time to data index
    current_day <- floor(t) + 1
    if (current_day > nrow(env_data)) current_day <- nrow(env_data)
    
    T_val <- env_data$Temp[current_day]
    I_val <- env_data$DLI[current_day]
    day_h   <- env_data$day_hours[current_day]
    night_h <- env_data$night_hours[current_day]    
    
    # --- Environmental Limitations from Pontier et al. (2024) ---
    #   TEMPERATURE: Peak at 10.15C, rapid 64% loss per 1C increase for both
    #     stipe and blade
    fT <- ifelse(T_val <= 10.15, 
                 exp(-0.05 * (T_val - 10.15)^2), 
                 exp(-0.64 * (T_val - 10.15)))
    
    # DLI effect on BLADE growth: Hump-shaped (Optimum 20-40)
    fL_blade <- ifelse(I_val < 20, 
                       I_val / 20, 
                       ifelse(I_val > 40, 
                              max(0, 1 - (I_val - 40) / 20),
                              1.0))
    
    # DLI effect on STIPE: Positive linear 
    fL_stipe <- function(I) {
      I_min <- 13.3 
      I_max <- 41.3
      if (I <= I_min) return(0)
      return(min(1, (I - I_min) / (I_max - I_min)))
    }
    
    # --- Carbon Dynamics (Weigel & Pfister 2021) ---
    # Gross Production (umol C / g DW / day)
    Day_Fixation     <- max_fixation  * fT  * fL_blade * day_h
    Night_Fixation   <- dark_fixation * fT *  night_h
    Gross_Production <- Day_Fixation + Night_Fixation
    
    # Losses (DOC and Respiration)
    # W & P suggests respiration is ~10% of gross fixation
    DOC_loss    <- Gross_Production * DOC_leakage
    Resp_loss   <- Gross_Production * respiration_rate 
    
    Net_C_Gain_umol <- Gross_Production - DOC_loss - Resp_loss
    
    # CONVERT: umol C/g DW -> g C gain / g DW
    # (1e-6 to get mol, * 12.01 to get grams)
    gC_per_gDW <- Net_C_Gain_umol * 1e-6 * 12.01
    
    # BIOMASS ACCUMULATION: Convert g Carbon to g Dry Weight
    specific_growth_rate <- gC_per_gDW / carbon_fraction
    
    # --- Net Growth (Biomass Change) ---
    dB_fr <- (specific_growth_rate * B_fr) - (slough_rate_fr * B_fr)
    
    # --- Stipe Logic from Potier with Multiplier Slowdown ---
    # We use a multiplier that stays near 1 until B_th gets close to K_th
    # This reflects the plant reaching the surface and shifting energy to blades
    stipe_multiplier <- (1 - (B_th / K_th) ^2) # Squared for a sharper drop-off near the end
    stipe_multiplier <- max(0, stipe_multiplier) # Ensure it never goes negative
    
    dB_th <- (r_max_th * fT * stipe_multiplier * B_th) - (slough_rate_th * B_th)
    
    # Return the derivatives AND the diagnostic variables
    # (The first elements of the list must be the derivatives)
    return(list(c(dB_th, dB_fr), 
                Gross_Prod = Gross_Production,
                Net_C_Gain = Net_C_Gain_umol,
                Growth_Rate_Fr = specific_growth_rate))
  })
}


#--------- Plotting Functions ----------
plot_blade_distribution <- function(output_df) {
  # 1. Extract the final time step
  final_step <- output_df[nrow(output_df), ]
  
  # 2. Identify blade columns (all columns except time, B_th, and any summary stats)
  # We assume blade columns are numeric IDs or started with "init_B_fr"
  # A safe way is to exclude known non-blade columns:
  blade_cols <- setdiff(names(final_step), c("time", "B_th", "Total_B_fr", "Total_Plant_DW", "mean_blade_size", "sd_blade_size"))
  
  # 3. Pivot the data to "long" format for plotting
  final_blades_long <- pivot_longer(final_step, 
                                    cols = all_of(blade_cols), 
                                    names_to = "blade_id", 
                                    values_to = "biomass")
  
  # 4. Create the plot
  p <- ggplot(final_blades_long, aes(x = biomass)) +
    geom_histogram(aes(y = ..density..), bins = 15, fill = "seagreen", color = "white", alpha = 0.7) +
    geom_density(color = "darkgreen", size = 1) +
    geom_vline(aes(xintercept = mean(biomass)), color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = mean(final_blades_long$biomass), y = 0, label = "Mean Size", 
             vjust = -1, color = "red", angle = 90) +
    labs(title = paste("Distribution of Final Blade Sizes (n =", length(blade_cols), "blades)"),
         subtitle = "Variability driven by Weigel & Pfister (2021) metabolic parameters",
         x = "Blade Biomass (g Dry Weight)",
         y = "Density") +
    theme_minimal()
  
  return(p)
}

plot_kelp_biomass <- function(bio_df) {
  # 1. Transform the data to "Long Format" for ggplot2
  # This moves the biomass types into a single column for easy color mapping
  plot_data <- bio_df %>%
    pivot_longer(cols = c(Frond_biomass, Stipe_biomass, Total_biomass),
                 names_to = "Component",
                 values_to = "Biomass")
  
  # 2. Create the plot
  p <- ggplot(plot_data, aes(x = time, y = Biomass, color = Component )) +
    geom_line(size = 1.1) +
    # Custom colors: Green for Fronds, Brown for Stipe, Black for Total
    scale_color_manual(values = c("Frond_biomass" = "green", 
                                  "Stipe_biomass" = "chocolate4", 
                                  "Total_biomass" = "black")) +
    
    labs(title = "Kelp Biomass Accumulation Over Time",
         subtitle = "Comparison of structural stipe and photosynthetic fronds",
         x = "Day of Year",
         y = "Dry Weight Biomass (g)" ) + #,
    #         color = "Plant Component",
    #         linetype = "Plant Component") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

plot_kelp_dynamics <- function(bio_df, env_daily, sen_day = 240) {
  
  # 1. Get the full Date sequence from env_daily
  # Assumes bio_df and env_daily have the same number of rows/steps
  bio_df$Date <- env_daily$Date[1:nrow(bio_df)]
  bio_df$DOY  <- yday(bio_df$Date) # Extract Day of Year for the math
  
  # 2. Define sloughing parameters locally
  slough_min <- 0.002
  slough_max <- 0.10
  
  # 3. Process data
  plot_data <- bio_df %>%
    mutate(
      # Calculate curve based on DOY from the Date column
      Slough_Rate = slough_min + (slough_max - slough_min) / 
        (1 + exp(-0.1 * (DOY - sen_day)))
    ) %>%
    select(Date, Total_biomass, Slough_Rate) %>%
    pivot_longer(cols = c(Total_biomass, Slough_Rate),
                 names_to = "Variable",
                 values_to = "Value") %>%
    mutate(Variable = factor(Variable, 
                             levels = c("Total_biomass", "Slough_Rate"),
                             labels = c("Total Biomass (g DW)", "Sloughing Rate (day^-1)")))
  
  # 4. Create Plot
  p <- ggplot(plot_data, aes(x = Date, y = Value)) +
    geom_line(aes(color = Variable), size = 1.2) +
    
    # Area fill
    geom_area(data = subset(plot_data, Variable == "Total Biomass (g DW)"),
              aes(fill = Variable), alpha = 0.3) +
    
    # Vertical Senescence Line (Need to convert DOY 240 to a Date for this year)
    # We grab the year from the first date in your data
    geom_vline(xintercept = as.Date(sen_day, origin = paste0(year(bio_df$Date[1])-1, "-12-31")), 
               linetype = "dashed", color = "grey50") +
    
    facet_wrap(~Variable, ncol = 1, scales = "free_y", strip.position = "left") +
    
    scale_color_manual(values = c("seagreen", "firebrick")) +
    scale_fill_manual(values = c("seagreen")) +
    
    # Use standard Date scaling on x-axis
    scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
    
    labs(title = "Seasonal Growth vs. Senescence",
         x = "Date",
         y = NULL) +
    
    theme_minimal() +
    theme(legend.position = "none",
          strip.placement = "outside",
          strip.text = element_text(face = "bold", size = 10))
  
  return(p)
}

plot_individual_tissues <- function(ode_output, env_daily) {
  
  # 1. Convert ODE output to Data Frame
  df <- as.data.frame(ode_output)
  
  # 2. Map 'time' to Dates
  #    Use match() to align simulation step with environment row
  df$Date <- env_daily$Date[match(df$time, env_daily$day)]
  
  # 3. Extract Stipe Data (Column 2)
  stipe_col_name <- colnames(df)[2]
  stipe_df <- df %>%
    select(Date, Biomass = all_of(stipe_col_name)) %>%
    mutate(Tissue = "Stipe")
  
  # 4. Extract Blade Data (Columns 3 to End)
  #    Pivot longer to get distinct lines for each blade
  blade_cols <- colnames(df)[3:(ncol(df)-1)] # Exclude 'Date'
  
  blade_df <- df %>%
    select(Date, all_of(blade_cols)) %>%
    pivot_longer(cols = -Date, names_to = "Blade_ID", values_to = "Biomass") %>%
    mutate(Tissue = "Blade")
  
  # 5. Plot with Log10 Scale
  p <- ggplot() +
    # A. Blades (Thin Green Lines)
    geom_line(data = blade_df, 
              aes(x = Date, y = Biomass, group = Blade_ID), 
              color = "forestgreen", size = 0.3, alpha = 0.5) +
    
    # B. Stipe (Thick Brown Line)
    geom_line(data = stipe_df, 
              aes(x = Date, y = Biomass), 
              color = "sienna", size = 1.5) +
    
    # C. Log Scale Transformation
    scale_y_log10(labels = scales::label_comma()) + 
    
    # D. Formatting
    scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
    labs(title = "Individual Tissue Biomass (Log Scale)",
         subtitle = "Visualizing the staggered emergence of blade cohorts",
         y = "Biomass (g DW) - Log10 Scale", x = "Date") +
    theme_minimal()
  
  return(p)
}
