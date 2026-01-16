#----------------------------------------------------------------------------
# Script:  LSSM_ODE.R
# Created: January 2026
## Purpose: Run a simple Dynamic Energy Budget-like model based on coupled, 
# ordinary differential equations. Written with support from Gemini.
# Updates:
#
#============================================================================

# Install the deSolve package if you don't have it
# install.packages("deSolve")
library(deSolve)
library(lubridate)
library(ggplot2)
library(tidyr)
library(dplyr)

#------------------------- Functions -----------------------------

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
    mutate(Datestamp = as.Date(Timestamp)) %>%
    # Group by each day
    group_by(Datestamp) %>%
    # Count hours where PAR is greater than 0
    summarize(
      day_hours   = sum(PAR > 0, na.rm = TRUE),
      night_hours = sum(PAR <= 0, na.rm = TRUE),
      .groups = "drop"
    )
}


# Define the kelp growth ODE system
# Adjusted ODE function based on Table 1 (Pontier et al. 2024)
grow_kelp_model <- function(t, state, pars, env_data) {
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

# First ODE model. One stipe and one blade
plot_kelp_growth <- function(ode_output_df) {
  
  # 1. Convert the ODE matrix output to a usable data frame
  df_plot <- as.data.frame(ode_output_df)
  
  # 2. Add total biomass column
  df_plot$B_total <- df_plot$B_th + df_plot$B_fr
  
  # 3. Reshape data from wide to long format for ggplot (Tidy Data)
  df_long <- tidyr::pivot_longer(
    df_plot,
    cols = starts_with("B_"),
    names_to = "Component",
    values_to = "Biomass_gDW"
  )
  
  # 4. Generate the plot
  ggplot(df_long, aes(x = time, y = Biomass_gDW, color = Component)) +
    geom_line(size = 1.2) +
    labs(
      title = "Simulated Bull Kelp Growth Dynamics (ODE Model)",
      x = "Time (Days)",
      y = "Biomass (g Dry Weight)",
      color = "Biomass Type"
    ) +
    scale_color_manual(values = c("B_fr" = "green", 
                                  "B_th" = "brown", 
                                  "B_total" = "black")) +
    theme_minimal()
}


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
    # Assumes your env_data has a column named 'day' or 'DOY'
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
    
    # --- Variable cloughing Logic (Sigmoid) ---
    # Create a smooth transition from low summer sloughing to high fall sloughing
    # Steepness is set to 0.1 (hardcoded) for a gradual 30-day shift
    current_slough <- slough_min + (slough_max - slough_min) / 
                      (1 + exp(-0.1 * (real_DOY - day_senescence)))
    
    # --- Stipe Growth (Still limited by structural capacity K_th) ---
    stipe_mult <- max(0, (1 - (B_th / K_th)^2))
    dB_th      <- r_max_th * fT * stipe_mult * B_th 
    
    # --- Stochastic Blade Logic  ---
    # GP = Gross Production (umol C / g DW / day) - GP_total is a vector of 42 values
    GP_day   <- (blade_pars$max_fix * fT * fL_blade) * day_h
    GP_night <- (dark_fixation * fT) * night_h 
    GP_total <- GP_day + GP_night
    
    # Exudation and Respiration losses based on measured Std Errors
    DOC_loss  <- (blade_pars$doc_day * fT * fL_blade * day_h) + (blade_pars$doc_night * fT * night_h)
    Resp_loss <- GP_total * blade_pars$resp_prop * fT_resp
    
    Net_C_Gain <- GP_total - DOC_loss - Resp_loss # umol C / gDW / day
    
    # Conversion to growth rate: (umol C -> mol C -> g C) / (g C / g DW)
    # 12.011 is the molar mass of Carbon
    tissue_growth <- (Net_C_Gain * 1e-6 * 12.011) / dry_to_C
    
    # Final Derivative: Growth - Sloughing
    # No r_max_fr cap applied here
    dB_fr_vector <- (tissue_growth * B_fr_vector) - (current_slough * B_fr_vector)
    
    return(list(c(dB_th, dB_fr_vector)))
  })
}


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

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

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

# --- Usage ---
# plot_kelp_dynamics(bio_df, env_daily, sen_day = 240)

#-----------Environmental data prep --------------------------------------------
temp_env <- data.frame( "Timestamp" = DST_focal1$DateTime, "Temp" = DST_focal1$Temp )
PAR_env  <- data.frame( "Timestamp" = par_df$Timestamp, "PAR" = par_df$mean )
env_hourly <- make_env_hourly( PAR_env, temp_env )
head(env_hourly)
  
DLI_env  <- data.frame( "Datestamp" = DLI_df$Date, "DLI" = DLI_df$mean )
temp_env <- data.frame( "Datestamp" = date( DST_focal1$DateTime ), "Temp" = DST_focal1$Temp )
env_daily <- make_env_daily( DLI_env, temp_env )
head(env_daily)

# Daylight hours for the days in env_daily 
day_hours <- get_daylight_hours( PAR_env )

# Join daylight hours to available daily environmental data 
env_daily <- env_daily %>%
  left_join(day_hours, by = "Datestamp")









#---- Parameterize and run the multiple blade ODE model -----
# 1. Individual plant setup
n_blades <- round(rnorm(1, mean = 41.5, sd = 13.2))

# 2. Draw unique metabolic parameters for EACH blade based on W&P
#   Reportes SEs are converted to SDs based on sample size (see report)
blade_population <- data.frame(
  max_fix   = rnorm(n_blades, 78.5, sd = 33.5),
  doc_day   = rnorm(n_blades, 10.8, sd = 4.51),
  doc_night = rnorm(n_blades, 3.65, sd = 3.84),
  resp_prop = rnorm(n_blades, 0.15, sd = 0.08)  # SE retained for stability
)

# 3. Biological "Guard Rails" (Truncation)
# Prevent negative physiology or impossible values
blade_population$max_fix   <- pmax(blade_population$max_fix, 10)  # Min fixation 10 umol
blade_population$doc_day   <- pmax(blade_population$doc_day, 0)   # No negative leakage
blade_population$doc_night <- pmax(blade_population$doc_night, 0) # No negative leakage

# Clamp respiration: 
# Min 0.05 (super efficient) to Max 0.50 (dying tissue)
blade_population$resp_prop <- pmin(pmax(blade_population$resp_prop, 0.05), 0.50)


model_params <- list(
  r_max_th       = 0.13,   # Max daily rate of stipe growth (in cm/day)
  dark_fixation  = 3.75,   # Light indepdt nighttime fixation - Weigel and Pfister
  DOC_leakage    = 0.162,  # Ave of daily fixed carbon lost as DOC - W & P
  resp_rate      = 0.10,   # Estimated 
  wet_to_dry     = 0.13,   # portion of kelp that is not water 
  dry_to_C       = 0.15,   # portion of dry kelp that is carbon
  K_th           = 2.0,    # The biomass (g DW) at with stipe stops growing
  slough_min     = 0.002,  # Summer baseline (0.2% loss/day)
  slough_max     = 0.10,   # Senescence peak (10% loss/day)
  day_senescence = 240     # The day (approx late August) when loss ramps up
)

# 3. Initial state
init_B_fr <- rep(0.1, n_blades) # Start with 0.01g DW per blade
initial_state <- c(B_th = 0.05, init_B_fr)

# 4. Define time sequence (from day 1 to day 100)
times <- env_daily$day

# d. Run the ODE
output <- ode(
  y = initial_state, 
  times = times, 
  func = grow_kelp_multiple_blades, 
  parms = model_params,
  env_data = env_daily,
  blade_pars = blade_population
)

#-------------- Visualization and diagnostics ---------------------------------

#Ensure output has multiple blades 
ncol( output )

first_blade_day1 <- output[2, 3]
last_blade_day1  <- output[2, ncol(output)]

final_row <- output[nrow(output), ]
manual_frond_sum <- sum(final_row[3:ncol(output)])
stipe_val <- final_row[2]

#----- Visualization using output as a df ----- 
out_df <- as.data.frame( output )

#--- 
bio_df <- get_biomass_timeseries( out_df )
kelp_plot <- plot_kelp_biomass(bio_df)
print(kelp_plot)


#--- Show blade distribution  ---
p_dist <- plot_blade_distribution( out_df )
print(p_dist)


# --- Plot kelp with sloughing ---
# Add real DOY to bio_df
start_day <- env_daily$day[1] # e.g., 135
bio_df$Real_DOY <- bio_df$time + start_day

plot_kelp_dynamics(bio_df, env_daily, sen_day = 240)



#---- Parameterize and run the single blade ODE model -----
# Define parameters (initial estimates, tune these with the paper data)
model_params <- list(
  r_max_th       = 0.13,    # Max daily rate of stipe growth (in cm/day)
  max_fixation   = 120.78,  # Peak umol C/g/h in June - Weigel and Pfister
  dark_fixation  = 3.75,    # Light indepdt nighttime fixation - Weigel and Pfister
  DOC_leakage    = 0.162,   # Ave of daily fixed carbon lost as DOC - W & P
  resp_rate      = 0.10,    # Estimated 
  carbon_fract   = 0.15,    # Estimated
  K_th           = 2.0,     # The biomass (g DW) at with stipe stops growing
  slough_rate_th = 0.001,
  slough_rate_fr = 0.01 
)

# 3. Define initial state (start with small kelp plant)
initial_state <- c(
  B_th = 0.005,  # initial stipe dry mass (g DW Carbon)
  B_fr = 0.5     # initial frond dry mass (g DW Carbon)
)

# 4. Define time sequence (from day 1 to day 100)
times <- env_daily$day

# 5. Solve the ODE system using the 'ode' function
output <- ode(
  times = times,
  y = initial_state,
  func = grow_kelp_model,
  parms = model_params,
  env_data = env_daily,
  method = "rk4" # Runge-Kutta 4th order method
)
#---------------------------------------------------------
str( output )
# View results
output[ 1:20,]

plot_kelp_growth( output )


