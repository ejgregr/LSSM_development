#----------------------------------------------------------------------------
# Script:  LSSM_config2.R
# Created: January 13 2026. 
# Purpose: Support second iteration of the nereo growth and chemistry work
#          Includes packages, paths, and functions.
# Notes:
#  - Feb 16: Added Chemistry functions. Lots of funky plotting. 
#
#================================== Load require packages =================================
# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "lubridate", "dplyr", "stringr",
                        "rmarkdown","knitr", "tinytex", "kableExtra", #currently causing trouble.
                        "seacarb", "gsw", "deSolve", "tidyr", "readxl", "akima" )

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
source_dir <- "C:/Data/Git/LSSM_water_analysis/Source_data"

# Load 2025 observational data. 
input_path <- file.path(DEB_dir, "kelp_project_data.RData" )
load(input_path)


# Bunch o code to load and prep the discrete water samples.
load_discrete_samples <- function(){
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
  
  return( inout_data )
}




#--------------------------- Light Support Functions ---------------------------

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


#---------------------------- DEB Plotting Functions ---------------------------

plot_kelp_biomass_WW <- function(out_df,
                              time_col="time", bth_col="B_th_WW", bfr_col="B_fr_WW",
                              log_y = FALSE) {
  stopifnot(all(c(time_col, bth_col, bfr_col) %in% names(out_df)))
  
  df_long <- data.frame(
    time    = rep(out_df[[time_col]], 3),
    biomass = c(out_df[[bfr_col]],
                out_df[[bth_col]],
                out_df[[bfr_col]] + out_df[[bth_col]]),
    series  = factor(rep(c("Frond", "Stipe", "Total"),
                         each = nrow(out_df)),
                     levels = c("Frond", "Stipe", "Total"))
  )
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = biomass,
                                             color = series)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = "Time (days)",
                  y = "Biomass (g WW)",
                  color = NULL,
                  linetype = NULL) +
    ggplot2::theme_classic()
  
  if (log_y) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  p
}

plot_kelp_biomass_DW <- function(out_df,
                              time_col="time", bth_col="B_th", bfr_col="B_fr",
                                 log_y = FALSE) {
  stopifnot(all(c(time_col, bth_col, bfr_col) %in% names(out_df)))
  
  df_long <- data.frame(
    time    = rep(out_df[[time_col]], 3),
    biomass = c(out_df[[bfr_col]],
                out_df[[bth_col]],
                out_df[[bfr_col]] + out_df[[bth_col]]),
    series  = factor(rep(c("Frond", "Stipe", "Total"),
                         each = nrow(out_df)),
                     levels = c("Frond", "Stipe", "Total"))
  )

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = biomass,
                                             color = series)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(x = "Time (days)",
                  y = "Biomass (g DW)",
                  color = NULL,
                  linetype = NULL) +
    ggplot2::theme_classic()
  
  if (log_y) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  p
}

# Shows all losses relative to gross C fixation
plot_C_losses <- function(df, time_col = "real_DOY", log_y = FALSE) {
  
  # Required columns
  req <- c( time_col, "GP_umolC", "DOC_umolC", "R_prod_umolC","R_maint_umolC","Slough_umolC" )
  missing <- setdiff(req, names(df))
  if (length(missing) > 0) {
    stop("Missing columns in data frame: ", paste(missing, collapse = ", "))
  }
  
  # Build long format
  long <- rbind(
#    data.frame(real_DOY = df[[time_col]], flux_type = "Gross production",        umolC = df$GP_umolC),
    data.frame(real_DOY = df[[time_col]], flux_type = "DOC leakage",             umolC = df$DOC_umolC),
    data.frame(real_DOY = df[[time_col]], flux_type = "Production respiration",  umolC = df$R_prod_umolC),
    data.frame(real_DOY = df[[time_col]], flux_type = "Maintenance respiration", umolC = df$R_maint_umolC),
    data.frame(real_DOY = df[[time_col]], flux_type = "Sloughing (as C)",        umolC = df$Slough_umolC)
#    data.frame(real_DOY = df[[time_col]], flux_type = "Light limitation",        umolC = df$fL),
#    data.frame(real_DOY = df[[time_col]], flux_type = "Temperature Limitation",  umolC = df$fT)
  )
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = real_DOY, y = umolC, color = flux_type)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Day of Year",
      y = expression(mu*"mol C d"^-1),
      color = NULL
    ) +
    ggplot2::theme_classic()
  
  if (log_y) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  return(p)
}

# Takes temp from env_daily df, calculates and plots associated scaling factor
plot_temperature_scaling <- function( env_dat ){
  
  # Plotting bits 
  DOY     <- env_dat$real_DOY
  T_daily <- env_dat$Temp
  fT_daily <- fT_reparam( T_daily )
  
  # Set up plot margins to allow second y-axis
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  # --- Left axis: Temperature ---
  plot( 
    DOY, T_daily,
    type = "l", col = "steelblue", lwd = 2,
    xlab = "Day of Year", ylab = "Temperature (°C)", bty = "l" )
  
  # --- Right axis: fT ---
  par(new = TRUE)
  
  plot(
    DOY, fT_daily,
    type = "l", col = "firebrick", lwd = 2,
    axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1) )
  
  axis(
    side=4, at=seq(0, 1, by=0.2), col = "firebrick", col.axis = "firebrick" )
  
  mtext(
    "Temperature scaling factor (fT)",
    side = 4,
    line = 3,
    col = "firebrick"
  )
} 

# Takes DLI from env_daily df, calculates and plots associated scaling factor
plot_DLI_scaling <- function( env_dat ){
  
  # Days to plot
  DOY <- env_dat$real_DOY
  DLI       <- env_dat$DLI
  fL_daily  <- fL_reparam( DLI )
  
  # Set up plot margins to allow second y-axis
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  # --- Left axis: DLI ---
  plot( 
    DOY, DLI,
    type = "l", col = "steelblue", lwd = 2,
    xlab = "Day of Year", ylab = "Daily Light Interval (mol m⁻² d⁻¹)", bty = "l" )
  
  # --- Right axis: fL ---
  par(new = TRUE)
  
  plot(
    DOY, fL_daily,
    type = "l", col = "firebrick", lwd = 2,
    axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1) )
  
  axis(
    side=4, at=seq(0, 1, by=0.2), col = "firebrick", col.axis = "firebrick" )
  
  mtext(
    "DLI scaling factor (fL)",
    side = 4,
    line = 3,
    col = "firebrick"
  )
} 


# Parameterless functions to plot the environmental curves.
plot_fT_curve <- function(T_opt = 10.0, sigma_warm = 1.4,
                          Tmin = 0, T_min = 7, T_max = 20, by = 0.1) {
  
  # Warm-side Gaussian penalty with cold-side plateau (Pontier-style)
  fT_gauss_warm <- function(T_C, T_opt = 10.0, sigma_warm = 1.4, Tmin = 0) {
    out <- ifelse(
      T_C <= T_opt,
      1.0,
      exp(- (T_C - T_opt)^2 / (2 * sigma_warm^2))
    )
    pmax(Tmin, pmin(1, out))
  }
  
  T_seq <- seq(T_min, T_max, by = by)
  df <- data.frame(
    T_C = T_seq,
    fT  = fT_gauss_warm(T_seq, T_opt = T_opt, sigma_warm = sigma_warm, Tmin = Tmin)
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = T_C, y = fT)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = expression("Temperature ("*degree*C*")"),
      y = expression(f[T])
    ) +
    ggplot2::theme_classic()
}


plot_fL_curve <- function(L_low = 20, L_high = 40, 
                          low_loss_per_10 = 0.23, high_loss_per_10 = 0.23,
                          DLI_min = 0, DLI_max = 70, DLI_step = 0.5
) {
  df <- data.frame(DLI = seq(DLI_min, DLI_max, by = DLI_step))
  df$fL <- fL_reparam(
    DLI = df$DLI,
    L_low = L_low,
    L_high = L_high,
    low_gain_per_10 = low_loss_per_10,
    high_loss_per_10 = high_loss_per_10
  )
  
  ggplot(df, aes(x = DLI, y = fL)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = c(L_low, L_high), linetype = "dashed") +
    labs(
      x = expression("Daily Light Integral (mol m"^-2*" d"^-1*")"),
      y = "Light scaling factor (fL)"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic()
}


#-------------------------- DEB Support Functions ---------------------------- 
# ChatGPT supported functions to scale T and DLI effects as per Pontier et al.


# Warm-side Gaussian penalty with a cold-side plateau. Simple representation of
# the warm temperature limitation identified by Pontier et al. 
fT_reparam <- function(T_C, T_opt = 10.0, sigma_warm = 1.4) {
  ifelse(
    T_C <= T_opt,
    1.0,
    exp(- (T_C - T_opt)^2 / (2 * sigma_warm^2))
  )
}

# JUST IN CASE want to add this extra penalty. There are some high temperatures. 
fT_gauss_warm_kink <- function(T_C, T_opt = 10.0, sigma_warm = 1.4,
                               T_kink = 12.0, extra_loss_per_C = 0.35) {
  base <- ifelse(
    T_C <= T_opt,
    1.0,
    exp(- (T_C - T_opt)^2 / (2 * sigma_warm^2))
  )
  # additional exponential penalty only above T_kink
  m_extra <- 1 - extra_loss_per_C
  k_extra <- -log(m_extra)
  kink <- ifelse(T_C <= T_kink, 1.0, exp(-k_extra * (T_C - T_kink)))
  pmax(0, pmin(1, base * kink))
}


# --- Helper for DLI scaling: convert "p% change per unit" into exponential rate constant ---
# If multiplier per unit is m (e.g., 0.77 per +1°C), then k = -log(m) / unit
k_from_multiplier <- function(multiplier, unit = 1) {
  -log(multiplier) / unit
}

# --- DLI scaling for blades: hump-shaped with low- and high-light penalties ---
# Approx ~23% gain per +10 DLI near ~20 (low-light side) and ~23% loss per +10 above ~40 (high-light side)
fL_reparam <- function(DLI, L_low = 20, L_high = 60, # was 40
                       low_gain_per_10 = 0.23,   # interpreted as "being 10 below optimal costs ~23%"
                       high_loss_per_10 = 0.10) {  # was 23
  
  # Convert to multipliers:
  # low side: at (L_low - 10) -> multiplier ~ (1 - 0.23) = 0.77
  m_low10  <- 1 - low_gain_per_10
  k_low    <- k_from_multiplier(m_low10, 10)   # per (mol m^-2 d^-1)
  
  # high side: at (L_high + 10) -> multiplier ~0.77
  m_high10 <- 1 - high_loss_per_10
  k_high   <- k_from_multiplier(m_high10, 10)  # per (mol m^-2 d^-1)
  
  out <- ifelse(
    DLI < L_low,
    exp(-k_low * (L_low - DLI)),           # decays below L_low
    ifelse(
      DLI <= L_high,
      1.0,                                 # optimum plateau
      exp(-k_high * (DLI - L_high))         # decays above L_high
    )
  )
  
  # Keep bounded [0,1]
  pmax(0, pmin(1, out))
}


#-------------------------- Chemistry Support Functions ------------------------ 
# Code provided by ChatGPT Feb 2026. Human digestion in progress. 
# The set of functions below is in the order that they are called to set up 
# and run the kelp incubation simulation. 


# ----------------- Carbonate chemistry Utilities -----------------------
# closure on discrete samples: (TA + pH) -> DIC
# Computes dissolved inorganic carbon (DIC, µmol kg⁻¹) from measured total 
# alkalinity and pH using seacarb::carb() given temperature and salinity. 
# Provides carbonate-system closure for discrete water samples.
add_dic_from_ta_ph <- function(df) {
  res <- seacarb::carb(
    flag = 8,
    var1 = df$pH,
    var2 = df$TA,      # umol/kg
    S    = df$Salt,
    T    = df$Temp,
    Patm = 1,
    P    = 1,          # ~1 m
    Pt   = 0,
    Sit  = 0
  )
  df$DIC <- res$DIC    # umol/kg
  df
}

# A wrapper to filter incomplete rows and append DIC to the full discrete 
# water dataset without aggregating by site or condition.
add_DIC_no_grouping <- function(df) {
  req_vars <- c("Station","Date","Time","Temp","Salt","TA","pCO2","pH","DateTime","Site","IO")
  df_clean <- df %>%
    dplyr::select(any_of(req_vars)) %>%
    tidyr::drop_na(Temp, Salt, TA, pH)
  
  if (nrow(df_clean) < nrow(df)) {
    message("Dropped ", nrow(df) - nrow(df_clean),
            " rows with missing carbonate inputs.")
  }
  
  add_dic_from_ta_ph(df_clean)
}


# Computes pH from total alkalinity and DIC (holding TA constant) using seacarb::carb().
# Translates modelled DIC changes into pH responses.
ph_from_ta_dic <- function(TA_umolkg, DIC_umolkg, S, T_C, P_dbar = 1) {
  # You confirmed this flag is correct in your setup.
  res <- seacarb::carb(
    flag = 15,            # TA + DIC -> pH (per your confirmation)
    var1 = TA_umolkg,
    var2 = DIC_umolkg,
    S    = S,
    T    = T_C,
    Patm = 1,
    P    = P_dbar,
    Pt   = 0,
    Sit  = 0
  )
  res$pH
}


# ----------------- Water-sample pairing and preparation  -----------------------

# Pairs Inside / Outside water samples by Site & Date, giving a single row 
# per paired event w clearly labeled _In / _Out carbonate and physical variables.
make_pairs_out_in <- function(water_with_dic) {
  water_with_dic %>%
    select(Station, Site, Date, DateTime, IO, Temp, Salt, TA, pH, pCO2, DIC) %>%
    pivot_wider(
      id_cols = c(Site, Date),
      names_from = IO,
      values_from = c(Station, DateTime, Temp, Salt, TA, pH, pCO2, DIC),
      names_sep = "_"
    ) %>%
    drop_na(DIC_In, DIC_Out, TA_Out, pH_In, Temp_Out, Salt_Out, DateTime_Out)
}

# =============================================================================
# 2) Compute biomass-specific day/night u_net from DEB + env_daily
#
# Assumptions (per your update):
# - DEB includes net_daily_gDW (net inorganic C removal rate per gDW per day)
#   Units expected: umol C gDW^-1 day^-1
# - Day/night partitioning of net follows the ratio of GP_day_umolC : GP_night_umolC
# - Use day_hours and night_hours from env_daily (no estimation)
# =============================================================================

# Partitions the DEB-provided net daily carbon uptake per gram dry weight (net_daily_gDW) 
# into day and night hourly uptake rates, using the ratio of gross primary production 
#during day and night and externally supplied day/night hours.
compute_unet_day_night_from_deb <- function(row,
                                            net_col   = "net_daily_gDW",
                                            gp_day_col   = "GP_day_umolC",
                                            gp_night_col = "GP_night_umolC",
                                            day_hours_col   = "day_hours",
                                            night_hours_col = "night_hours") {
  
  net_daily_gdw <- row[[net_col]]            # umol C gDW^-1 day^-1
  gp_day <- row[[gp_day_col]]
  gp_night <- row[[gp_night_col]]
  
  day_h <- row[[day_hours_col]]
  night_h <- row[[night_hours_col]]
  
  gp_sum <- gp_day + gp_night
  # If GP totals are missing/zero, fall back to proportional-to-hours partition
  if (is.na(gp_sum) || gp_sum <= 0) {
    frac_day <- day_h / (day_h + night_h)
  } else {
    frac_day <- gp_day / gp_sum
  }
  frac_night <- 1 - frac_day
  
  # Phase-specific net (umol C gDW^-1 per phase)
  net_day_gdw   <- net_daily_gdw * frac_day
  net_night_gdw <- net_daily_gdw * frac_night
  
  # Convert to hourly rates (umol C gDW^-1 h^-1)
  u_day   <- net_day_gdw   / day_h
  u_night <- net_night_gdw / night_h
  
  list(u_day = u_day, u_night = u_night, frac_day = frac_day, day_h = day_h, night_h = night_h)
}


#---------------------- Tea box forward model  ------------------------------
# Assumptions: - No gas exchange, Well mixed, TA constant
#              - DIC updated by u_net(t)*biomass/water_mass
#---------------------- Tea box forward model  ------------------------------

# Classify a timestamp as day or night based on clock hour. Used to switch 
# between day and night carbon uptake rates during the incubation simulation.
is_day_hour <- function(datetime) {
  h <- lubridate::hour(datetime)
  h >= 6 & h < 18
}


# Simulates the time evolution of DIC and pH in a well-mixed, fixed-mass water 
# parcel exposed to kelp biomass for a specified duration, w above assumptions. 
run_tea <- function(init_TA_umolkg,
                    init_DIC_umolkg,
                    S, T_C,
                    water_mass_kg,
                    kelp_biomass_gdw,
                    exposure_hours,
                    u_day_umolC_gdw_h,
                    u_night_umolC_gdw_h,
                    start_datetime) {
  
  stopifnot(exposure_hours <= 24)
  
  times_h <- c(0, seq(1, floor(exposure_hours), by = 1))
  if (tail(times_h, 1) < exposure_hours) times_h <- c(times_h, exposure_hours)
  
  DIC <- numeric(length(times_h))
  pH  <- numeric(length(times_h))
  
  DIC[1] <- init_DIC_umolkg
  pH[1]  <- ph_from_ta_dic(init_TA_umolkg, DIC[1], S, T_C)
  
  for (i in 2:length(times_h)) {
    dt <- times_h[i] - times_h[i - 1]
    now_dt <- start_datetime + hours(times_h[i - 1])
    
    u_net <- if (is_day_hour(now_dt)) u_day_umolC_gdw_h else u_night_umolC_gdw_h
    
    dDIC_umolkg <- -(u_net * kelp_biomass_gdw * dt) / water_mass_kg
    DIC[i] <- DIC[i - 1] + dDIC_umolkg
    pH[i]  <- ph_from_ta_dic(init_TA_umolkg, DIC[i], S, T_C)
  }
  
  tibble(
    t_h = times_h,
    DateTime = start_datetime + hours(times_h),
    TA_umolkg = init_TA_umolkg,
    DIC_umolkg = DIC,
    pH = pH
  )
}

# Joins paired water samples to DEB day & env_daily day/night hrs. Select by DOY
# pulled from water sample. DEB provides kelp carbon uptake rates and day/night 
# photoperiod info to each sampling event.
prepare_pairs_with_deb_env <- function(pairs_df,
                                       kelp_df,
                                       env_daily,
                                       deb_doy_col = "real_DOY",
                                       env_doy_col = "real_DOY") {
  
  pairs_df %>%
    mutate(DOY = yday(Date)) %>%
    left_join(kelp_df %>% mutate(DOY = .data[[deb_doy_col]]), by = "DOY") %>%
    left_join(env_daily %>% mutate(DOY = .data[[env_doy_col]]) %>%
                select(DOY, day_hours, night_hours),
              by = "DOY")
}


# --------------------- landscape simulations ------------------------

# Runs the kelp-tea box model across a Cartesian grid of water mass, kelp biomass, 
# and exposure time for each valid Site × Date pairing, returning modeled pH changes, 
# DIC changes, and diagnostics while capturing and reporting failures for debugging.
simulate_landscape <- function(pairs_df,
                               kelp_df,
                               env_daily,
                               water_mass_kg_vec,
                               kelp_biomass_gdw_vec,
                               exposure_hours_vec,
                               verbose = TRUE) {
  
  msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # ---- 1) Join paired water to DEB + env ----
  msg("[1/6] Joining paired water samples to DEB + env_daily ...")
  pairs2 <- prepare_pairs_with_deb_env(pairs_df, kelp_df, env_daily)
  
  # ---- 2) Guard rails + basic validity filtering (drop rows that will break) ----
  msg("[2/6] Applying guard rails & dropping invalid rows ...")
  n0 <- nrow(pairs2)
  
  pairs2 <- pairs2 %>%
    mutate(
      gp_sum    = GP_day_umolC + GP_night_umolC,
      hours_sum = day_hours + night_hours
    ) %>%
    filter(
      # Required to run tea + carbonate calc
      !is.na(TA_Out), !is.na(DIC_Out),
      !is.na(Salt_Out), !is.na(Temp_Out),
      !is.na(DateTime_Out),
      !is.na(pH_Out), !is.na(pH_In),
      !is.na(DIC_In),
      
      # Required for DEB partition
      !is.na(net_daily_gDW),
      !is.na(day_hours), !is.na(night_hours),
      hours_sum > 0, day_hours > 0, night_hours > 0,
      
      # We'll allow gp_sum <= 0 (falls back to hours partition),
      # but we do NOT allow gp_sum = NA because it can propagate NAs in seacarb checks.
      !is.na(gp_sum)
    ) %>%
    select(-gp_sum, -hours_sum)
  
  msg("  Kept ", nrow(pairs2), " / ", n0, " paired events after filtering (dropped ", n0 - nrow(pairs2), ").")
  
  if (nrow(pairs2) == 0) stop("No valid paired events after join + guard rails.")
  
  # ---- 3) Expand grid of experimental treatments ----
  msg("[3/6] Expanding parameter grid (water_mass × kelp_biomass × exposure_hours) ...")
  grid <- tidyr::expand_grid(
    Site = unique(pairs2$Site),
    Date = unique(pairs2$Date),
    water_mass_kg = water_mass_kg_vec,
    kelp_biomass_gdw = kelp_biomass_gdw_vec,
    exposure_hours = exposure_hours_vec
  ) %>%
    left_join(pairs2, by = c("Site", "Date"))
  
  msg("  Grid rows: ", nrow(grid))
  
  # ---- 4) Compute u_day/u_night *outside* the big mutate, as an intermediate table ----
  msg("[4/6] Computing u_day/u_night from DEB net_daily_gDW + GP ratios ...")
  
  # Identify rows that still have missing required fields after join to grid
  req <- c("net_daily_gDW","GP_day_umolC","GP_night_umolC","day_hours","night_hours",
           "TA_Out","DIC_Out","Salt_Out","Temp_Out","DateTime_Out","pH_Out","pH_In","DIC_In")
  missing_req <- grid %>%
    mutate(.row_id = row_number()) %>%
    filter(if_any(all_of(req), is.na)) %>%
    select(.row_id, Site, Date, any_of(req))
  
  if (nrow(missing_req) > 0) {
    msg("  NOTE: ", nrow(missing_req), " grid rows have missing required inputs (will be dropped).")
  }
  
  grid_ok <- grid %>%
    mutate(.row_id = row_number()) %>%
    filter(!if_any(all_of(req), is.na))
  
  msg("  Rows with complete inputs: ", nrow(grid_ok), " / ", nrow(grid))
  
  # Compute rates row-by-row, but store them in a clean intermediate
  rates_tbl <- grid_ok %>%
    rowwise() %>%
    mutate(
      .rates = list(compute_unet_day_night_from_deb(cur_data())),
      u_day   = .rates$u_day,
      u_night = .rates$u_night,
      frac_day = .rates$frac_day,
      day_h   = .rates$day_h,
      night_h = .rates$night_h
    ) %>%
    ungroup() %>%
    select(.row_id, u_day, u_night, frac_day, day_h, night_h)
  
  # Join rates back
  grid_ok <- grid_ok %>%
    left_join(rates_tbl, by = ".row_id")
  
  # ---- 5) Run tea in a controlled loop with progress messages & per-row error trapping ----
  msg("[5/6] Running tea simulations row-by-row (with error capture) ...")
  
  n_sim <- nrow(grid_ok)
  if (n_sim == 0) stop("No simulations to run after filtering.")
  
  # Preallocate result list (faster + easier to debug)
  out_list <- vector("list", n_sim)
  
  for (i in seq_len(n_sim)) {
    
    if (verbose && (i %% 250 == 0 || i == 1 || i == n_sim)) {
      msg("  Tea run ", i, " / ", n_sim)
    }
    
    r <- grid_ok[i, ]
    
    # Extra debug guard: seacarb sometimes errors if S/T are NA or out-of-range
    # Your error indicates a missing value in a logical check (often NA T or S).
    if (is.na(r$Temp_Out) || is.na(r$Salt_Out)) {
      out_list[[i]] <- tibble(.row_id = r$.row_id, ok = FALSE,
                              err = "Temp_Out or Salt_Out is NA (should have been filtered).")
      next
    }
    
    # Run tea with tryCatch so you can inspect failures instead of crashing
    tea_res <- tryCatch(
      run_tea(
        init_TA_umolkg   = r$TA_Out,
        init_DIC_umolkg  = r$DIC_Out,
        S = r$Salt_Out, T_C = r$Temp_Out,
        water_mass_kg = r$water_mass_kg,
        kelp_biomass_gdw = r$kelp_biomass_gdw,
        exposure_hours = r$exposure_hours,
        u_day_umolC_gdw_h = r$u_day,
        u_night_umolC_gdw_h = r$u_night,
        start_datetime = r$DateTime_Out
      ),
      error = function(e) e
    )
    
    if (inherits(tea_res, "error")) {
      out_list[[i]] <- tibble(
        .row_id = r$.row_id,
        ok = FALSE,
        err = conditionMessage(tea_res)
      )
      next
    }
    
    # Successful run: compute endpoints
    pH_end  <- tea_res$pH[nrow(tea_res)]
    DIC_end <- tea_res$DIC_umolkg[nrow(tea_res)]
    
    out_list[[i]] <- tibble(
      .row_id = r$.row_id,
      ok = TRUE,
      pH_end = pH_end,
      DIC_end = DIC_end,
      dpH_end = pH_end - r$pH_Out,
      dpH_needed = r$pH_In - r$pH_Out,
      dDIC_removed = r$DIC_Out - DIC_end,
      dDIC_needed  = r$DIC_Out - r$DIC_In
    )
  }
  
  res_tbl <- bind_rows(out_list)
  
  # ---- 6) Assemble final output tables ----
  msg("[6/6] Assembling outputs ...")
  
  # Attach results back to grid_ok (successful + failed)
  sim_out <- grid_ok %>%
    left_join(res_tbl, by = ".row_id")
  
  # A clean table of failures for debugging
  failures <- sim_out %>%
    filter(!ok) %>%
    select(.row_id, Site, Date, water_mass_kg, kelp_biomass_gdw, exposure_hours,
           Temp_Out, Salt_Out, TA_Out, DIC_Out, pH_Out,
           net_daily_gDW, GP_day_umolC, GP_night_umolC, day_hours, night_hours,
           u_day, u_night,
           any_of("err") )  #any_of() added to fix a select error as 'err' may not propagate
  
  if (nrow(failures) > 0) {
    msg("  WARN: ", nrow(failures), " simulations failed. Inspect the returned `$failures` table.")
  } else {
    msg("  All simulations succeeded.")
  }
  
  # Keep only successful sims in the main output (typical use),
  # but return both so you can debug.
  sims <- sim_out %>%
    filter(ok) %>%
    select(Site, Date, DOY,
           water_mass_kg, kelp_biomass_gdw, exposure_hours,
           TA_Out, DIC_Out, pH_Out, DIC_In, pH_In,
           dpH_needed, pH_end, dpH_end,
           dDIC_needed, dDIC_removed,
           net_daily_gDW, frac_day, day_h, night_h, u_day, u_night)
  
  list(
    sims = sims,
    failures = failures,
    grid_ok = grid_ok,
    missing_required = missing_req
  )
}

# As above, but working with kelp density
# Feb 16: Currently using this version. 
simulate_dens_landscape <- function(pairs_df,
                                    kelp_df,
                                    env_daily,
                                    kelp_density_gdwkg_vec,     # <-- NEW: gDW / kg
                                    exposure_hours_vec,
                                    water_mass_kg = 1,          # keep fixed (or pass a vector if you want)
                                    verbose = TRUE) {
  
  msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # ---- 1) Join paired water samples to DEB + env ----
  msg("[1/6] Joining paired water samples to DEB + env_daily ...")
  pairs2 <- prepare_pairs_with_deb_env(pairs_df, kelp_df, env_daily)
  
  # ---- 2) Guard rails + validity filtering ----
  msg("[2/6] Applying guard rails & dropping invalid rows ...")
  n0 <- nrow(pairs2)
  
  pairs2 <- pairs2 %>%
    mutate(
      gp_sum    = GP_day_umolC + GP_night_umolC,
      hours_sum = day_hours + night_hours
    ) %>%
    filter(
      # Required to run tea + carbonate calc
      !is.na(TA_Out), !is.na(DIC_Out),
      !is.na(Salt_Out), !is.na(Temp_Out),
      !is.na(DateTime_Out),
      !is.na(pH_Out), !is.na(pH_In),
      !is.na(DIC_In),
      
      # Required for DEB partition
      !is.na(net_daily_gDW),
      !is.na(day_hours), !is.na(night_hours),
      hours_sum > 0, day_hours > 0, night_hours > 0,
      
      # Allow gp_sum <= 0 (fallback to hours), but don't allow NA
      !is.na(gp_sum)
    ) %>%
    select(-gp_sum, -hours_sum)
  
  msg("  Kept ", nrow(pairs2), " / ", n0, " paired events after filtering (dropped ", n0 - nrow(pairs2), ").")
  if (nrow(pairs2) == 0) stop("No valid paired events after join + guard rails.")
  
  # ---- 3) Expand grid in density space ----
  msg("[3/6] Expanding parameter grid (density × exposure × water_mass) ...")
  
  # allow water_mass_kg to be scalar or vector
  water_mass_kg_vec <- water_mass_kg
  
  grid <- tidyr::expand_grid(
    Site = unique(pairs2$Site),
    Date = unique(pairs2$Date),
    water_mass_kg = water_mass_kg_vec,
    kelp_density_gdwkg = kelp_density_gdwkg_vec,
    exposure_hours = exposure_hours_vec
  ) %>%
    left_join(pairs2, by = c("Site", "Date")) %>%
    mutate(
      kelp_biomass_gdw = kelp_density_gdwkg * water_mass_kg
    )
  
  msg("  Grid rows: ", nrow(grid))
  
  # ---- 4) Compute u_day/u_night as an intermediate table ----
  msg("[4/6] Computing u_day/u_night from DEB net_daily_gDW + GP ratios + env photoperiod ...")
  
  req <- c("net_daily_gDW","GP_day_umolC","GP_night_umolC","day_hours","night_hours",
           "TA_Out","DIC_Out","Salt_Out","Temp_Out","DateTime_Out","pH_Out","pH_In","DIC_In",
           "kelp_density_gdwkg","water_mass_kg","kelp_biomass_gdw","exposure_hours")
  
  missing_req <- grid %>%
    mutate(.row_id = row_number()) %>%
    filter(if_any(all_of(req), is.na)) %>%
    select(.row_id, Site, Date, any_of(req))
  
  if (nrow(missing_req) > 0) {
    msg("  NOTE: ", nrow(missing_req), " grid rows have missing required inputs (will be dropped).")
  }
  
  grid_ok <- grid %>%
    mutate(.row_id = row_number()) %>%
    filter(!if_any(all_of(req), is.na))
  
  msg("  Rows with complete inputs: ", nrow(grid_ok), " / ", nrow(grid))
  if (nrow(grid_ok) == 0) stop("No simulations to run after filtering missing inputs.")
  
  rates_tbl <- grid_ok %>%
    rowwise() %>%
    mutate(
      .rates = list(compute_unet_day_night_from_deb(cur_data())),
      u_day   = .rates$u_day,
      u_night = .rates$u_night,
      frac_day = .rates$frac_day,
      day_h   = .rates$day_h,
      night_h = .rates$night_h
    ) %>%
    ungroup() %>%
    select(.row_id, u_day, u_night, frac_day, day_h, night_h)
  
  grid_ok <- grid_ok %>%
    left_join(rates_tbl, by = ".row_id")
  
  # ---- 5) Run tea simulations row-by-row with error capture ----
  msg("[5/6] Running tea simulations row-by-row (with error capture) ...")
  
  n_sim <- nrow(grid_ok)
  out_list <- vector("list", n_sim)
  
  for (i in seq_len(n_sim)) {
    
    if (verbose && (i %% 250 == 0 || i == 1 || i == n_sim)) {
      msg("  Tea run ", i, " / ", n_sim)
    }
    
    r <- grid_ok[i, ]
    
    tea_res <- tryCatch(
      run_tea(
        init_TA_umolkg   = r$TA_Out,
        init_DIC_umolkg  = r$DIC_Out,
        S = r$Salt_Out, T_C = r$Temp_Out,
        water_mass_kg = r$water_mass_kg,
        kelp_biomass_gdw = r$kelp_biomass_gdw,
        exposure_hours = r$exposure_hours,
        u_day_umolC_gdw_h = r$u_day,
        u_night_umolC_gdw_h = r$u_night,
        start_datetime = r$DateTime_Out
      ),
      error = function(e) e
    )
    
    if (inherits(tea_res, "error")) {
      out_list[[i]] <- tibble(
        .row_id = r$.row_id,
        ok = FALSE,
        err = conditionMessage(tea_res)
      )
      next
    }
    
    pH_end  <- tea_res$pH[nrow(tea_res)]
    DIC_end <- tea_res$DIC_umolkg[nrow(tea_res)]
    
    out_list[[i]] <- tibble(
      .row_id = r$.row_id,
      ok = TRUE,
      err = NA_character_,
      pH_end = pH_end,
      DIC_end = DIC_end,
      dpH_end = pH_end - r$pH_Out,
      dpH_needed = r$pH_In - r$pH_Out,
      dDIC_removed = r$DIC_Out - DIC_end,
      dDIC_needed  = r$DIC_Out - r$DIC_In
    )
  }
  
  res_tbl <- bind_rows(out_list)
  
  # ---- 6) Assemble outputs ----
  msg("[6/6] Assembling outputs ...")
  
  sim_out <- grid_ok %>%
    left_join(res_tbl, by = ".row_id")
  
  failures <- sim_out %>%
    filter(!ok) %>%
    select(.row_id, Site, Date,
           water_mass_kg, kelp_density_gdwkg, kelp_biomass_gdw, exposure_hours,
           Temp_Out, Salt_Out, TA_Out, DIC_Out, pH_Out,
           net_daily_gDW, GP_day_umolC, GP_night_umolC, day_hours, night_hours,
           u_day, u_night,
           err)
  
  if (nrow(failures) > 0) {
    msg("  WARN: ", nrow(failures), " simulations failed. Inspect `$failures`.")
  } else {
    msg("  All simulations succeeded.")
  }
  
  sims <- sim_out %>%
    filter(ok) %>%
    select(Site, Date, DOY,
           water_mass_kg, kelp_density_gdwkg, kelp_biomass_gdw, exposure_hours,
           TA_Out, DIC_Out, pH_Out, DIC_In, pH_In,
           dpH_needed, pH_end, dpH_end,
           dDIC_needed, dDIC_removed,
           net_daily_gDW, frac_day, day_h, night_h, u_day, u_night)
  
  list(
    sims = sims,
    failures = failures,
    grid_ok = grid_ok,
    missing_required = missing_req
  )
}



# --------------- Visualization helpers for landscape slices -------------------

# Show the modelled change in pH due to kelp carbon uptake over the specified 
# exposure time, as a function of water mass and kelp biomass.
# Starts with the outside-bed chemistry for a given site and date. 
# Warmer colors indicate larger pH elevation relative to the outside condition.
plot_landscape_dpH <- function(sim_df, site_pick, date_pick, exposure_pick) {
  df <- sim_df %>% filter(Site == site_pick, Date == date_pick, exposure_hours == exposure_pick)
  
  ggplot(df, aes(x = water_mass_kg, y = kelp_biomass_gdw, fill = dpH_end)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = hcl.colors(100, "RdYlBu", rev = TRUE),
      na.value = "grey90"
    ) +
    labs(
      title = paste0(site_pick, " ", date_pick, " (", exposure_pick, " h)"),
      x = "Water mass (kg)",
      y = "Kelp biomass (g DW)",
      fill = expression(Delta*"pH (Sim - Out)")
    ) +
    theme_classic()
}

# Shows the difference between modelled and observed pH change inside vs outside
# a kelp bed for the same site, date, and exposure time. 
# Values near zero indicate combinations of water mass and kelp biomass for which 
# the kelp-tea model reproduces the observed pH elevation, while positive 
# or negative values indicate over- or under-prediction, respectively.
plot_landscape_match_error <- function(sim_df, site_pick, date_pick, exposure_pick) {
  df <- sim_df %>%
    filter(Site == site_pick, Date == date_pick, exposure_hours == exposure_pick) %>%
    mutate(match_err = dpH_end - dpH_needed)
  
  ggplot(df, aes(x = water_mass_kg, y = kelp_biomass_gdw, fill = match_err)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = hcl.colors(100, "RdYlBu", rev = TRUE),
      na.value = "grey90"
    ) +
    labs(
      title = paste0(site_pick, " ", date_pick, " (", exposure_pick, " h)"),
      x = "Water mass (kg)",
      y = "Kelp biomass (g DW)",
      fill = expression("pH error (model - needed)")
    ) +
    theme_classic()
}


# Shows the modelled change in pH due to kelp carbon uptake over the specified 
# exposure time, as a function of water mass and kelp biomass. AND does this for
# every site, creating a facet plot
plot_landscape_dpH_facets <- function(sim_df,
                                      date_pick,
                                      exposure_pick) {
  
  df <- sim_df %>%
    dplyr::filter(
      Date == date_pick,
      exposure_hours == exposure_pick
    )
  
  if (nrow(df) == 0) {
    stop("No rows found for this Date × exposure_hours combination.")
  }
  
  ggplot(df, aes(
    x = water_mass_kg,
    y = kelp_biomass_gdw,
    fill = dpH_end
  )) +
    geom_tile() +
    facet_wrap(~ Site) +
    scale_fill_gradientn(
      colours = hcl.colors(100, "RdYlBu", rev = TRUE),
      na.value = "grey90",
      name = expression(Delta * "pH (end − out)")
    ) +
    labs(
      title = paste0(
        "Modeled kelp-driven pH change (", exposure_pick, " h exposure)\n",
        "Date: ", date_pick
      ),
      x = "Water mass (kg)",
      y = "Kelp biomass (g DW)"
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1, "lines")
    )
}


#--------- Functions for interpolating simulations and merging data  -----------

# Interpolate a continuous surface in (exposure_hours, log10(density)) space
# for selected site and date.
make_interp_surface <- function(sims, site_pick, date_pick, z = "dpH_end",
                                nx = 200, ny = 200) {
  
  # Filter + prep (interpolate in log-density space to preserve resolution at low rho)
  df <- sims %>%
    dplyr::filter(Site == site_pick, Date == as.Date(date_pick)) %>%
    dplyr::filter(is.finite(kelp_density_gdwkg), kelp_density_gdwkg > 0) %>%
    dplyr::mutate(log10_rho = log10(kelp_density_gdwkg)) %>%
    dplyr::select(exposure_hours, log10_rho, dplyr::all_of(z)) %>%
    dplyr::distinct() %>%
    tidyr::drop_na()
  
  if (nrow(df) < 4) {
    stop("Need at least 4 points to interpolate for this Site × Date (after filtering).")
  }
  
  # Interpolate
  itp <- akima::interp(
    x = df$exposure_hours,
    y = df$log10_rho,
    z = df[[z]],
    nx = nx, ny = ny,
    duplicate = "mean"
  )
  
  # Build surface tibble
  surf <- expand.grid(
    exposure_hours = itp$x,
    log10_rho = itp$y
  )
  
  surf[[z]] <- as.vector(itp$z)
  surf$kelp_density_gdwkg <- 10^(surf$log10_rho)
  
  tibble::as_tibble(surf)
}


# Pull the delta pH observation for selected site and date (ΔpH = In − Out)
make_obs_dpH <- function(sims, site_pick, date_pick) {
  
  x <- sims %>%
    dplyr::filter(Site == site_pick, Date == as.Date(date_pick)) %>%
    dplyr::distinct(dpH_needed) %>%
    dplyr::pull(dpH_needed)
  
  if (length(x) == 0) stop("No dpH_needed found for this Site × Date.")
  if (length(x) > 1) warning("Multiple dpH_needed values found; using the first.")
  
  x[1]
}


# 3) Plot interpolated surface + observed ΔpH contour
plot_interp_density_exposure_dpH <- function(surf,
                                             dpH_obs,
                                             site_pick = NULL,
                                             date_pick = NULL,
                                             z_col = "dpH_end",
                                             log10_density = TRUE) {
  
  # Required columns
  req <- c("exposure_hours", "kelp_density_gdwkg", z_col)
  missing <- setdiff(req, names(surf))
  if (length(missing) > 0) {
    stop("`surf` is missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # dpH_obs sanity
  if (!is.numeric(dpH_obs) || length(dpH_obs) != 1 || is.na(dpH_obs)) {
    stop("`dpH_obs` must be a single, non-NA numeric value (observed ΔpH = In − Out).")
  }
  
  # Filter to finite z (prevents contour/raster warnings)
  surf2 <- surf %>%
    dplyr::filter(is.finite(.data[[z_col]]),
                  is.finite(exposure_hours),
                  is.finite(kelp_density_gdwkg),
                  kelp_density_gdwkg > 0)
  
  if (nrow(surf2) == 0) stop("No finite values in `surf` after filtering.")
  
  # Range check: only draw observed contour if within range
  zmin <- min(surf2[[z_col]], na.rm = TRUE)
  zmax <- max(surf2[[z_col]], na.rm = TRUE)
  has_obs_contour <- (dpH_obs >= zmin && dpH_obs <= zmax)
  
  # Title
  title_txt <- "Interpolated kelp-driven pH response"
  if (!is.null(site_pick) || !is.null(date_pick)) {
    title_txt <- paste0(
      title_txt,
      "\n",
      if (!is.null(site_pick)) paste0("Site: ", site_pick) else "",
      if (!is.null(site_pick) && !is.null(date_pick)) " | " else "",
      if (!is.null(date_pick)) paste0("Date: ", as.Date(date_pick)) else ""
    )
  }
  
  p <- ggplot2::ggplot(
    surf2,
    ggplot2::aes(
      x = exposure_hours,
      y = kelp_density_gdwkg,
      fill = .data[[z_col]]
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colours = hcl.colors(200, "RdYlBu", rev = TRUE),
      na.value = "grey90",
      name = expression(Delta * "pH (end − out)")
    ) +
    ggplot2::labs(
      title = title_txt,
      subtitle = paste0("Dashed = observed ΔpH (In − Out) = ", signif(dpH_obs, 3)),
      x = "Exposure time (h)",
      y = expression("Kelp density (g DW kg"^-1*")")
    ) +
    ggplot2::theme_classic()
  
  if (has_obs_contour) {
    p <- p +
      ggplot2::geom_contour(
        aes(
          x = exposure_hours,
          y = kelp_density_gdwkg,
          z = .data[[z_col]]
        ),
        breaks = dpH_obs,
        colour = "black",
        linetype = "dashed",
        linewidth = 0.6
      )
  }
  
  if (log10_density) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  p
}


# This loops over all available dates for a given site, builds the interpolated
# surface for each date, and attaches the observed ΔpH. It then combines everything 
# into one long data frame and produces one faceted plot (facet by Date)
# This avoids trying to facet separate ggplot objects (which doesn’t work cleanly in ggplot2).
plot_interp_density_exposure_by_date <- function(sims,
                                                 site_pick,
                                                 nx = 200,
                                                 ny = 200,
                                                 z_col = "dpH_end",
                                                 log10_density = TRUE) {
  
  dates <- sims %>%
    dplyr::filter(Site == site_pick) %>%
    dplyr::distinct(Date) %>%
    dplyr::arrange(Date) %>%
    dplyr::pull(Date)
  
  if (length(dates) == 0) stop("No dates available for site: ", site_pick)
  
  surf_list <- lapply(dates, function(d) {
    
    surf <- make_interp_surface(
      sims,
      site_pick = site_pick,
      date_pick = d,
      z = z_col,
      nx = nx,
      ny = ny
    )
    
    dpH_obs <- make_obs_dpH(sims, site_pick = site_pick, date_pick = d)
    
    surf %>%
      dplyr::mutate(Site = site_pick, Date = d, dpH_obs = dpH_obs)
  })
  
  surf_all <- dplyr::bind_rows(surf_list)
  
  # Base plot
  p <- ggplot2::ggplot(
    surf_all,
    ggplot2::aes(
      x = exposure_hours,
      y = kelp_density_gdwkg,
      fill = .data[[z_col]]
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::facet_wrap(~ Date) +
    ggplot2::scale_fill_gradientn(
      colours = hcl.colors(200, "RdYlBu", rev = TRUE),
      na.value = "grey90",
      name = expression(Delta * "pH (end − out)")
    ) +
    ggplot2::labs(
      title = paste0(
        "Interpolated kelp-driven pH response across dates\n",
        "Site: ", site_pick
      ),
      subtitle = "Dashed line = observed ΔpH (In − Out)",
      x = "Exposure time (h)",
      y = expression("Kelp density (g DW kg"^-1*")")
    ) +
    ggplot2::theme_classic()
  
  if (log10_density) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  # Add one contour layer per date at that date's dpH_obs
  for (d in dates) {
    
    s_d <- dplyr::filter(surf_all, Date == d)
    dpH_d <- s_d$dpH_obs[1]
    
    # Skip if dpH_d is NA or outside the surface range (prevents "zero contours" warnings)
    zmin <- min(s_d[[z_col]], na.rm = TRUE)
    zmax <- max(s_d[[z_col]], na.rm = TRUE)
    if (is.na(dpH_d) || dpH_d < zmin || dpH_d > zmax) next
    
    p <- p +
      ggplot2::geom_contour(
        data = s_d,
        inherit.aes = FALSE,
        ggplot2::aes(
          x = exposure_hours,
          y = kelp_density_gdwkg,
          z = .data[[z_col]]
        ),
        breaks = dpH_d,
        colour = "black",
        linetype = "dashed",
        linewidth = 0.75
      )
  }
  
  p
}


# Diagnostic plot produces table of available contours from descrete samples ----
summarise_contour_availability <- function(sims,
                                           site_pick,
                                           nx = 200,
                                           ny = 200,
                                           z_col = "dpH_end") {
  
  dates <- sims %>%
    dplyr::filter(Site == site_pick) %>%
    dplyr::distinct(Date) %>%
    dplyr::arrange(Date) %>%
    dplyr::pull(Date)
  
  if (length(dates) == 0) {
    stop("No dates available for site: ", site_pick)
  }
  
  dplyr::bind_rows(lapply(dates, function(d) {
    
    # Interpolated surface for this date (same as plotting workflow)
    surf <- make_interp_surface(
      sims,
      site_pick = site_pick,
      date_pick = d,
      z = z_col,
      nx = nx,
      ny = ny
    )
    
    # Observed ΔpH for this date
    dpH_obs <- make_obs_dpH(
      sims,
      site_pick = site_pick,
      date_pick = d
    )
    
    # Surface range (finite only)
    zvals <- surf[[z_col]]
    zmin <- min(zvals, na.rm = TRUE)
    zmax <- max(zvals, na.rm = TRUE)
    
    tibble::tibble(
      Site = site_pick,
      Date = as.Date(d),
      dpH_obs = dpH_obs,
      zmin = zmin,
      zmax = zmax,
      has_contour = is.finite(dpH_obs) && dpH_obs >= zmin && dpH_obs <= zmax
    )
  }))
}



#--------------------------- Density-based plot -------------------------------

# A density–time analogue of plot_landscape_dpH() above.
plot_density_exposure_dpH <- function(sims,
                                      site_pick,
                                      date_pick,
                                      log10_density = TRUE) {
  
  df <- sims %>%
    dplyr::filter(
      Site == site_pick,
      Date == date_pick
    )
  
  if (nrow(df) == 0) {
    stop("No simulations found for this Site × Date combination.")
  }
  
  # Observed inside–outside ΔpH (single value per Site × Date)
  dpH_obs <- unique(df$dpH_needed)
  if (length(dpH_obs) != 1 || is.na(dpH_obs)) {
    stop("Expected exactly one observed ΔpH for this Site × Date.")
  }
  
  zmin <- min(df$dpH_end, na.rm = TRUE)
  zmax <- max(df$dpH_end, na.rm = TRUE)
  has_obs_contour <- (dpH_obs >= zmin && dpH_obs <= zmax)
  
  p <- ggplot(
    df,
    aes(
      x = exposure_hours,
      y = kelp_density_gdwkg,
      fill = dpH_end
    )
  ) +
    geom_tile() +
    
    # Observed ΔpH contour
    {
      if (has_obs_contour)
        geom_contour(aes(z = dpH_end),
                     breaks = dpH_obs,
                     colour = "black",
                     linetype = "dashed",
                     linewidth = 0.6)
    } +
    
    scale_fill_gradientn(
      colours = hcl.colors(100, "RdYlBu", rev = TRUE),
      name = expression(Delta * "pH (end − out)"),
      na.value = "grey90"
    ) +
    labs(
      title = paste0(
        "Modeled kelp-driven pH change\n",
        "Site: ", site_pick, " | Date: ", date_pick
      ),
      x = "Exposure time (h)",
      y = expression("Kelp density (g DW kg"^-1*")")
    ) +
    theme_classic()
  
  if (log10_density) {
    p <- p +
      scale_y_log10()
  }
  
  p
}



#------------ Set of functions to plot simulations with in-out data ----------------
#--- MAKES ONE BIG plot ... 

# 1) Build a list of all the interpolated density/exposure surfaces
build_interp_grid_data <- function(sims,
                                   nx = 200,
                                   ny = 200,
                                   z_col = "dpH_end") {
  
  combos <- sims %>%
    dplyr::distinct(Site, Date) %>%
    dplyr::arrange(Date, Site)
  
  if (nrow(combos) == 0) stop("No Site × Date pairs found in `sims`.")
  
  surf_list <- vector("list", nrow(combos))
  
  for (i in seq_len(nrow(combos))) {
    
    s <- combos$Site[i]
    d <- combos$Date[i]
    
    surf <- make_interp_surface(
      sims,
      site_pick = s,
      date_pick = d,
      z = z_col,
      nx = nx,
      ny = ny
    )
    
    dpH_obs <- make_obs_dpH(
      sims,
      site_pick = s,
      date_pick = d
    )
    
    surf_list[[i]] <- surf %>%
      dplyr::mutate(
        Site = s,
        Date = as.Date(d),
        dpH_obs = dpH_obs
      )
  }
  
  dplyr::bind_rows(surf_list)
}

# 2) A wrapper around summarise_contour_availability() to do all sites in sims.
summarise_contour_availability_all_sites <- function(sims,
                                                     nx = 200,
                                                     ny = 200,
                                                     z_col = "dpH_end") {
  
  sites <- sims %>%
    dplyr::distinct(Site) %>%
    dplyr::arrange(Site) %>%
    dplyr::pull(Site)
  
  if (length(sites) == 0) stop("No sites found in `sims`.")
  
  dplyr::bind_rows(lapply(sites, function(s) {
    
    summarise_contour_availability(
      sims      = sims,
      site_pick = s,
      nx        = nx,
      ny        = ny,
      z_col     = z_col
    )
    
  }))
}

# 2) Plot full grid with date- and site-specific observed contours
plot_interp_grid_sites_dates <- function(surfaces,
                                         z_col = "dpH_end",
                                         log10_density = TRUE) {
  
  req <- c("Site", "Date", "exposure_hours", "kelp_density_gdwkg", "dpH_obs", z_col)
  miss <- setdiff(req, names(surfaces))
  if (length(miss) > 0) stop("`surf_all` missing: ", paste(miss, collapse = ", "))
  
  # Ensure consistent facet ordering
  surfaces <- surfaces %>%
    dplyr::mutate(
      Site = factor(Site, levels = sort(unique(Site))),
      Date = factor(Date, levels = sort(unique(Date)))
    )
  
  # Base raster grid
  p <- ggplot2::ggplot(
    surfaces,
    ggplot2::aes(
      x = exposure_hours,
      y = kelp_density_gdwkg,
      fill = .data[[z_col]]
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(Date ~ Site) +
    ggplot2::scale_fill_gradientn(
      colours = hcl.colors(200, "RdYlBu", rev = TRUE),
      na.value = "grey90",
      name = expression(Delta * "pH (end − out)")
    ) +
    ggplot2::labs(
      title = "Interpolated kelp-driven pH response surface",
      subtitle = "Dashed contour = observed ΔpH (In − Out) for that Site × Date",
      x = "Exposure time (h)",
      y = expression("Kelp density (g DW kg"^-1*")")
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.spacing = grid::unit(0.8, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9)
    )
  
  if (log10_density) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  # Add observed contour per facet (Site × Date-specific break)
  combos <- surfaces %>%
    dplyr::distinct(Site, Date, dpH_obs) %>%
    dplyr::arrange(Date, Site)
  
  for (i in seq_len(nrow(combos))) {
    
    s <- combos$Site[i]
    d <- combos$Date[i]
    dpH <- combos$dpH_obs[i]
    
    # Data slice for this facet
    s_df <- surfaces %>%
      dplyr::filter(Site == s, Date == d) %>%
      dplyr::filter(is.finite(.data[[z_col]]),
                    is.finite(kelp_density_gdwkg),
                    kelp_density_gdwkg > 0,
                    is.finite(exposure_hours))
    
    if (nrow(s_df) == 0 || !is.finite(dpH)) next
    
    # Only draw if dpH is within the modeled range (otherwise no contour exists)
    zmin <- min(s_df[[z_col]], na.rm = TRUE)
    zmax <- max(s_df[[z_col]], na.rm = TRUE)
    if (dpH < zmin || dpH > zmax) next
    
    p <- p +
      ggplot2::geom_contour(
        data = s_df,
        inherit.aes = FALSE,
        ggplot2::aes(
          x = exposure_hours,
          y = kelp_density_gdwkg,
          z = .data[[z_col]]
        ),
        breaks = dpH,
        colour = "grey",
        linetype = "dashed",
        linewidth = 0.5
      )
  }
  
  # Add dpH as text to each panel
  p <- p +
    ggplot2::geom_text(
      data = label_df,
      inherit.aes = FALSE,
      ggplot2::aes(
        x = x_lab,
        y = y_lab,
        label = sprintf("ΔpHobs=%.3f", dpH_obs)
      ),
      colour = "white",
      hjust = -0.5,
      vjust = 20,
      size = 3
    )
  #hjust 1=left of plot; 0.5 = middle of plot
  #vjust 1=top of plot; 0=missing; -0.3=missing, 0.3=just above top. Ah. increasing
  #   values move text down today. hjust has a different scale.
  
  return( p )
}


#----
# FIN.
