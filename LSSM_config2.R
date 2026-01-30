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


#----
# FIN.
