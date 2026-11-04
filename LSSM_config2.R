#################################################################################
# Script:  LSSM_config2.R
# Created: January 13 2026. 
# Purpose: Support second iteration of the nereo growth and chemistry work
#          Includes packages, paths, and functions.
# Notes:
#  - Early development supported by ChatGPT AI
#  - Feb 16: Added Chemistry functions. Lots of funky plotting. 
#----------------------------------------------------------------------------
# Updates:
# 26Oct18: 
# - Removed all data loading functions (these are now in LSSM_Water_Analysis.proj)
#  26Nov: 
# - Major revisions as part of re-factoring the code. This time with Claude AI
#################################################################################

#---- Load require packages ----
# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "lubridate", "dplyr", "stringr",
                        "rmarkdown","knitr", "tinytex", "kableExtra", #currently causing trouble.
                        "seacarb", "gsw", "deSolve", "tidyr", "readxl", "akima" )

# Other packages that might be useful. 
# "tidyr", "raster", "stringr", "rasterVis",
# "diffeR", "vegan", "ranger", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
# "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr", "binr", 'gwxtab'

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)
version$version.string

#---- Configuration ----
bc_time <- 'America/Vancouver'

data_dir <- "C:/Data/Git/LSSM_prototype/Data"
PDF_dir  <- "C:/Data/Git/LSSM_prototype/PDF" # for Markdown script
source_dir <- "C:/Data/Git/LSSM_water_analysis/Results"


#-------------------- SUPPORT FUNCTIONS ----------------------------------
#           Nothing should be created below here 
# NOTE: Discrete data now loaded in water_analysis.proj. 

#-------------------------- DEB Support Functions ---------------------------- 
# Functions to scale T and DLI effects as per Pontier et al.

# Warm-side Gaussian penalty with a cold-side plateau. Simple representation of
# the warm temperature limitation identified by Pontier et al. 
fT_reparam <- function(T_C, T_opt = 10.0, sigma_warm = 1.4) {
  ifelse(
    T_C <= T_opt,
    1.0,
    exp(- (T_C - T_opt)^2 / (2 * sigma_warm^2))
  )
}

# --- DLI scaling for blades: hump-shaped with low- and high-light penalties ---
# Approx ~23% gain per +10 DLI near ~20 (low-light side) and ~23% loss per +10 above ~40 (high-light side)
fL_reparam <- function(DLI, L_low = 20, L_high = 40, 
                       low_gain_per_10 = 0.23,   # interpreted as "being 10 below optimal costs ~23%"
                       high_loss_per_10 = 0.23) {
  
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

# --- Helper for DLI scaling: convert "p% change per unit" into exponential rate constant ---
# If multiplier per unit is m (e.g., 0.77 per +1°C), then k = -log(m) / unit
k_from_multiplier <- function(multiplier, unit = 1) {
  -log(multiplier) / unit
}


#---------------------------- DEB Plotting Functions ---------------------------

# Updated to ensure that DW->WW conversion is timely and updated x axis label with date
plot_kelp_biomass_WW <- function(out_df,
                                 wet_to_dry = model_params$wet_to_dry,
                                 time_col = "real_DOY",
                                 months = TRUE,
                                 bth_col  = "B_th",
                                 bfr_col  = "B_fr",
                                 year = 2025,          # sets month start days (leap-year safe)
                                 log_y = FALSE) {
  stopifnot(all(c(time_col, bth_col, bfr_col) %in% names(out_df)))
  
  th <- out_df[[bth_col]] / wet_to_dry   # g DW -> g WW
  fr <- out_df[[bfr_col]] / wet_to_dry
  
  df_long <- data.frame(
    time    = rep(out_df[[time_col]], 3),
    biomass = c(fr, th, fr + th),
    series  = factor(rep(c("Frond", "Stipe", "Total"), each = nrow(out_df)),
                     levels = c("Frond", "Stipe", "Total"))
  )
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = biomass, color = series)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(y = "Biomass (g WW)", color = NULL) +
    ggplot2::theme_classic()
  
  if (time_col == "real_DOY" && months == TRUE) {
    # DOY of the 1st of each month; keep those inside the data range
    month_start <- as.numeric(format(as.Date(sprintf("%d-%02d-01", year, 1:12)), "%j"))
    rng  <- range(out_df[[time_col]], na.rm = TRUE)
    keep <- month_start >= rng[1] & month_start <= rng[2]
    p <- p +
      ggplot2::scale_x_continuous(breaks = month_start[keep], labels = month.abb[keep]) +
      ggplot2::labs(x = NULL)
  } else {
    p <- p + ggplot2::labs(x = "Time (days since start)")
  }
  
  if (log_y) p <- p + ggplot2::scale_y_log10()
  p
}

# Plot Daily carbon budget: gross production partitioned into DOC leakage, the two
# respiration terms and the net carbon available for frond growth. The stack
# height equals GP. Sloughing is drawn as a separate line because it is removed
# from standing frond biomass, not from that day's GP.
# Daily carbon budget per kelp: gross production partitioned into DOC leakage,
# the two respiration terms and the net carbon available for frond growth.
# The stack height equals GP. Sloughing is drawn as a separate line because it
# is removed from standing frond biomass, not from that day's GP.
plot_C_budget <- function(df, time_col = "real_DOY", months = TRUE, year = 2025,
                          prop = FALSE,        # TRUE = fraction of GP instead of mmol C
                          slough = TRUE) {     # draw sloughing line (ignored if prop = TRUE)
  req <- c(time_col, "GP_umolC", "DOC_umolC", "R_prod_umolC", "R_maint_umolC", "Slough_umolC")
  missing <- setdiff(req, names(df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))
  
  to_mmol <- 1 / 1000                       # umol C -> mmol C
  net <- df$GP_umolC - df$DOC_umolC - df$R_prod_umolC - df$R_maint_umolC
  if (any(net < 0, na.rm = TRUE)) {
    warning("Net carbon < 0 on ", sum(net < 0, na.rm = TRUE),
            " days; negative values are shown as zero, so the stack is shorter than GP on those days.")
    net <- pmax(net, 0)
  }
  
  lev <- c("DOC leakage", "Production respiration", "Maintenance respiration", "Net carbon gain")
  long <- data.frame(
    day   = rep(df[[time_col]], 4),
    mmolC = c(df$DOC_umolC, df$R_prod_umolC, df$R_maint_umolC, net) * to_mmol,
    flux  = factor(rep(lev, each = nrow(df)), levels = lev)   # first level plots on top
  )
  
  cols <- c("DOC leakage" = "#E69F00", "Production respiration" = "#56B4E9",
            "Maintenance respiration" = "#0072B2", "Net carbon gain" = "#009E73")
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = day, y = mmolC, fill = flux)) +
    ggplot2::geom_area(position = if (prop) "fill" else "stack", alpha = 0.9) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(y = if (prop) "Fraction of gross production"
                  else expression("mmol C kelp"^-1*" d"^-1),
                  fill = NULL) +
    ggplot2::theme_classic()
  
  if (slough && !prop) {
    sl <- data.frame(day = df[[time_col]], mmolC = df$Slough_umolC * to_mmol)
    p <- p + ggplot2::geom_line(data = sl,
                                ggplot2::aes(x = day, y = mmolC, linetype = "Sloughing (as C)"),
                                inherit.aes = FALSE, colour = "black", linewidth = 0.8) +
      ggplot2::scale_linetype_manual(values = c("Sloughing (as C)" = "dashed")) +
      ggplot2::labs(linetype = NULL)
  }
  
  if (time_col == "real_DOY" && months) {
    month_start <- as.numeric(format(as.Date(sprintf("%d-%02d-01", year, 1:12)), "%j"))
    rng  <- range(df[[time_col]], na.rm = TRUE)
    keep <- month_start >= rng[1] & month_start <= rng[2]
    p <- p + ggplot2::scale_x_continuous(breaks = month_start[keep], labels = month.abb[keep],
                                         expand = c(0, 0)) +
      ggplot2::labs(x = NULL)
  } else {
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::labs(x = "Day of year")
  }
  
    # Legend inside the plot, top left
  inside <- utils::packageVersion("ggplot2") >= "3.5.0"
  p <- p + ggplot2::theme(
    legend.justification = c(0, 1),
    legend.background    = ggplot2::element_blank(),
    legend.key.size      = ggplot2::unit(0.4, "cm"),
    legend.spacing.y     = ggplot2::unit(0.1, "cm"))
  p <- if (inside) {
    p + ggplot2::theme(legend.position = "inside",
                       legend.position.inside = c(0.02, 0.98))
  } else {
    p + ggplot2::theme(legend.position = c(0.02, 0.98))
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

# Stipe-frond allometry: the two-piece stipe target (Eq. 17), the Stekoll regression it
# is built from, and the model's own trajectory. Plotted in g WW; log-log by default
# so the small-plant branch is visible.
# Stipe-frond allometry: the two-piece stipe target (Eq. 17), the Stekoll regression it
# is built from, and the model's own trajectory. Plotted in g WW; log-log by default
# so the small-plant branch is visible.
plot_allometry <- function(kelp_df = NULL,
                           pars = model_params,
                           log = "xy") {
  wd  <- pars$wet_to_dry
  a   <- pars$th_slope                 # slope (ratio, unit-free)
  b   <- pars$th_int / wd              # intercept, g WW (105.2)
  Bx  <- pars$th_B_cross / wd          # g WW (50): limit of Stekoll data
  rho <- pars$th_max_ratio             # 2.43
  
  lo <- if (grepl("x", log)) 1 else 0
  hi <- 1200
  x  <- if (lo > 0) exp(seq(log(lo), log(hi), length.out = 400)) else seq(lo, hi, length.out = 400)
  reg    <- a * x + b                  # Stekoll regression, extended below Bx
  target <- pmin(rho * x, a * x + b)   # model stipe target
  ylo <- if (grepl("y", log)) 1 else 0
  
  op <- par(mar = c(4.5, 4.5, 1, 1) + 0.1); on.exit(par(op))
  plot(x, target, type = "n", log = log, xlim = c(lo, hi),
       ylim = c(max(ylo, min(target)), max(c(reg, target)) * 1.05),
       xlab = "Frond biomass (g WW)", ylab = "Stipe and bulb biomass (g WW)", bty = "l")
  u <- par("usr")
  xl <- if (grepl("x", log)) 10^u[1] else u[1]
  yl <- if (grepl("y", log)) 10^u[3] else u[3]
  yh <- if (grepl("y", log)) 10^u[4] else u[4]
  
  # region below Bx: no data behind the regression
  rect(xl, yl, Bx, yh, col = adjustcolor("grey90", 0.5), border = NA)
  abline(v = Bx, lty = 2, col = "grey40")
  
  lines(x, reg, lty = 3, lwd = 2, col = "grey40")        # regression, extrapolated
  lines(x, target, lwd = 3, col = "black")               # model target
  
  if (!is.null(kelp_df)) {                               # model trajectory
    fr <- kelp_df$B_fr / wd; th <- kelp_df$B_th / wd
    lines(fr, th, col = "steelblue", lwd = 2)
    points(fr[1], th[1], pch = 21, bg = "white", col = "steelblue", cex = 1.6, lwd = 2)
    n <- length(fr)
    points(fr[n], th[n], pch = 16, col = "steelblue", cex = 1.6)
    text(fr[1], th[1], "start", pos = 4, col = "steelblue", cex = 0.85)
    text(fr[n], th[n], "end of run", pos = 1, col = "steelblue", cex = 0.85)
  }
  
  # label for B_x, hung below the top edge (adj, not pos, so it isn't clipped)
  text(Bx, yh, bquote(B[x] == .(round(Bx)) ~ "g WW"),
       adj = c(-0.1, 1.5), cex = 0.8, col = "grey30")
  text(xl, yh, expression(atop("no data below " * B[x], "(proportional branch assumed)")),
       adj = c(-0.05, 1.3), cex = 0.8, col = "grey30")
  
  labs <- list(bquote("Proportional branch: " * B[th] == .(signif(rho, 3)) ~ B[fr]),
               bquote("Stekoll et al. 2006 line: " * B[th] == .(a) ~ B[fr] + .(round(b, 1)) ~ "g WW"),
               bquote("Regression extended below " * B[x]))
  lty <- c(1, 1, 3); lwd <- c(3, 3, 2); col <- c("black", "black", "grey40")
  if (!is.null(kelp_df)) {
    labs <- c(labs, list("Model trajectory"))
    lty <- c(lty, 1); lwd <- c(lwd, 2); col <- c(col, "steelblue")
  }
  legend("bottomright", bty = "n", cex = 0.85, legend = as.expression(labs),
         lty = lty, lwd = lwd, col = col)
}

# Daily temperature and fT over the season, using the model's parameters.
plot_temperature_scaling <- function(env_dat,
                                     T_opt = model_params$T_opt,
                                     sigma_warm = model_params$sigma_warm) {
  DOY      <- env_dat$real_DOY
  T_daily  <- env_dat$Temp
  fT_daily <- fT_reparam(T_daily, T_opt = T_opt, sigma_warm = sigma_warm)
  
  op <- par(mar = c(5, 4, 4, 4) + 0.1); on.exit(par(op))
  plot(DOY, T_daily, type = "l", col = "steelblue", lwd = 2,
       xlab = "Day of Year", ylab = "Temperature (°C)", bty = "l")
  par(new = TRUE)
  plot(DOY, fT_daily, type = "l", col = "firebrick", lwd = 2,
       axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
  axis(4, at = seq(0, 1, by = 0.2), col = "firebrick", col.axis = "firebrick")
  mtext("Temperature scaling factor (fT)", side = 4, line = 3, col = "firebrick")
}

# Daily DLI and fL over the season. Applies the same surface-reflectance 
# adjustment the model applies before calling fL_reparam().
plot_DLI_scaling <- function(env_dat, alpha_water = model_params$alpha_water) {
  DOY      <- env_dat$real_DOY
  DLI      <- (1 - alpha_water) * env_dat$DLI
  fL_daily <- fL_reparam(DLI)
  
  op <- par(mar = c(5, 4, 4, 4) + 0.1); on.exit(par(op))
  plot(DOY, DLI, type = "l", col = "steelblue", lwd = 2,
       xlab = "Day of Year", ylab = "DLI after reflectance (mol m⁻² d⁻¹)", bty = "l")
  par(new = TRUE)
  plot(DOY, fL_daily, type = "l", col = "firebrick", lwd = 2,
       axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
  axis(4, at = seq(0, 1, by = 0.2), col = "firebrick", col.axis = "firebrick")
  mtext("DLI scaling factor (fL)", side = 4, line = 3, col = "firebrick")
}

# Parameterless functions to plot the environmental curves.
plot_fT_curve <- function(T_opt = model_params$T_opt,
                          sigma_warm = model_params$sigma_warm,
                          T_min = 5, T_max = 20, by = 0.1) {
  
  T_seq <- seq(T_min, T_max, by = by)
  df <- data.frame(T_C = T_seq,
                   fT  = fT_reparam(T_seq, T_opt = T_opt, sigma_warm = sigma_warm))
  
  ggplot2::ggplot(df, ggplot2::aes(x = T_C, y = fT)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = T_opt, linetype = "dashed") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = expression("Temperature ("*degree*C*")"),
                  y = expression(f[T])) +
    ggplot2::theme_classic()
}

# Light scaling curve. Defaults come from fL_reparam()'s own formals, so the 
# plot shows exactly what the model uses.
plot_fL_curve <- function(L_low  = formals(fL_reparam)$L_low,
                          L_high = formals(fL_reparam)$L_high,
                          low_gain_per_10  = formals(fL_reparam)$low_gain_per_10,
                          high_loss_per_10 = formals(fL_reparam)$high_loss_per_10,
                          DLI_min = 0, DLI_max = 70, DLI_step = 0.5) {
  
  df <- data.frame(DLI = seq(DLI_min, DLI_max, by = DLI_step))
  df$fL <- fL_reparam(df$DLI, L_low = L_low, L_high = L_high,
                      low_gain_per_10 = low_gain_per_10,
                      high_loss_per_10 = high_loss_per_10)
  
  ggplot2::ggplot(df, ggplot2::aes(x = DLI, y = fL)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = c(L_low, L_high), linetype = "dashed") +
    ggplot2::labs(x = expression("Daily Light Integral (mol m"^-2*" d"^-1*")"),
                  y = expression(f[L])) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_classic()
}


# ----------------- Water-sample pairing and preparation  -----------------------

# Combines Inside / Outside water samples by Site & Date, giving a single row 
# per paired event w clearly labeled _In / _Out carbonate and physical variables.
# Returns one row per Site x Date, with variables suffixed _In / _Out. 
# IO and Site come from Station (e.g. CP_In); Drops stations without _In/_Out suffix.
make_pairs_out_in <- function(chem_water) {
  chem_water %>%
    mutate(IO   = case_when(grepl("_in$",  Station, ignore.case = TRUE) ~ "In",
                            grepl("_out$", Station, ignore.case = TRUE) ~ "Out"),
           Site = sub("_(in|out)$", "", Station, ignore.case = TRUE),
           Date = as.Date(DateTime, tz = bc_time) ) %>%
    filter(!is.na(IO)) %>%
    average_replicates() %>%
    pivot_wider(id_cols = c(Site, Date), names_from = IO,
                values_from = c(DateTime, Temp, Salinity, TA, pH, pCO2, DIC, OArag),
                names_sep = "_") %>%
    drop_na()
}

# Joins paired water samples to DEB day & env_daily day/night hrs. Select by DOY
# pulled from water sample. DEB provides kelp carbon uptake rates and day/night 
# photoperiod info to each sampling event.
merge_pairs_with_deb_days <- function(pairs_df,
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

# Helper function for make_pairs_out_in()
# Average replicate samples (same Site x Date x IO), printing the records involved.
average_replicates <- function(d) {
  reps <- d %>%
    group_by(Site, Date, IO) %>% filter(n() > 1) %>% ungroup() %>%
    arrange(Site, Date, IO, DateTime) %>%
    select(Site, Date, IO, DateTime, Temp, Salinity, TA, pH, pCO2, DIC)
  if (nrow(reps) > 0) {
    message("Averaging ", nrow(reps), " replicate records:")
    print(as.data.frame(reps), row.names = FALSE)
  }
  d %>%
    group_by(Site, Date, IO) %>%
    summarise(across(c(DateTime, Temp, Salinity, TA, pH, pCO2, DIC, OArag), mean), .groups = "drop")
}

# NOTE DLI is already daily. Temperature is averaged.
make_env_daily <- function(dli_df, temp_df) {
  # 1. Process temperature: Convert to Date and calculate daily mean
  temp_daily <- temp_df %>%
    mutate(Datestamp = as.Date(Datestamp, tz = bc_time)) %>%
    group_by(Datestamp) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE), .groups = "drop")
  
  # 2. Process DLI: Ensure Datestamp is Date type
  dli_daily <- dli_df %>%
    mutate(Datestamp = as.Date(Datestamp, tz = bc_time))
  
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
    mutate(Datestamp = as.Date(Datestamp, tz = bc_time)) %>%
    # Group by each day
    group_by(Datestamp) %>%
    # Count hours where PAR is greater than 0
    summarize(
      day_hours   = sum(PAR > 0, na.rm = TRUE),
      night_hours = sum(PAR <= 0, na.rm = TRUE),
      .groups = "drop"
    )
}


# ----------------- Carbonate chemistry Utilities -----------------------
# A wrapper to filter incomplete rows and append DIC to the full discrete 
# water dataset without aggregating by site or condition.
# Adds total alkalinity (TA) and DIC, both in umol/kg, calculated from pCO2 and pH.
add_DIC_Alk <- function(df) {
  ok <- complete.cases(df[, c("Temp", "Salinity", "pCO2", "pH")])
  if (any(!ok)) message("Dropped ", sum(!ok), " rows with missing T, S, pCO2 or pH.")
  df <- df[ok, ]
  
  res <- seacarb::carb(flag = 21, var1 = df$pCO2, var2 = df$pH,   # pCO2 (uatm), pH (total)
                       S = df$Salinity, T = df$Temp, Patm = 1, P = 1, Pt = 0, Sit = 0)
  df$TA  <- res$ALK * 1e6    # seacarb returns mol/kg
  df$DIC <- res$DIC * 1e6
  df
}

# TA vs salinity: scatter, linear fit with 95% CI band, and residuals vs salinity.
# df needs Salinity and TA (umol/kg). Returns the lm object.
fit_TA_salinity <- function(df) {
  df <- df[complete.cases(df[, c("Salinity", "TA")]), ]
  m  <- lm(TA ~ Salinity, data = df)
  print(summary(m)); print(confint(m))
  
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1)); on.exit(par(op))
  
  sx <- seq(min(df$Salinity), max(df$Salinity), length.out = 100)
  pr <- predict(m, data.frame(Salinity = sx), interval = "confidence")
  plot(df$Salinity, df$TA, pch = 16, col = adjustcolor("steelblue", 0.6), bty = "l",
       xlab = "Salinity", ylab = expression("TA ("*mu*"mol kg"^-1*")"))
  polygon(c(sx, rev(sx)), c(pr[, 2], rev(pr[, 3])), col = adjustcolor("grey", 0.4), border = NA)
  lines(sx, pr[, 1], lwd = 2)
  s <- summary(m)
  legend("topleft", bty = "n", cex = 0.85,
         legend = c(sprintf("TA = %.1f + %.1f \u00D7 S", coef(m)[1], coef(m)[2]),
                    sprintf("R\u00B2 = %.2f,  n = %d,  residual SD = %.1f", s$r.squared, nrow(df), s$sigma)))
  
  plot(df$Salinity, resid(m), pch = 16, col = adjustcolor("steelblue", 0.6), bty = "l",
       xlab = "Salinity", ylab = expression("Residual TA ("*mu*"mol kg"^-1*")"))
  abline(h = 0, lty = 2)
  invisible(m)
}

# Computes pH from total alkalinity and DIC (holding TA constant) using seacarb::carb().
# Translates modelled DIC changes into pH responses.
ph_from_ta_dic <- function(TA_umolkg, DIC_umolkg, S, T_C, P_dbar = 1) {
  # You confirmed this flag is correct in your setup.
  res <- seacarb::carb(
    flag = 15,            # TA + DIC -> pH (per your confirmation)
    var1 = TA_umolkg / 1e6,
    var2 = DIC_umolkg /1e6,
    S    = S,
    T    = T_C,
    Patm = 1,
    P    = P_dbar,
    Pt   = 0,
    Sit  = 0
  )
  res$pH
}


#----- Tea box model -----
# Simulates the time evolution of DIC and pH in a closed, well-mixed water parcel 
# exposed to kelp at a given density (plants per kg water) for up to 24 h.
# Assumptions: TA constant, no gas exchange or flushing
# Starting ocean conditions are the "Outside" samples. Returns DIC (umol/kg) and pH at each step.
run_tea <- function(init_TA_umolkg, init_DIC_umolkg, S, T_C,
                    kelp_density_plants_kg,
                    exposure_hours,
                    u_day_umolC_plant_h, u_night_umolC_plant_h,
                    start_datetime, day_hours,
                    lon = -126.6,
                    utc_offset = as.numeric(format(start_datetime, "%z")) / 100) {
  
  stopifnot(exposure_hours <= 24)
  
  # Hourly steps, plus a final partial step if exposure is not a whole number of hours
  times_h <- c(0, seq_len(floor(exposure_hours)))
  if (tail(times_h, 1) < exposure_hours) times_h <- c(times_h, exposure_hours)
  
  # Daylight window: day_hours (from the DEB) centred on solar noon, in clock hours
  noon <- 12 + (15 * utc_offset - lon) / 15
  win  <- noon + c(-1, 1) * day_hours / 2
  # Fraction of the clock interval [a, b] in daylight (shifts handle runs past midnight)
  day_frac <- function(a, b) {
    sum(sapply(c(-24, 0, 24), function(k)
      max(0, min(b, win[2] + k) - max(a, win[1] + k)))) / (b - a)
  }
  start_clock <- lubridate::hour(start_datetime) + lubridate::minute(start_datetime) / 60
  
  DIC <- numeric(length(times_h)); pH <- numeric(length(times_h))
  DIC[1] <- init_DIC_umolkg
  pH[1]  <- ph_from_ta_dic(init_TA_umolkg, DIC[1], S, T_C)
  
  for (i in 2:length(times_h)) {
    dt <- times_h[i] - times_h[i - 1]
    # Weight day and night rates by the daylight fraction of this step
    f <- day_frac(start_clock + times_h[i - 1], start_clock + times_h[i])
    u <- f * u_day_umolC_plant_h + (1 - f) * u_night_umolC_plant_h   # umol C plant-1 h-1
    # Per-plant drawdown scaled by plants per kg of water (negative u at night raises DIC)
    DIC[i] <- DIC[i - 1] - u * kelp_density_plants_kg * dt           # umol kg-1
    # TA is held constant, so pH follows from the new DIC
    pH[i]  <- ph_from_ta_dic(init_TA_umolkg, DIC[i], S, T_C)
  }
  
  tibble(t_h = times_h, DateTime = start_datetime + hours(times_h),
         TA_umolkg = init_TA_umolkg, DIC_umolkg = DIC, pH = pH)
}


# --------------------- landscape simulations ------------------------
# Runs the kelp-tea box model across a Cartesian grid of kelp densities 
# (kelp biomass per water mass) and exposure times for each valid Site × Date pairing.
# Returns modeled pH and DIC changes, along with diagnostics for debugging,
# including: successful runs ($sims), failed runs ($failures), and every 
# attempted run with inputs ($grid_ok).
simulate_dens_landscape <- function(pairs_df, kelp_df, env_daily,
                                    kelp_density_plants_m2_vec,
                                    exposure_hours_vec,
                                    box_area_m2  = 1,
                                    box_depth_m  = 10,
                                    seawater_kg_m3 = 1025,
                                    verbose = TRUE) {
  
  msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # ---- 1) Join pairs to DEB fluxes and photoperiod (pairs outside the DEB run are dropped) ----
  msg("[1/4] Joining paired water samples to DEB + env_daily ...")
  pairs2 <- merge_pairs_with_deb_days(pairs_df, kelp_df, env_daily)
  
  req <- c("TA_Out", "DIC_Out", "Salinity_Out", "Temp_Out", "DateTime_Out",
           "pH_Out", "pH_In", "DIC_In",
           "GP_day_umolC", "GP_night_umolC", "R_prod_umolC", "R_maint_umolC",
           "day_hours", "night_hours")
  n0 <- nrow(pairs2)
  pairs2 <- pairs2 %>%
    filter(!if_any(all_of(req), is.na), day_hours > 0, night_hours > 0)
  msg("  Kept ", nrow(pairs2), " / ", n0, " paired events with complete inputs.")
  if (nrow(pairs2) == 0) stop("No valid paired events after join + filtering.")
  
  # ---- 2) Build grid (paired Site x Date x kelp density x exposure time) ----
  # With one row per simulation containing water chemistry, DEB fluxes, 
  # photo-period and per-plant uptake rates passed to run_tea(). The version 
  # returned by run_tea() has the results attached (pH_end, DIC_end, dpH_end, ok/err).
  msg("[2/4] Building grid (pairs x density x exposure) and per-plant uptake rates ...")
  
  water_mass_kg <- box_area_m2 * box_depth_m * seawater_kg_m3   # 10,250 kg for the default box
  grid_ok <- pairs2 %>%
    tidyr::crossing(kelp_density_plants_m2 = kelp_density_plants_m2_vec,
                    exposure_hours     = exposure_hours_vec) %>%
    mutate(
      .row_id = row_number(),
      kelp_density_plants_kg = kelp_density_plants_m2 * box_area_m2 / water_mass_kg,
      
      # DIC drawdown = GP - respiration; DOC counts as removed. Production respiration is
      # daytime only; maintenance respiration is split by hours of day and night.
      maint_day   = R_maint_umolC * day_hours / (day_hours + night_hours),
      maint_night = R_maint_umolC - maint_day,
      u_day       = (GP_day_umolC - R_prod_umolC - maint_day) / day_hours,     # umol C plant-1 h-1
      u_night     = (GP_night_umolC - maint_night) / night_hours
    )
  msg("  Simulations to run: ", nrow(grid_ok))
  
  # ---- 3) Run the tea box for each row, capturing errors ----
  msg("[3/4] Running tea simulations ...")
  n_sim <- nrow(grid_ok)
  out_list <- vector("list", n_sim)
  
  for (i in seq_len(n_sim)) {
    if (verbose && (i == 1 || i %% 250 == 0 || i == n_sim)) msg("  Tea run ", i, " / ", n_sim)
    r <- grid_ok[i, ]
    
    tea_res <- tryCatch(
      run_tea(init_TA_umolkg = r$TA_Out, init_DIC_umolkg = r$DIC_Out,
              S = r$Salinity_Out, T_C = r$Temp_Out,
              kelp_density_plants_kg = r$kelp_density_plants_kg,
              exposure_hours = r$exposure_hours,
              u_day_umolC_plant_h = r$u_day, u_night_umolC_plant_h = r$u_night,
              start_datetime = r$DateTime_Out, day_hours = r$day_hours),
      error = function(e) e)
    
    if (inherits(tea_res, "error")) {
      out_list[[i]] <- tibble(.row_id = r$.row_id, ok = FALSE, err = conditionMessage(tea_res))
      next
    }
    pH_end  <- tea_res$pH[nrow(tea_res)]
    DIC_end <- tea_res$DIC_umolkg[nrow(tea_res)]
    out_list[[i]] <- tibble(
      .row_id = r$.row_id, ok = TRUE, err = NA_character_,
      pH_end = pH_end, DIC_end = DIC_end,
      dpH_end = pH_end - r$pH_Out,       dpH_needed  = r$pH_In - r$pH_Out,
      dDIC_removed = r$DIC_Out - DIC_end, dDIC_needed = r$DIC_Out - r$DIC_In
    )
  }
  
  # ---- 4) Assemble outputs ----
  msg("[4/4] Assembling outputs ...")
  sim_out <- grid_ok %>% left_join(bind_rows(out_list), by = ".row_id")
  
  failures <- sim_out %>% filter(!ok) %>%
    select(.row_id, Site, Date, kelp_density_plants_m2, exposure_hours,
           Temp_Out, Salinity_Out, TA_Out, DIC_Out, u_day, u_night, err)
  msg(if (nrow(failures) > 0) paste0("  WARN: ", nrow(failures), " simulations failed; see $failures.")
      else "  All simulations succeeded.")
  
  sims <- sim_out %>% filter(ok) %>%
    select(Site, Date, DOY, kelp_density_plants_m2, exposure_hours,
           TA_Out, DIC_Out, pH_Out, DIC_In, pH_In,
           dpH_needed, pH_end, dpH_end, dDIC_needed, dDIC_removed,
           day_hours, night_hours, u_day, u_night)
  
  list(sims = sims, failures = failures, grid_ok = sim_out)
}

make_dpH_df <- function(sim_out) {
  
  sim_out$sims %>%
    dplyr::mutate(
      Date = as.Date(Date),
      dpH_end    = pH_end - pH_Out,   # modeled change
      dpH_needed = pH_In  - pH_Out,   # observed change
      match_error = dpH_end - dpH_needed
    ) %>%
    dplyr::select(
      Site,
      Date,
      dpH_end,
      dpH_needed,
      match_error
    )
}

make_dpH_summary_df <- function(sim_out) {
  
  sim_out$sims %>%
    dplyr::mutate(
      Date = as.Date(Date),
      dpH_end    = pH_end - pH_Out,   # modeled ΔpH
      dpH_needed = pH_In  - pH_Out,   # observed ΔpH
      match_error = dpH_end - dpH_needed
    ) %>%
    dplyr::group_by(Site, Date) %>%
    dplyr::slice_min(
      order_by = abs(match_error),
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Site,
      Date,
      dpH_needed,
      dpH_end,
      match_error,
      exposure_hours
    )
}

# Time ends of a data distribution. For removing outliers
filter_outliers <- function(df, max_abs) {
  df %>%
    dplyr::filter(abs(match_error) <= max_abs)
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
    dplyr::filter(is.finite(kelp_density_plants_m2), kelp_density_plants_m2 > 0) %>%
    dplyr::mutate(log10_rho = log10(kelp_density_plants_m2)) %>%
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
  surf$kelp_density_plants_m2 <- 10^(surf$log10_rho)
  
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

# Plot interpolated surface + observed ΔpH contour
plot_interp_density_exposure_dpH <- function(surf,
                                             dpH_obs,
                                             site_pick = NULL,
                                             date_pick = NULL,
                                             z_col = "dpH_end",
                                             log10_density = TRUE) {
  
  # Required columns
  req <- c("exposure_hours", "kelp_density_plants_m2", z_col)
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
                  is.finite(kelp_density_plants_m2),
                  kelp_density_plants_m2 > 0)
  
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
      y = kelp_density_plants_m2,
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
      y = expression("Kelp density (plants m"^-2*")")
    ) +
    ggplot2::theme_classic()
  
  if (has_obs_contour) {
    p <- p +
      ggplot2::geom_contour(
        aes(
          x = exposure_hours,
          y = kelp_density_plants_m2,
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
      y = kelp_density_plants_m2,
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
      y = expression("Kelp density (plants m"^-2*")")
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
          y = kelp_density_plants_m2,
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

# Diagnostic plot produces table of available contours from discrete samples ----
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
      y = kelp_density_plants_m2,
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
      y = expression("Kelp density (plants m"^-2*")")
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
        Date = as.Date(d, tz="America/Vancouver"),
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
  
  req <- c("Site", "Date", "exposure_hours", "kelp_density_plants_m2", "dpH_obs", z_col)
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
      y = kelp_density_plants_m2,
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
      y = expression("Kelp density (plants m"^-2*")")
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
                    is.finite(kelp_density_plants_m2),
                    kelp_density_plants_m2 > 0,
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
          y = kelp_density_plants_m2,
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
