# Code provided by ChatGPT Feb02 2026. Digestion in progress. 


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
# Assumptions:
#    - No gas exchange, Well mixed, TA constant
#    - DIC updated by u_net(t)*biomass/water_mass
#---------------------- Tea box forward model  ------------------------------

# Classifies a timestamp as day or night based on clock hour. Used to switch 
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


# --------------- landscape simulations and diagnostics ----------------------

# Runs the kelp-tea box model across a Cartesian grid of water mass, kelp biomass, 
# and exposure time for each valid Site × Date pairing, returning modeled pH changes, 
# DIC changes, and diagnosticswhile capturing and reporting failures for debugging.
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




# =============================================================================
# 8) Example workflow (adapt object names to yours)
# =============================================================================

# (A) Add DIC to discrete water samples
inout_dic <- add_DIC_no_grouping(inout_data)

# (B) Pair In/Out by Site + Date

# NB: CP IN on May is a replicate and an odd duck in the pairs DF.
pairs <- make_pairs_out_in(inout_dic)

# (C) Run response landscape (requires kelp_df + env_daily)
sim_out <- simulate_landscape(
   pairs_df = pairs,
   kelp_df = kelp_df,
   env_daily = env_daily,
   water_mass_kg_vec = c(10, 25, 50),
   kelp_biomass_gdw_vec = c(5, 10, 15),
   exposure_hours_vec = c(6, 12, 24)
#   water_mass_kg_vec = c(5, 10, 20, 30, 40, 50, 60, 80, 100),
#   kelp_biomass_gdw_vec = c(5, 10, 20, 30, 40, 50),
#   exposure_hours_vec = c(2, 6, 12, 24)
 )

#--- Check out simulation components 
pairs2 <- prepare_pairs_with_deb_env(pairs, kelp_df, env_daily)
# Combines DEB growth day with the sampling day
# Fist CP paired sample has no growth assoc. Collected before moorings in?

avail_plots <-pairs2 %>% distinct(Site, Date) %>% arrange(Site, Date)

# Available dates by site
x$sims %>%
  dplyr::filter(Site == 'CP') %>%
  dplyr::distinct(Date) %>%
  dplyr::arrange(Date)


# (Da) Pick a slice
site_pick <- "CP"
date_pick <- as.Date("2025-07-22")
exposure_pick <- 12

df_slice <- x$sims %>%
  dplyr::filter(Site == site_pick, Date == date_pick, exposure_hours == exposure_pick)

nrow(df_slice)


# (Da) Plot a slice
plot_landscape_dpH(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22"), exposure_pick = 2)
plot_landscape_dpH(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22"), exposure_pick = 12)
plot_landscape_dpH(sim_out$sims, site_pick = "CP", date_pick = as.Date("2025-07-22"), exposure_pick = 24)


# Diagnostics 

summary(sim_out$sims$net_daily_gDW)
summary(sim_out$sims$u_day)
summary(sim_out$sims$u_night)

sim_out$sims %>%
  summarise(
    frac_net_positive = mean(net_daily_gDW > 0),
    frac_dDIC_pos     = mean(dDIC_removed > 0),
    frac_dpH_pos      = mean(dpH_end > 0)
  )





#---- FIN. ----
