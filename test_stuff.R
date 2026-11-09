
# Single density, high exposure resolution for ceiling test
# Important to fix density for testing

### SIM CALL!!
sim_ceil <- simulate_dens_landscape(
  pairs, kelp_df, env_daily,
  kelp_density_plants_m2_vec = 7.66,           # W&P mean
  exposure_hours_vec         = seq(1, 24, 1),
  verbose = TRUE )



# Ceiling test (all Site × Date pairs)
ceiling_test <- function(sims, tol = 0.016) {
  ct <- sims %>%
    dplyr::group_by(Site, Date) %>%
    dplyr::summarise(dpH_obs = dpH_needed[1],
                     ceiling = max(dpH_end, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(exceeds = dpH_obs > ceiling + tol)
  
  cat("Pairs:", nrow(ct), "\n",
      "Observed above ceiling + tol:", sum(ct$exceeds), "\n",
      "Spearman (ceiling vs observed):",
      round(cor(ct$ceiling, ct$dpH_obs, method = "spearman"), 2), "\n")
  ct
}

ct <- ceiling_test(sim_ceil$sims, tol = 0.016)


plot_ceiling_test <- function(ct, tol = 0.016) {
  plot(ct$ceiling, ct$dpH_obs, pch = 16, col = ifelse(ct$exceeds, "red", "grey40"),
       xlab = expression("Closed-box ceiling "*Delta*"pH (7.66 plants m"^-2*")"),
       ylab = expression("Observed "*Delta*"pH (In - Out)"))
  abline(0, 1)                          # observed = ceiling
  abline(tol, 1, lty = 2)               # ceiling + tolerance
  abline(h = 0, lty = 3)
  legend("topleft", bty = "n", pch = c(16, 16), col = c("grey40", "red"),
         legend = c("within ceiling + tol", "above ceiling + tol"))
}
plot_ceiling_test(ct)



# Back to the ceiling test after adjusting chem_data df for Hakai discrepancy ... 

pairs <- make_pairs_out_in(chem_data)

### SIM CALL!!
sim_ceil <- simulate_dens_landscape(pairs, kelp_df, env_daily,
                                    kelp_density_plants_m2_vec = 7.66,
                                    exposure_hours_vec         = seq(1, 24, 1))
ct <- ceiling_test(sim_ceil$sims, tol = 0.016)
plot_ceiling_test(ct)



b <- chem_data[as.Date(chem_data$DateTime, tz = bc_time) != as.Date("2025-09-15"), ]

d <- b$pH - b$pH_lab # departure from the 1:1 line

plot(b$pH_lab, d, pch = 16)
abline(h = 0, lty = 2)

tapply(d, grp, function(x) c(n = length(x), bias = mean(x), sd = sd(x), rmse = sqrt(mean(x^2))))

t.test(d ~ grp)          # does the average departure (bias) differ between the sets?
fligner.test(d ~ grp)    # does the scatter around the line differ?






# Plots recovered DIC (from pCO2_lab and pH_lab) minus the lab's adjusted TCO2, against TCO2_lab.
# A check on add_DIC_Alk(): points on the dashed line mean the recovered and lab values agree.
# df needs T_lab, S_lab, pCO2_lab, pH_lab, TCO2_lab, DateTime. Optional highlight_date marks one day in red.
plot_DIC_check <- function(df, highlight_date = NULL) {
  rec  <- add_DIC_Alk(df)                    # works on a copy; df is not modified
  dDIC <- rec$DIC - rec$TCO2_lab
  
  col <- rep("grey40", nrow(rec))
  if (!is.null(highlight_date))
    col[as.Date(rec$DateTime, tz = bc_time) == as.Date(highlight_date)] <- "red"
  
  plot(rec$TCO2_lab, dDIC, pch = 16, col = col, bty = "l",
       xlab = "Lab Adj. TCO2 (umol/kg)", ylab = "Recovered DIC - lab TCO2 (umol/kg)")
  abline(h = 0, lty = 2)
  if (!is.null(highlight_date))
    legend("bottomright", bty = "n", pch = 16, col = c("grey40", "red"),
           legend = c("other dates", format(as.Date(highlight_date))))
  invisible(data.frame(Station = rec$Station, DateTime = rec$DateTime, dDIC = dDIC))
}
a <- plot_DIC_check(chem_data)                              # basic

# Plots recovered DIC (from pCO2_lab and pH_lab) minus the lab's adjusted TCO2, against TCO2_lab.
# Points are coloured by the three patterns we saw:
#   1) a batch problem (batch_date), 2) a systematic offset at high TCO2, 3) a small offset at low TCO2.
# df needs T_lab, S_lab, pCO2_lab, pH_lab, TCO2_lab, DateTime.
plot_DIC_check <- function(df, batch_date = "2025-09-15", split = 1500) {
  rec  <- add_DIC_Alk(df)                    # works on a copy; df is not modified
  dDIC <- rec$DIC - rec$TCO2_lab
  
  is_batch <- as.Date(rec$DateTime, tz = bc_time) == as.Date(batch_date)
  grp <- ifelse(is_batch, "batch",
                ifelse(rec$TCO2_lab >= split, "high", "low"))
  cols <- c(batch = "#D55E00", high = "#0072B2", low = "#009E73")
  labs <- c(batch = paste("Batch", format(as.Date(batch_date))),
            high  = paste0("High TCO2 offset (>= ", split, ")"),
            low   = paste0("Low TCO2 (< ", split, ")"))
  
  plot(rec$TCO2_lab, dDIC, pch = 16, col = cols[grp], bty = "l",
       xlab = "Lab Adj. TCO2 (umol/kg)", ylab = "Recovered DIC - lab TCO2 (umol/kg)")
  abline(h = 0, lty = 2)
  legend("topleft", bty = "n", pch = 16, col = cols, legend = labs)
  invisible(data.frame(Station = rec$Station, DateTime = rec$DateTime,
                       dDIC = dDIC, group = grp))
}
a <- plot_DIC_check(chem_data)

# Plots recovered DIC (from pCO2_lab and pH_lab) minus the lab's adjusted TCO2, against TCO2_lab.
# Two reference lines: y = 0 (recovered agrees with lab) and a fitted sloped line (suspected
# constants / pH-scale difference). Each point goes to the nearer line; points more than `tol`
# from both are "neither". The sloped line is seeded from the high-TCO2 band (>= seed_min).
# df needs T_lab, S_lab, pCO2_lab, pH_lab, TCO2_lab, DateTime.
plot_DIC_check <- function(df, tol = 10, seed_min = 1650) {
  rec <- add_DIC_Alk(df)                      # works on a copy; df is not modified
  d   <- data.frame(x = rec$TCO2_lab, y = rec$DIC - rec$TCO2_lab)
  
  # Fit the sloped line, refitting on the points nearest to it (3 passes)
  keep <- d$x >= seed_min & d$y > 0
  for (i in 1:3) {
    fit  <- lm(y ~ x, data = d[keep, ])
    pred <- predict(fit, d)
    keep <- abs(d$y - pred) < tol & abs(d$y - pred) < abs(d$y)
  }
  message(sprintf("Sloped line: dDIC = %.1f + %.4f x TCO2 (zero at %.0f umol/kg), n = %d",
                  coef(fit)[1], coef(fit)[2], -coef(fit)[1] / coef(fit)[2], sum(keep)))
  
  d_s <- abs(d$y - pred); d_0 <- abs(d$y)
  grp <- ifelse(pmin(d_s, d_0) > tol, "neither", ifelse(d_s < d_0, "slope", "zero"))
  
  cols <- c(slope = "#0072B2", zero = "#009E73", neither = "#D55E00")
  labs <- c(slope   = "Near sloped line (constants / pH scale?)",
            zero    = "Near y = 0 (agrees with lab)",
            neither = "Far from both lines")
  n <- table(factor(grp, levels = names(cols)))
  
  plot(d$x, d$y, pch = 16, col = cols[grp], bty = "l",
       xlab = "Lab Adj. TCO2 (umol/kg)", ylab = "Recovered DIC - lab TCO2 (umol/kg)")
  abline(h = 0, lty = 2)
  abline(fit, lty = 2, col = cols["slope"])
  legend("topleft", bty = "n", pch = 16, col = cols, legend = paste0(labs, " (n = ", n, ")"))
  
  invisible(data.frame(Station = rec$Station, DateTime = rec$DateTime,
                       TCO2_lab = d$x, dDIC = d$y, group = grp))
}

# Plots recovered DIC (from pCO2_lab and pH_lab) minus the lab's adjusted TCO2, against TCO2_lab.
# Step 1: samples from batch_date are split off first (all of that day).
# Step 2: the rest are grouped by proximity to a fitted sloped line (suspected constants /
# pH-scale difference) or to y = 0 (agrees with lab); points more than `tol` from both are "neither".
# The sloped line is seeded from the high-TCO2 band (>= seed_min), excluding the batch.
# df needs T_lab, S_lab, pCO2_lab, pH_lab, TCO2_lab, DateTime.
plot_DIC_check <- function(df, batch_date = "2025-09-15", tol = 10, seed_min = 1650) {
  rec <- add_DIC_Alk(df)                      # works on a copy; df is not modified
  d   <- data.frame(x = rec$TCO2_lab, y = rec$DIC - rec$TCO2_lab)
  is_batch <- as.Date(rec$DateTime, tz = bc_time) == as.Date(batch_date)
  
  # Fit the sloped line on non-batch points, refitting on the points nearest to it (3 passes)
  keep <- !is_batch & d$x >= seed_min & d$y > 0
  for (i in 1:3) {
    fit  <- lm(y ~ x, data = d[keep, ])
    pred <- predict(fit, d)
    keep <- !is_batch & abs(d$y - pred) < tol & abs(d$y - pred) < abs(d$y)
  }
  message(sprintf("Sloped line: dDIC = %.1f + %.4f x TCO2 (zero at %.0f umol/kg), n = %d",
                  coef(fit)[1], coef(fit)[2], -coef(fit)[1] / coef(fit)[2], sum(keep)))
  
  d_s <- abs(d$y - pred); d_0 <- abs(d$y)
  grp <- ifelse(is_batch, "batch",
                ifelse(pmin(d_s, d_0) > tol, "neither",
                       ifelse(d_s < d_0, "slope", "zero")))
  
  cols <- c(batch = "#D55E00", slope = "#0072B2", zero = "#009E73", neither = "grey")
  labs <- c(batch   = paste("Batch", format(as.Date(batch_date))),
            slope   = "Near sloped line (constants / pH scale?)",
            zero    = "Near y = 0 (agrees with lab)",
            neither = "Other, far from both lines")
  n <- table(factor(grp, levels = names(cols)))
  
  plot(d$x, d$y, pch = 16, col = cols[grp], bty = "l",
       xlab = "Lab Adj. TCO2 (umol/kg)", ylab = "Recovered DIC - lab TCO2 (umol/kg)")
  abline(h = 0, lty = 2)
  abline(fit, lty = 2, col = cols["slope"])
  legend("topleft", bty = "n", pch = 16, col = cols, legend = paste0(labs, " (n = ", n, ")"))
  
  invisible(data.frame(Station = rec$Station, DateTime = rec$DateTime,
                       TCO2_lab = d$x, dDIC = d$y, group = grp))
}

plot_DIC_check(chem_data)


#------------------------------------------------------------------------

#----  Testing here to ensure our chemistry matches Hakai's for the Co2Pro processing ----
chem_data$dDIC <- chem_data$DIC - chem_data$TCO2_lab
chem_data$dTA  <- chem_data$TA  - chem_data$TA_lab

summary(chem_data[, c("dDIC", "dTA")])

# The plot shows 3 different patterns in the deviance including a parameter problem, 
# a conversion problem, and a regular bias. 
#---> NEED parameters used by Hakai, review of "batch of 15" by Wiley.
plot(chem_data$TCO2_lab, chem_data$dDIC, pch = 16, col = "grey40",
     xlab = "Lab Adj. TCO2 (umol/kg)", ylab = "Recovered DIC - lab TCO2")
abline(h = 0, lty = 2)

# Check consistency of Hakai calculations by comparing lab reported pH to 
# pH calcd from lab TCO2 & TA (using seacarb).  
chk <- seacarb::carb(flag = 15, var1 = chem_data$TA_lab / 1e6,
                     var2 = chem_data$TCO2_lab / 1e6,
                     S = chem_data$Salinity, T = chem_data$Temp,
                     Patm = 1, P = 1, Pt = 0, Sit = 0)
chem_data$dpH_chk <- chk$pH - chem_data$pH
summary(chem_data$dpH_chk)
bad <- chem_data[abs(chem_data$dpH_chk) > 0.02, c("Station", "DateTime", "TCO2_lab", "dDIC", "dpH_chk")]
bad

# Bad batch sent to Wiley - 2026/11/05.
write.csv(bad, file.path(results_dir, "bad_26Nov06.csv"), row.names = FALSE)







