library(dplyr)
library(tidyr)
library(ggplot2)
library(metR)  # for text on contours




#------------------------------------------------

library(dplyr)
library(ggplot2)

# Ensure Date is Date and orderable
df_err <- sim_out$sims %>%
  dplyr::filter(!is.na(match_error)) %>%
  dplyr::mutate(Date = as.Date(Date))

ggplot(df_err, aes(x = Site, y = match_error)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  
  # distribution
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    fill = "grey90",
    colour = "grey30"
  ) +
  
  # chronological points
  geom_jitter(
    aes(colour = Date),
    width = 0.15,
    height = 0,
    size = 2,
    alpha = 0.9
  ) +
  
  scale_colour_viridis_c(
    option = "plasma",
    name = "Date"
  ) +
  
  labs(
    x = "Site",
    y = expression("Match error (modeled − observed "*Delta*"pH)"),
    title = "Distribution of model–observation pH mismatch",
    subtitle = "Points colored by sampling date; dashed line = perfect match"
  ) +
  
  theme_classic()



#-----------------------------------------------

# ---- Assumptions ----
dw_frac <- 0.12     # gDW / gWW
rho_sw  <- 1.025    # kg seawater / L

# ---- Ranges ----
wet_kg_vec <- seq(0.1, 10, length.out = 200)
vol_L_vec  <- seq(1, 300, length.out = 200)

# ---- Compute density ----
df <- expand_grid(
  kelp_wet_kg = wet_kg_vec,
  water_L = vol_L_vec
) %>%
  mutate(
    kelp_dw_g  = kelp_wet_kg * 1000 * dw_frac,
    water_kg   = water_L * rho_sw,
    dens_gdwkg = kelp_dw_g / water_kg
  )

# ---- Log-spaced contour levels ----
levels_log <- 10^seq(
  log10(min(df$dens_gdwkg)),
  log10(max(df$dens_gdwkg)),
  length.out = 8
)

# or just use specific values ... 
levels_log <- c(1, 5, 10, 50, 100)


x <- geom_contour(breaks = levels_log, linewidth = 0.6)


# ---- Contour plot (log space) ----
p <- ggplot(df, aes(x = water_L, y = kelp_wet_kg, z = dens_gdwkg)) +
  geom_contour(breaks = levels_log, linewidth = 0.6) +
  metR::geom_text_contour(
    breaks = levels_log,
    aes(label = signif(after_stat(level), 2)),
    size = 3,
    stroke = 0.2,
    label.placer = metR::label_placer_flattest(),
    min.size = 0
#    label.placer = metR::label_placer_n(n = 1),
#    check_overlap = FALSE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Water volume (L)",
    y = "Kelp wet weight (kg WW)",
    title = "Contours of equal kelp density (log-spaced)",
    subtitle = expression("Density = gDW kg"^-1*" seawater")
  ) +
  theme_classic()

print(p)
