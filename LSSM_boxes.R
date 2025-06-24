#################################################################################
# Script:  LSSM_boxes.R
# Created: June 3, 2025. EJG
# 
# Working with ChatGPT to produce a proof-of-concept model for simulating the 
# movement of DIC from kelp forests to the broader Broughton region. 

# ## Updates: 
# 2025/06/09: A working 12-box tracer model now does vertical exchange between surface and deep waters, 
#   and horizontal exchange between Klukuane, its nearshore, and its 3 main connections to the outside.
#   It also allows for freshwater effects from Knight Inlet. 
#   Runs on estimates of volume, tracer (e.g., DIC) concentrations, and flow matrices between the boxes.

# TO DO:
# - Adjust initial DIC values so season is stable (Rather than approach stability from low values)
# - Include vertical mixing so there is some DIC in the bottom layer
# - Add role of kelp in DIC drawdown
# - Adjust boundary boxes for (daily? weekly?) change in DIC? 
#
#################################################################################

setwd( "c:/Data/Git/LSSM_development")
source( "LSSM_box_functions.R" )

#---- Constants and global variables ----
box_dir <- "C:/Data/Git/LSSM_development/Boxes"

# Load configuration files (flow matrix and box volumes)
volume_df  <- read.csv( paste(box_dir, "box_volumes.csv", sep="/"), stringsAsFactors = FALSE)
flow_df    <- read.csv( paste(box_dir, "flow_matrix_active.csv", sep="/"), stringsAsFactors = FALSE)
#climate_df <- read.csv( paste(box_dir, "box_climate.csv", sep="/"), stringsAsFactors = FALSE)

# Get box names and details; create an index
boxes <- volume_df$BoxName
n_boxes <- length(boxes)
box_index <- setNames(seq_along(boxes), boxes)

# identify internal and boundary boxes
iboxes <- boxes[ volume_df$BoxType =="internal" ]
bboxes <- boxes[ volume_df$BoxType =="boundary" ]

# Create volume vectors 
vol_surf <- setNames(volume_df$SurfaceVolume_m3, volume_df$BoxName)[boxes]
vol_bot  <- setNames(volume_df$BottomVolume_m3, volume_df$BoxName)[boxes]

# Initialize flow matrices
ebb_surf <- matrix(0, n_boxes, n_boxes, dimnames = list(boxes, boxes))
ebb_bot  <- matrix(0, n_boxes, n_boxes, dimnames = list(boxes, boxes))
flood_surf <- matrix(0, n_boxes, n_boxes, dimnames = list(boxes, boxes))
flood_bot  <- matrix(0, n_boxes, n_boxes, dimnames = list(boxes, boxes))

x <- populateFlowMatrices()
ebb_surf   <- x$ebb_surf
ebb_bot    <- x$ebb_bot
flood_surf <- x$flood_surf
flood_bot  <- x$flood_bot

# Need some time. 
n_steps <- 150  # total tidal cycles (e.g., ~50 days)


# Load vertical mixing constants for boxes
alpha_mix <- setNames(volume_df$VerticalMixing, volume_df$BoxName)


#---- Initialise DIC concentration arrays ----

# Assign fixed values for DIC to boundary boxes
# Drop 2 oceanic boundaries by 100 ... 
fixed_surf_conc <- c(QCS = 2100, JS_M = 2000, KI_F = 1800)
fixed_bot_conc  <- c(QCS = 2150, JS_M = 2100, KI_F = 1900)



# Initialise tracer matrices (DIC example)
conc_surf <- matrix(NA, nrow = n_boxes, ncol = n_steps, dimnames = list(boxes, NULL))
conc_bot  <- matrix(NA, nrow = n_boxes, ncol = n_steps, dimnames = list(boxes, NULL))

### DIC HERE FOR BBOXES SHOULD BE SET TO FIXED VALUES ABOVE
# Example initial DIC conditions (µmol/kg)
# Box order:      GK_Shore, GK, KI_M, KI_F, JS_M, QCS
conc_surf[, 1] <- c(1900, 1950, 2000, 1800, 2100, 2200 )[match(boxes, boxes)]
conc_bot[, 1]  <- c(1920, 1970, 2020, 1900, 2200, 2250 )[match(boxes, boxes)]

# Grow some kelp for the KK_shore box
kelp_params <- data.frame(
  N_plants = 1000000,
  B_init = 0.025,
  B_max = 9.23,
  R_max = 0.065,
  wet_to_dry = 0.13,
  dry_to_C = 0.25
)

kelp_uptake_vector <- generate_kelp_DIC_uptake(kelp_params, day_stamps)
dim(kelp_uptake_vector )
head(kelp_uptake_vector)
#plot(kelp_uptake$B_plant)

cum_DIC_loss <- cumsum(kelp_uptake_vector$DIC_uptake_mol)
sum( cum_DIC_loss )

###########################
#---- Simulation loop ----

for (t in 2:n_steps) {
  # Choose ebb or flood
  if (t %% 2 == 0) {
    flow_surf <- ebb_surf
    flow_bot  <- ebb_bot
  } else {
    flow_surf <- flood_surf
    flow_bot  <- flood_bot
  }
  
  # Surface layer update
  for (j in boxes) {
    inflow_surf <- sum(flow_surf[, j] * conc_surf[, t - 1] * vol_surf)
    outflow_surf <- sum(flow_surf[j, ] * conc_surf[j, t - 1] * vol_surf[j])
    conc_surf[j, t] <- (conc_surf[j, t - 1] * vol_surf[j] + inflow_surf - outflow_surf) / vol_surf[j]
  }
  
  # Bottom layer update
  for (j in boxes) {
    inflow_bot <- sum(flow_bot[, j] * conc_bot[, t - 1] * vol_bot)
    outflow_bot <- sum(flow_bot[j, ] * conc_bot[j, t - 1] * vol_bot[j])
    conc_bot[j, t] <- (conc_bot[j, t - 1] * vol_bot[j] + inflow_bot - outflow_bot) / vol_bot[j]
  }
  
  # Vertical mixing after tracer advection
  
  for (j in boxes) {
    mix <- applyVerticalMixing(
      box = j,
      conc_surf = conc_surf[, t],
      conc_bot = conc_bot[, t],
      vol_surf = vol_surf,
      vol_bot = vol_bot,
      alpha_mix = alpha_mix
    )
    conc_surf[j, t] <- mix$conc_surf
    conc_bot[j, t]  <- mix$conc_bot
  }
  
  # Kelp drawdown in KK_shore during update
  if ("GK_shore" %in% boxes) {
    # mol → µmol/kg
    loss_umol_kg <- kelp_uptake_vector$DIC_uptake_mol[t] * 1e6 / vol_surf["GK_shore"]
    conc_surf["GK_shore", t] <- conc_surf["GK_shore", t] - loss_umol_kg
  }
  
  # Reset boundary conditions
  conc_surf[bboxes, t] <- fixed_surf_conc[bboxes]
  conc_bot[bboxes, t]  <- fixed_bot_conc[bboxes]
  
}

conc_surf

#----- Troubleshooting
# sum(flood_surf)
# sum(ebb_surf)
# sum(flood_bot)
# sum(ebb_bot)
# 
# flow_df
# 
# cat("vol_bot[M_KI] =", vol_bot["M_KI"], "\n")
# sum(flow_bot[, "M_KI"] * conc_bot[, 1] * vol_bot)


#---- Plot and inspect results ----
# Combine in plotting order
all_boxes <- c(iboxes, bboxes)
line_types <- ifelse(all_boxes %in% iboxes, 1, 2)
colors <- 1:length(all_boxes)

# Plot surface layer
matplot(t(conc_surf[all_boxes, ]), type = 'l', lty = line_types, col = colors,
        main = "Surface Layer Tracer Concentration", xlab = "Time step", ylab = "µmol/kg")
legend("topright", legend = all_boxes, col = colors, lty = line_types, bty = "n")

# Plot bottom layer
matplot(t(conc_bot[all_boxes, ]), type = 'l', lty = line_types, col = colors,
        main = "Bottom Layer Tracer Concentration", xlab = "Time step", ylab = "µmol/kg")
legend("topright", legend = all_boxes, col = colors, lty = line_types, bty = "n")

# plot( kelp_uptake_vector$DIC_uptake_mol )
### Fin
