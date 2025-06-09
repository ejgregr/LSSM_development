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
flow_df   <- read.csv( paste(box_dir, "flow_matrix.csv", sep="/"), stringsAsFactors = FALSE)
volume_df <- read.csv( paste(box_dir, "box_volumes.csv", sep="/"))

# Get box names and number; create an index
boxes <- volume_df$BoxName
n_boxes <- length(boxes)
box_index <- setNames(seq_along(boxes), boxes)

# identify internal and boundary boxes
iboxes <- c("KK", "KK_shore", "M_KI")
bboxes <- c("QCS", "M_JS", "F_KI")

# Assign fixed values for DIC to bboxes
fixed_surf_conc <- c(QCS = 2250, M_JS = 2200, F_KI = 1900)
fixed_bot_conc  <- c(QCS = 2280, M_JS = 2300, F_KI = 2000)

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

# Add vertical mixing within boxes
alpha_mix <- c(
  KK        = 0.1,   # weak stratification
  KK_shore  = .9,   # shallow, well-mixed
  M_KI      = 0.02,  # turbulent, moderate mixing
  M_JS      = 0.02,  # turbulent, moderate mixing
  QCS       = 0.005, # deep, well stratified
  F_KI      = 0.005  # well stratified, density-based
)

# Initialise concentration arrays
n_steps <- 100  # total tidal cycles (e.g., 50 days)

# Initialise tracer matrices (DIC example)
conc_surf <- matrix(NA, nrow = n_boxes, ncol = n_steps, dimnames = list(boxes, NULL))
conc_bot  <- matrix(NA, nrow = n_boxes, ncol = n_steps, dimnames = list(boxes, NULL))

# Example initial DIC conditions (µmol/kg)
conc_surf[, 1] <- c(2100, 2000, 2150, 2200, 2250, 1900)[match(boxes, boxes)]
conc_bot[, 1]  <- c(2200, 2000, 2180, 2300, 2280, 2000)[match(boxes, boxes)]


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
    
  # for (j in boxes) {
  #   if (vol_surf[j] > 0 && vol_bot[j] > 0) {
  #     
  #     # Lookup box-specific alpha
  #     alpha <- alpha_mix[j]
  # 
  #     # Masses in each layer
  #     mass_s <- conc_surf[j, t] * vol_surf[j]
  #     mass_b <- conc_bot[j, t]  * vol_bot[j]
  # 
  #     # Exchange between layers
  #     delta_s_to_b <- alpha * mass_s
  #     delta_b_to_s <- alpha * mass_b
  # 
  #     # Update masses
  #     mass_s_new <- mass_s - delta_s_to_b + delta_b_to_s
  #     mass_b_new <- mass_b - delta_b_to_s + delta_s_to_b
  # 
  #     # Convert to concentration
  #     conc_surf[j, t] <- mass_s_new / vol_surf[j]
  #     conc_bot[j, t]  <- mass_b_new / vol_bot[j]
  #   }
  # }
    
    conc_surf[j, t] <- mix$conc_surf
    conc_bot[j, t]  <- mix$conc_bot
  }
  
  
  
  # Reset boundary conditions
  conc_surf[bboxes, t] <- fixed_surf_conc[bboxes]
  conc_bot[bboxes, t]  <- fixed_bot_conc[bboxes]
}


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


### Fin
