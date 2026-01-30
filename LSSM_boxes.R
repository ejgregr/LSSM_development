#################################################################################
# Script:  LSSM_boxes.R
# Created: June 3, 2025. EJG
# 
# Working with ChatGPT to produce a proof-of-concept model for simulating the 
# movement of DIC from kelp forests to the broader Broughton region. 

# ## Updates: 
# 2025/06/09: A working 12-box tracer model (6 horizontal x 2 vertical) now does 
#   vertical exchange between surface and deep waters, and horizontal exchange between 
#   Village Sea, its nearshore, and its 3 main connections to the outside.
#   It also allows for freshwater effects from Knight Inlet. 
#   Runs on estimates of volume, tracer (e.g., DIC) concentrations, and flow matrices between the boxes.

# TO DO:
# - Adjust initial DIC values so season is stable (Rather than approach stability from low values)
# - Include vertical mixing so there is some DIC in the bottom layer
# - Add role of kelp in DIC drawdown
# - Adjust boundary boxes for (daily? weekly?) change in DIC? 
#
#################################################################################

#---- Constants and global variables ----
setwd( "c:/Data/Git/LSSM_development")

box_dir <- "C:/Data/Git/LSSM_development/Boxes"
DEB_dir  <- "C:/Data/Git/LSSM_development/DEB"

source( "LSSM_box_functions.R" )

#---- Prepare the box model ----
# Load configuration files (flow matrix and box volumes)
volume_df  <- read.csv( paste(box_dir, "box_volumes.csv", sep="/"), stringsAsFactors = FALSE)
#flow_df    <- read.csv( paste(box_dir, "flow_matrix_active.csv", sep="/"), stringsAsFactors = FALSE)
flow_df    <- read.csv( paste(box_dir, "flow_matrix_Nov16_almost.csv", sep="/"), stringsAsFactors = FALSE)
#climate_df <- read.csv( paste(box_dir, "box_climate.csv", sep="/"), stringsAsFactors = FALSE)

# Trim the input data (cuz notes in csv file)
volume_df <- volume_df[,-ncol(volume_df)]

# Get box names and details; create an index
boxes <- volume_df$BoxName
n_boxes <- length(boxes)
box_index <- setNames(seq_along(boxes), boxes)

# identify internal and boundary boxes
iboxes <- boxes[ volume_df$BoxType =="internal" ]
bboxes <- boxes[ volume_df$BoxType =="boundary" ]

# Create volume vectors 
vol_surf <- setNames(volume_df$SurfVol_m3, volume_df$BoxName)[boxes]
vol_bot  <- setNames(volume_df$BottVol_m3, volume_df$BoxName)[boxes]

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

# Load vertical mixing constants for boxes
alpha_mix <- setNames(volume_df$VerticalMix, volume_df$BoxName)
alpha_mix <- setNames( c(0, 0, 0, 0, 0, 0), volume_df$BoxName )

#---- Load kelp growth ---- 
# Use saved kelp growth time series (see LSSM_DEB.r, LSSM_growth.R vestigial but still some good pieces?)
kelp_grow <- read.csv(paste0( DEB_dir, "/kelp_grow_Nov15.csv") )

kelp_uptake <- generate_kelp_DIC_uptake(kelp_grow, 10000)
head( kelp_uptake)
# For longer 'burn in' time series when kelp growth doesn't matter
#kelp_uptake <- generate_kelp_DIC_uptake(kelp_params, c(kelp_grow$days,kelp_grow$days,kelp_grow$days) )

n_steps <- dim( kelp_uptake )[[1]] 

#---- Initialise DIC concentration arrays ----
# Assign fixed values for DIC to boundary boxes
# fixed_surf_conc <- c(QCS = 1900, KI_F = 1900)
fixed_bot_conc  <- c(QCS = 1920, KI_F = 1920)

# Initialise tracer matrices (DIC example)
conc_surf <- matrix(NA, nrow = n_boxes, ncol = n_steps, dimnames = list(boxes, NULL))
conc_bot  <- matrix(NA, nrow = n_boxes, ncol = n_steps, dimnames = list(boxes, NULL))

### DIC HERE FOR BBOXES SHOULD BE SET TO FIXED VALUES ABOVE
# Example initial DIC conditions (µmol/kg)

fixed_surf_conc <- c(QCS = 2000, KI_F = 1800)
fixed_bot_conc  <- c(QCS = 2100, KI_F = 1900)
# Box order:      VS_Shore, CP, VS, KI_M, KI_F, QCS
conc_surf[, 1] <- c(1900, 1900, 1900, 1900, fixed_surf_conc[2], fixed_surf_conc[1] )[match(boxes, boxes)]
conc_bot[, 1]  <- c(1920, 1920, 1920, 1920, fixed_bot_conc[2],  fixed_bot_conc[1] )[match(boxes, boxes)]

# Update boxes with balanced values (i.e., conc_surf after 'burn in')
conc_surf[, 1] <- c(4163, 2706, 2155, 365, fixed_surf_conc[2], fixed_surf_conc[1] )[match(boxes, boxes)]
conc_bot[, 1]  <- c(4163, 2706, 2155, 365, fixed_bot_conc[2],  fixed_bot_conc[1] )[match(boxes, boxes)]

cum_DIC_loss <- cumsum(kelp_uptake$DIC_uptake_mol)
sum( cum_DIC_loss )

#--------------------------- Simulation loop -------------------------------
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
  
  # Kelp drawdown in VS_shore during update
  kelp_DIC <- as.numeric( kelp_uptake$DIC_uptake_mol[t] )
  if ("VS_shore" %in% boxes) {
    # mol → µmol/kg
    loss_umol_kg <- kelp_DIC * 1e6 / vol_surf["VS_shore"]
    conc_surf["VS_shore", t] <- conc_surf["VS_shore", t] - loss_umol_kg
  }

  # Reset boundary conditions
  conc_surf[bboxes, t] <- fixed_surf_conc[bboxes]
  conc_bot[bboxes, t]  <- fixed_bot_conc[bboxes]
  
}
#--- End model loop 

#---- Plot and inspect results ----
#---- DIC concentrations All boxes ----
par(cex=1.2)
all_boxes <- c(iboxes, bboxes)
line_types <- ifelse(all_boxes %in% iboxes, 1, 2)
colors <- 1:length(all_boxes)

# Plot surface layer
matplot(t(conc_surf[all_boxes, ]), type = 'l', lty = line_types, lwd = 2, col = colors,
        main = "Surface Layer Tracer Concentration", xlab = "Time step", ylab = "µmol/kg")
legend("topright", legend = all_boxes, col = colors, lty = line_types, , lwd = 2, bty = "n")

# Plot bottom layer
#matplot(t(conc_bot[all_boxes, ]), type = 'l', lty = line_types, col = colors,
#        main = "Bottom Layer Tracer Concentration", xlab = "Time step", ylab = "µmol/kg")
#legend("topright", legend = all_boxes, col = colors, lty = line_types, bty = "n")

#---- DIC concentrations of interior boxes ----
par(cex = 1.2) # For clipping plots to ppt ... 
interior <- setdiff(boxes, bboxes)
matplot(
  t(conc_surf[interior, ]),
  type = "l", lwd = 2,
  lty = 1,                     # <-- force all lines to be solid
  xlab = "Time step",
  ylab = "Surface DIC (µmol/kg)",
  main = "Interior Surface Concentrations",
  col = 1:length(interior)
)
legend(
  "topright",
  legend = interior,
  col = 1:length(interior),
  lwd = 2,
  lty = 1                     # <-- solid in legend
)

#---- DIC concentrations of boundary boxes ----
matplot(
  t(conc_surf[bboxes, ]),
  type = "l", lwd = 2,
  lty = 1,                     # <-- force all lines to be solid
  xlab = "Time step",
  ylab = "Surface DIC (µmol/kg)",
  main = "Boundary Box Surface Concentrations",
  col = 1:length(bboxes)
)
legend(
  "topright", legend = bboxes, col = 1:length(bboxes), lwd = 2
)

#---- Total DIC over time ----
total_DIC <- colSums(conc_surf * vol_surf)
plot(
  total_DIC,
  type = "l", lwd = 2,
  lty = 1,                     # <-- force all lines to be solid
  xlab = "Time step",
  ylab = "Total Surface DIC (µmol/kg × m³)",
  main = "Total Surface DIC Over Time"
)


#---- testing drawndown effect on concentrations: each needs a separate run of the model  ----
# not sure this is needed ... :\
conc_no_kelp <- conc_surf
conc_1000_kelp <- conc_surf
conc_100000_kelp <- conc_surf

#saved runs
sum( conc_no_kelp['VS_shore',100:110] )
sum( conc_1000_kelp['VS_shore',100:110] )
sum( conc_100000_kelp['VS_shore',100:110] )

plot_box_series("VS_shore", conc_no_kelp, conc_1000_kelp, conc_100000_kelp)

plot_box_series <- function(box_name, c1, c2, c3) {
  
  # extract the row for the chosen box
  x <- rbind(
    c1[box_name, ],
    c2[box_name, ],
    c3[box_name, ]
  )
  
  # make the plot
  matplot(
    t(x),
    type = "l", lwd = 2, lty = 1,
    col = c("black", "red", "blue"),
    xlab = "Time", ylab = "DIC",
    main = box_name
  )
  
  legend(
    "topright",
    legend = c("c1", "c2", "c3"),
    col = c("black", "red", "blue"),
    lwd = 2, lty = 1
  )
}



### Fin


