#################################################################################
# Script:  LSSM_growth.R
# Created: January 21, 2024. EJG
# 
# Functions to support various aspects of Nereo growth and scaling. 
# Descriptions of source data and model design things in the RMD file.
# Assumptions:
#   - Currently assuming temperature-mediated logistic growth, to a maximum. 
#   - Daily growth seems reasonable. 
# 
# ## Updates: 
# 2025/xx/xx: 
#################################################################################


#---- Utility growth functions ----
# Logistic growth, with option for temperature and light factors.
logistic_growth <- function(B0, K, r, step, temp=1, lite=1) {
  B_predicted <- K / (1 + ((K - B0) / B0) * exp(-r * step * temp * lite))
  return(B_predicted)
}

# Calculate temperature effect on growth rate
# Pontier shows growth rate relatively stable below 10C, declines to about 1/2 by 14C
t_scale <- function(T) {
  scaled <- ifelse(T <= 10, 1, 0.5 / (T / T_max * 0.7))
  return( scaled )
}

# Light-dependent growth rate function
# Pontier shows growth rate peaks ~30 DLI. Simplify their decline to be 1/2 on both sides  
DLI_scale <- function(lite) {
  #(0.75*((lite - DLI_opt)^.5) / (DLI_range^2))
   1 - (((lite - DLI_opt) / DLI_range)^2) * 0.5
  
}


`#----- Grow a kelp plant during PRIMARY growing season (MAY to SEPT) -----

# Step 1 - simple logistic growth
y <- NULL
for (t in 1:round(365/2) ) {
  one <- logistic_growth( B_init, B_max, r_max, t)
  y <- c(y, one)
}
log_simp <- y

# Step 2 - Add temperature and light variability
temp_fact <- t_scale( t_daily )
DLI_fact  <- DLI_scale( light_DLI )

#Grow!
y <- NULL
for (t in 1:length(temp_fact) ) {
  one <- logistic_growth( B_init, B_max, r_max*temp_fact[t], t)
#  one <- logistic_growth( B_init, B_max, t, r_max*temp_fact)
  y <- c(y, one)
}
log_t <- y

y <- NULL
for (t in 1:length(temp_fact) ) {
  one <- logistic_growth( B_init, B_max, r_max*DLI_fact[t], t)
  #  one <- logistic_growth( B_init, B_max, t, r_max*temp_fact)
  y <- c(y, one)
}
log_DLI <- y

y <- NULL
for (t in 1:length(temp_fact) ) {
  one <- logistic_growth( B_init, B_max, r_max*DLI_fact[t]*temp_fact[t], t)
  #  one <- logistic_growth( B_init, B_max, t, r_max*temp_fact)
  y <- c(y, one)
}
log_env <- y

#----- Plot results -----
dev.off()
par(mfcol=c(4,2))

xvals <- timestamps

plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
plot( temp_fact, type = 'l' )
plot( DLI_fact, type = 'l' )
plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")


plot(log_simp*1000, type='l',xlab = 'Days', ylab = 'grams', main='A Plant')
plot(log_t*1000,    type='l',xlab = 'Days', ylab = 'grams', main='A Plant with Temp')
plot(log_DLI*1000,  type='l',xlab = 'Days', ylab = 'grams', main='A Plant with Light')
plot(log_env*1000,  type='l',xlab = 'Days', ylab = 'grams', main='A Plant with both')


#-------------
# 6 month daily growth rate, R found manually. 
r <- 0.065
x <- NULL
for (t in 1:round(365/2) ) {
  y <- logistic_growth( 0.001, 9.23, t, r)
  x <- c( x, y)
}
plot(x)
