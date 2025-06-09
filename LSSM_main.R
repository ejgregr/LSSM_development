#################################################################################
# Script:  LSSM_main.R
# Created: June 10 2024. EJG
# 
# Control script for the Local Seaweed Services Model. 
# Will eventually draw extensively on the regional model developed by Cam Bullen. 
#
## Updates: 
# 2025/01/06: Happy New Year. Revisiting after a few month hiatus.
# 2025/01/24: After a few weeks of development, have made some progress on 
# the parameter-based model. Reviewed Pontier, and Weigel and Pfister. And have
# conceptual design of what to do with water chemistry and dilution. 
# --> Deferring the life-cycle model until we sort out daily chemistry.
# 2025/05/27: Recent progress on using carb() and can show how a kelp plant can 
#   change water chemistry. Today, re-factored code to separate steps into 
#   different scripts. Runs clean. 
#
# TO DO: 
# - Fix day length so can remove kludges
# - Fix the q_mult fudge factor in PrepLightSim() used to match
#   calc'd DLI to reported by Pontier
# - Daily change in water parcel chemistry
#
#################################################################################

rm(list=ls(all=T)) # Erase environment.
set.seed(42)       # Setting seed for reproducibility

options(warn = -1) # CAREFUL. But the group messages are really annoying. :\

# Load packages and functions ... 
setwd( "c:/Data/Git/LSSM_development")
source( "LSSM_configuration.R" )


#---- Part 1: Load environmental drivers ----
source( "LSSM_drivers.R" )

#---- Part 2 Grow a kelp plant during main growing season (MAY to SEPT) ----
# Includes:
#   1) A simple Parametric model based on literature using simple logistic growth
#   2) Example using an unparameterised version of DEB.  
# Also exploring temperature and light inhibition factors.

source( "LSSM_growth.R" )

#---- Part 3: Calculate chemistry changes during  growing season ----

source( "LSSM_chemistry.R" )


#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 

today <- format(Sys.Date(), "%Y-%m-%d")

rmarkdown::render( "LSSM_documentation.Rmd",   
                   output_format = 'pdf_document',
                   output_dir = PDF_dir,
                   output_file = paste0( "LSSM_testing_", today ))

# FIN.
