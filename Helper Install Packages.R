#----------------------------
#--- Install Required Packages
#============================

## Install Packages ======================================================================
# install.packages(c(
#   "shiny",
#   "shinyAce",
#   "rsconnect",
#   "MASS",
#   "lazyeval",
#   "tidyr",
#   "dplyr",
#   "purrr",
#   "broom",
#   "ggplot2",
#   "RColorBrewer",
#   "reshape2",
#   "directlabels"
# ))

## Load Packages ======================================================================
LoadPackages <- function() {
  library(shiny,shinyAce)
  library(rsconnect)
  library(MASS)
  library(lazyeval)
  library(tidyr)
  library(dplyr)
  library(purrr)
  library(broom)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
  library(directlabels)
}

suppressWarnings(suppressMessages(LoadPackages()))

## Run Scripts ======================================================================
source("Helper Standard Atmosphere.R")
source("Helper Numerical Methods.R")
source("Helper Dataframes.R")
source("Helper Aerodynamic Functions.R")
source("Helper Initial Values.R")
source("Helper Main Calculation Functions.R")
# source("Helper Calculation Functions.R")

# theme_set(theme_linedraw())
theme_set(theme_bw())
options(scipen = 10)