#######################################
# Shoo the Flu evaluation bias analysis 

# configure data directories
# source base functions
# load libraries
#######################################
# devtools::install_github("wch/webshot")

library(dplyr)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
library(tidyr)
library(ggplot2)
library(assertthat)
library(GGally)
library(gridExtra)

# load base functions
source(paste0(here::here(), "/bias-analysis/0-base-functions.R"))

# define file paths 
data_path <- "/Volumes/Data/Temp/ceip-flu-data-age-sex-race-zip.RDS"

# note: the repo must contain folders called "results" and "figures", 
# which will be populated during analyses. 