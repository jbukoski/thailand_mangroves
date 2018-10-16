##############
## runall.R ##
##############
#
# A driver script for running the full remote sensing analysis
#

library(raster)
library(rgdal)
library(RStoolbox)
library(sf)
library(sp)
library(tidyverse)

source("/home/jbukoski/research/scripts/thailand_stocks/src/helper_funcs.R")

# Specify input/output directories

raw_dir <- "/home/jbukoski/research/data/thailand_stocks/raw/"
in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

# Specify meta variables for processing of data

site <- "krabi_"
roi <- "kre_"
year <- "1987_"
#years <- c("1987_", "1997_", "2007_", "2017_")

#-----------------------------------
# Run Landsat histogram matching:
# Performs histogram matching on the raw Landsat data
# Runs on all of the data, so no need to repeatedly perform

# source("./src/lsat_hist_matching.R")

#-----------------------------------
# Run Landsat preprocessing:
# Performs transformation of landsat data.
# Outputs "lsat" variable, as well as writes the processed raster
# to the out_dir

source("./src/lsat_preprocessing.R")

#-----------------------------------
# Run 

# Water masks

# SVM classification

# Accuracy assessment

source("./src/w")
