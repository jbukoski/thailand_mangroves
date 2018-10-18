## Accuracy assessment

library(raster)
library(tidyverse)

#---------------------------------

source("./src/helper_funcs.R")

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "krabi_"
year <- "2017_"

#---------------------------------
# Script for assessing the accuracy of the classified imagery

img <- brick(paste0(out_dir, site, year, "svm.tif"))
