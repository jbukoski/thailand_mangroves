# Accuracy assessment

library(raster)
library(rgdal)
library(sf)
library(tidyverse)

#---------------------------------

source("./src/helper_funcs.R")

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "krabi_"
year <- "2017_"

#---------------------------------

acc <- readOGR(dsn = paste0(out_dir), layer = paste0(site, year, "accuracy"))

table(acc$class, acc$ref_class)

