### Produces a water mask for the Krabi region

library(e1071)
library(getlandsat)
library(RStoolbox)
library(raster)
library(sf)
library(sp)
library(tidyverse)

source("/home/jbukoski/research/scripts/thailand_stocks/src/helper_funcs.R")

#----------------------
# Specify in/out directories

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

#----------------------
# Adjust the below values to output the correct water mask

site <- "nakorn_"
year <- "1987_"

#----------------------
# Load in data

band_names <- c("blue", "green", "red", "nir", "swir1", "swir2", "srtm",
                "brightness", "greenness", "wetness", "ndvi", "ndwi",
                "avg3", "avg5", "avg9", "var3", "var5", "var9")

rast <- brick(paste0(in_dir, year, site, "lsat.tif"))
names(rast) <- band_names

ndwi <- rast[["ndwi"]]

#--------------------
# Reclassify based on water as < 0

rc_mat <- matrix(c(-1, 0.0, NA,
                   0.0, 1, 1), ncol = 3, byrow = TRUE)

water_mask <- reclassify(ndwi, rc_mat)

water_clump <- clump(water_mask, directions = 8, gaps = TRUE)

clump_freq <- data.frame(freq(water_clump)) %>%
  arrange(-count)

idx <- clump_freq[2, ]

rc_mat_2 <- matrix(c(idx$value-0.1, idx$value+0.1, 1,
                     2, 50000, NA), ncol = 3, byrow = TRUE)

water_mask <- reclassify(water_clump, rc_mat_2)

plot(water_mask)

#-----------------------------------
# Write map to file

writeRaster(water_mask, paste0(out_dir, site, year, "water_mask.tif"), format = "GTiff", overwrite = TRUE)
