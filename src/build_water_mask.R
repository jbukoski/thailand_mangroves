## Produces water masks for the two sites
## The script must be run manually, as some years require
## an interative production of a water mask (i.e., the
## classification of water-masked imagery, followed by use of
## additional water bodies to update water-mask).

library(e1071)
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

site <- "krabi_"
year <- "2017_"

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

cutoff <- -0.25

rc_mat <- matrix(c(-1, cutoff, NA,
                   cutoff, 1, 1), ncol = 3, byrow = TRUE)

water_mask <- reclassify(ndwi, rc_mat)

water_clump <- clump(water_mask, directions = 8, gaps = TRUE)

clump_freq <- data.frame(freq(water_clump)) %>%
  arrange(-count)

idx <- clump_freq[2, ]

rc_mat_2 <- matrix(c(0, idx$value - 0.10, NA,
                     idx$value - 0.10, idx$value + 0.10, 1,
                     idx$value + 0.10, 50000, NA), ncol = 3, byrow = TRUE)

water_mask <- reclassify(water_clump, rc_mat_2)

rc_mat_3 <- matrix(c(NA, NA, 1,
                     0.9, 1.1, NA), ncol = 3, byrow = TRUE)

water_mask <- reclassify(water_mask, rc_mat_3)

plot(water_mask)

#---------------------------------------------------------------
# Additional steps for iterative masking 
# (in the case that the original water mask is insufficient)
# Needed this section for nakorn_1987 and nakorn_1997.

# water <- lsat_final == 5
# 
# water <-  projectRaster(water, crs = projection(water_mask), method = "ngb") # SVM classified raster to update mask with
# water <- resample(water, water_mask, method = "ngb")
# 
# water_clump <- clump(water, directions = 8, gaps = TRUE)
# 
# clump_df <- data.frame(freq(water_clump)) %>%
#   arrange(-count)
# 
# idx <- clump_df[clump_df$count > 1200, ]
# 
# water_add <- water_clump %in% idx$value[2:length(idx$value)]
# 
# mat_rc_add <- matrix(c(0.9, 1.1, 0,
#                        NA, NA, 1), ncol = 3, byrow = TRUE)
# water_mask <- reclassify(water_mask, mat_rc_add)
# 
# water_mask_new <- water_add + water_mask
# 
# final_mask_clump <- clump((water_mask_new != 1), directions = 8, gaps = TRUE)
# water_mask <- focal((final_mask_clump == 1), w = threes, fun = modal)
# 
# # water_mask <- focal(water_mask, w = threes, fun = modal) # In case of extra filtering
# 
# plot(water_mask)

#-----------------------------------
# Write map to file

writeRaster(water_mask, paste0(out_dir, site, year, "water_mask.tif"), format = "GTiff", overwrite = TRUE)
