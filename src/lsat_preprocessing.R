# Preprocessing of Landsat data
# Data transformations to be examined include:
# 1. Raw Landsat bands
# 2. Tasselled cap transformation

library(raster)
library(rgdal)
library(RStoolbox)
library(sf)
library(sp)

#----------------------
# Specify in/out directories

raw_dir <- "/home/jbukoski/research/data/thailand_stocks/raw/"
in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

#----------------------
# Adjust the below values to output the correct datasets

year <- "1987_"
site <- "nakorn_"

# Enter 4TM, 5TM or OLI depending on input data sensor

if(year == "1987_") {
  sensor <- "Landsat4TM"
} else if(year == "2017_") {
  sensor <- "Landsat8OLI"
} else {
  sensor <- "Landsat5TM"
}

#----------------------
# Load in data:
# 1. LSAT data for site and year of interest. Output from GEE script.
# 2. SRTM data for the site of interest.

lsat <- stack(paste0(in_dir, "matched_", year, site, "lsat.tif"))

names(lsat) <- c("blue", "green", "red", "nir", "swir1", "swir2", "srtm")

#----------------------
# Generate Tasseled Cap transformation, NDVI and NDWI indices
# Note - coefficients are sensor specific, be careful of which LSAT mission

lsat_tc <- tasseledCap(lsat[[c(1:6)]], sensor)

ndvi <- ( lsat[[4]] - lsat[[3]] ) / ( lsat[[4]] + lsat[[3]] )
rc_mat <- matrix(c(1, 100, 1, -100, -1, -1), nrow = 2, ncol = 3, byrow = T)
ndvi <- reclassify(ndvi, rc_mat)

ndwi <- ( lsat[[2]] - lsat[[4]] ) / ( lsat[[2]] + lsat[[4]] )  # McFeeters NDWI 1996
ndwi <- reclassify(ndwi, rc_mat)

lsat <- stack(lsat, lsat_tc, ndvi, ndwi)
names(lsat) <- c(names(lsat)[1:10], 'ndvi', 'ndwi')

#-------------------------
# Generate PCA transformation



#--------------------------
# Generate texture metrics

# Specify window sizes
threes <- matrix(1, nrow = 3, ncol = 3)
fives <- matrix(1, nrow = 5, ncol = 5)
nines <- matrix(1, nrow = 9, ncol = 9)

# Run for different metrics
ndvi_avg_3 <- focal(ndvi, w = threes, fun = mean)
ndvi_avg_5 <- focal(ndvi, w = fives, fun = mean)
ndvi_avg_9 <- focal(ndvi, w = nines, fun = mean)

ndvi_var_3 <- focal(ndvi, w = threes, fun = var)
ndvi_var_5 <- focal(ndvi, w = fives, fun = var)
ndvi_var_9 <- focal(ndvi, w = nines, fun = var)

texture_names <- c(names(lsat), "avg3", "avg5", "avg9", "var3", "var5", "var9")

lsat <- stack(lsat, ndvi_avg_3, ndvi_avg_5, ndvi_avg_9, 
              ndvi_var_3, ndvi_var_5, ndvi_var_9)
names(lsat) <- texture_names

#----------------------
# Write out files to in_dir for feeding into build_classification scripts

writeRaster(lsat, paste0(in_dir, year, site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
