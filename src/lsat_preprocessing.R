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

year <- "2017_"
site <- "krabi_"
sensor <- "OLI" # Enter 4TM, 5TM or OLI depending on input data sensor

#----------------------
# Load in data:
# 1. LSAT data for site and year of interest. Output from GEE script.
# 2. SRTM data for the site of interest.

lsat <- brick(paste0(in_dir, "matched_", year, site, "lsat.tif"))

names(lsat) <- c("blue", "green", "red", "nir", "swir1", "swir2", "srtm")

#----------------------
# Tasselled cap processing (B1: Brightness, B2: Greenness, B3: Wetness)
# Note - coefficients are sensor specific, be careful of which LSAT mission

if(sensor == "4TM") {
  lsat_tc <- tasseledCap(lsat[[c(1:6)]], "Landsat4TM")
  
  ndvi <- ( lsat[[4]] - lsat[[3]] ) / ( lsat[[4]] + lsat[[3]] )
  rc_mat <- matrix(c(1, 100, 1, -100, -1, -1), nrow = 2, ncol = 3, byrow = T)
  ndvi <- reclassify(ndvi, rc_mat)
  
  ndwi <- ( lsat[[2]] - lsat[[4]] ) / ( lsat[[2]] + lsat[[4]] )  # McFeeters NDWI 1996
  ndwi <- reclassify(ndwi, rc_mat)
  
  lsat_tc[[4]] <- lsat[[7]]
  lsat_tc[[5]] <- ndvi
  lsat_tc[[6]] <- ndwi

} else if(sensor == "5TM") {
  lsat_tc <- tasseledCap(lsat[[c(1:6)]], "Landsat5TM")
  
  ndvi <- ( lsat[[4]] - lsat[[3]] ) / ( lsat[[4]] + lsat[[3]] )
  rc_mat <- matrix(c(1, 100, 1, -100, -1, -1), nrow = 2, ncol = 3, byrow = T)
  ndvi <- reclassify(ndvi, rc_mat)
  
  ndwi <- ( lsat[[2]] - lsat[[4]] ) / ( lsat[[2]] + lsat[[4]] )  # McFeeters NDWI 1996
  ndwi <- reclassify(ndwi, rc_mat)
  
  lsat_tc[[4]] <- lsat[[7]]
  lsat_tc[[5]] <- ndvi
  lsat_tc[[6]] <- ndwi
  
} else if(sensor == "OLI") {
  lsat_tc <- tasseledCap(lsat[[c(1:6)]], "Landsat8OLI")
  
  ndvi <- ( lsat[[5]] - lsat[[4]] ) / ( lsat[[5]] + lsat[[4]] )
  rc_mat <- matrix(c(1, 100, 1, -100, -1, -1), nrow = 2, ncol = 3, byrow = T)
  ndvi <- reclassify(ndvi, rc_mat)
  
  ndwi <- ( lsat[[3]] - lsat[[5]] ) / ( lsat[[3]] + lsat[[5]] )  # McFeeters NDWI 1996
  ndwi <- reclassify(ndwi, rc_mat)
  
  lsat_tc[[4]] <- lsat[[7]]
  lsat_tc[[5]] <- ndvi
  lsat_tc[[6]] <- ndwi
  
}

#--------------------------
# Texture metrics

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

lsat_tc[[7]] <- ndvi_avg_3
lsat_tc[[8]] <- ndvi_avg_5
lsat_tc[[9]] <- ndvi_avg_9
lsat_tc[[10]] <- ndvi_var_3
lsat_tc[[11]] <- ndvi_var_5
lsat_tc[[12]] <- ndvi_var_9

names(lsat_tc) <- c("brightness", "greenness", "wetness", "srtm", "ndvi", "ndwi",
                    "avg3", "avg5", "avg9", "var3", "var5", "var9")

#----------------------
# Write out files to in_dir for feeding into build_classification scripts

writeRaster(lsat, paste0(in_dir, year, site, "lsat.tif"),
            format = "GTiff", overwrite = TRUE)

writeRaster(lsat_tc, paste0(in_dir, year, site, "lsat_tc.tif"), 
            format = "GTiff", overwrite = TRUE)
