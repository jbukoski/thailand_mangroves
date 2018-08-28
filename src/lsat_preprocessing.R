# Preprocessing of Landsat data
# Data transformations to be examined include:
# 1. Raw Landsat bands
# 2. Tasselled cap transformation

library(raster)
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
sensor <- "TM" # Enter either TM or OLI depending on input data sensor

#----------------------
# Load in data:
# 1. LSAT data for site and year of interest. Output from GEE script.
# 2. SRTM data for the site of interest.

lsat <- brick(paste0(raw_dir, "raw_", year, site, "lsat.tif"))

srtm <- raster(paste0(raw_dir, site, "srtm.tif")) %>%
  projectRaster(crs = crs(lsat)) %>%
  resample(lsat, method = "bilinear")

# Processing LSAT (appending SRTM band)

lsat[[7]] <- srtm

names(lsat) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

#----------------------
# Tasselled cap processing (B1: Brightness, B2: Greenness, B3: Wetness)
# Note - coefficients are sensor specific, be careful of which LSAT mission

if(sensor == "TM") {
  lsat_tc <- tasseledCap(lsat[[c(1:5, 7)]], "Landsat5TM")
  lsat_tc[[4]] <- srtm
} else if(sensor == "OLI") {
  lsat_tc <- tasseledCap(lsat[[c(1:5, 7)]], "Landsat8OLI")
  lsat_tc[[4]] <- srtm
}

names(lsat_tc) <- c("brightness", "greenness", "wetness", "srtm")

#----------------------
# Write out files to in_dir for feeding into build_classification scripts

writeRaster(lsat, paste0(in_dir, year, site, "lsat.tif"),
            format = "GTiff", overwrite = TRUE)

writeRaster(lsat_tc, paste0(in_dir, year, site, "lsat_tc.tif"), 
            format = "GTiff", overwrite = TRUE)
