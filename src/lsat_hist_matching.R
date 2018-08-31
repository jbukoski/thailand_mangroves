# Histogram matching of raw landsat data
# Matches all LSAT data to the 2017 imagery

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

site <- "krabi_"
sensor <- "TM" # Enter either TM or OLI depending on input data sensor

#----------------------
# Load in data:
# 1. LSAT data for site and year of interest. Output from GEE script.
# 2. SRTM data for the site of interest.

lsat2017 <- brick(paste0(raw_dir, "raw_", "2017_", site, "lsat.tif"))
lsatCRS <- crs(lsat2017)

lsat1987 <- brick(paste0(raw_dir, "raw_", "1987_", site, "lsat.tif")) %>%
  projectRaster(crs = lsatCRS)
lsat1997 <- brick(paste0(raw_dir, "raw_", "1997_", site, "lsat.tif"))  %>%
  projectRaster(crs = lsatCRS)
lsat2007 <- brick(paste0(raw_dir, "raw_", "2007_", site, "lsat.tif")) %>%
  projectRaster(crs = lsatCRS)

#------------------------
# Relative correction using histogram matching

lsat1987match <- histMatch(lsat1987, lsat2017)
lsat1997match <- histMatch(lsat1997, lsat2017)
lsat2007match <- histMatch(lsat2007, lsat2017)

#--------------------------
# Visualize

band_num <- 4

par(mfrow = c(3,3))
plot(lsat1987[[band_num]], main = '1987')
plot(lsat1987match[[band_num]], main = '1987 match')
plot(lsat2017[[band_num]], main = '2017')

plot(lsat1997[[band_num]], main = '1997')
plot(lsat1997match[[band_num]], main = '1997 match')
plot(lsat2017[[band_num]], main = '2017')

plot(lsat2007[[band_num]], main = '2007')
plot(lsat2007match[[band_num]], main = '2007 match')
plot(lsat2017[[band_num]], main = '2017')

#------------------------
# Write out histogram matched rasters

writeRaster(lsat1987match, paste0(in_dir, "matched_", "1987_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(lsat1997match, paste0(in_dir, "matched_", "1997_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(lsat2007match, paste0(in_dir, "matched_", "2007_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(lsat2017, paste0(in_dir, "matched_", "2017_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)

