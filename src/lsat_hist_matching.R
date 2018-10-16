# Histogram matching of raw landsat data
# Matches all LSAT data to the 2017 imagery

library(raster)
library(rgdal)
library(RStoolbox)
library(sf)
library(sp)

#----------------------
# Specify in/out directories

# raw_dir <- "/home/jbukoski/research/data/thailand_stocks/raw/"
# in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
# out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

#----------------------
# Adjust the below values to output the correct datasets

# site <- "nakorn_"
# roi <- "ppm_" # kre_ or ppm_

#---------------------------
# Load 2017 LSAT as a template

lsat2017 <- brick(paste0(raw_dir, "raw_", year, site, "lsat.tif")) 
lsatCRS <- crs(lsat2017)

#---------------------------
# Load site mask

mask_vec <- readOGR(paste0(in_dir, roi, "boundary.shp")) %>%
  spTransform(lsatCRS)

empty_rast <- raster(nrow = lsat2017@nrows, ncols = lsat2017@ncols, 
                     ext = lsat2017@extent, crs = lsat2017@crs)

mask_rast <- rasterize(mask_vec, empty_rast, field = "id")

#----------------------
# Load in data:
# 1. LSAT data for site and year of interest. Output from GEE script.
# 2. SRTM data for the site of interest.

lsat2017 <- mask(lsat2017, mask_rast)

lsat1987 <- stack(paste0(raw_dir, "raw_", "1987_", site, "lsat.tif")) %>%
  projectRaster(crs = lsatCRS) %>%
  resample(lsat2017, method = "bilinear") %>%
  mask(mask_rast)

lsat1997 <- brick(paste0(raw_dir, "raw_", "1997_", site, "lsat.tif")) %>%
  projectRaster(crs = lsatCRS) %>%
  resample(lsat2017, method = "bilinear") %>%
  mask(mask_rast)

lsat2007 <- brick(paste0(raw_dir, "raw_", "2007_", site, "lsat.tif")) %>%
  projectRaster(crs = lsatCRS) %>%
  resample(lsat2017, method = "bilinear") %>%
  mask(mask_rast)

srtm <- raster(paste0(raw_dir, site, "srtm.tif")) %>%
  projectRaster(crs = lsatCRS) %>%
  resample(lsat2017, method = "bilinear") %>%
  mask(mask_rast)

#------------------------
# Relative correction using histogram matching

lsat1987match <- histMatch(lsat1987, lsat2017)
lsat1997match <- histMatch(lsat1997, lsat2017)
lsat2007match <- histMatch(lsat2007, lsat2017)

#--------------------------
# Visualize

band_num <- 1

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
# Append SRTM

lsat1987match[[7]] <- srtm
lsat1997match[[7]] <- srtm
lsat2007match[[7]] <- srtm
lsat2017[[7]] <- srtm

#------------------------
# Write out histogram matched rasters

writeRaster(lsat1987match, paste0(in_dir, "matched_", "1987_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(lsat1997match, paste0(in_dir, "matched_", "1997_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(lsat2007match, paste0(in_dir, "matched_", "2007_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(lsat2017, paste0(in_dir, "matched_", "2017_", site, "lsat.tif"), format = "GTiff", overwrite = TRUE)
