### Produces a water mask for the Krabi region
## This is just a change

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

site <- "krabi_"
year <- "2015_"

#----------------------
# Load in data

rast <- brick(paste0(in_dir, year, site, "lsat.tif"))

srtm <- raster(paste0(in_dir, site, "srtm.tif")) %>%
  projectRaster(crs = crs(rast)) %>%
  resample(rast, method = "bilinear")

names(srtm) <- "srtm"

rast[[8]] <- srtm

polys <- read_sf(paste0(in_dir, site, "water_training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(rast)) %>%
  st_as_sf()

#----------------------------------
# Split polygons into training and validation datasets

set.seed(04021989)
train_idx <- sample.int(nrow(polys), nrow(polys)*.8)

training <- polys[train_idx, ]
valid <- polys[-train_idx, ]

#-----------------------
# Tasseled cap transformation (B1: Brightness, B2: Greenness, B3: Wetness)

lsat_tc <- tasseledCap(rast[[c(1:5, 7)]], "Landsat5TM")
lsat_tc[[4]] <- srtm

#--------------------------
# Extract training data and train SVM

tc_rast <- raster::extract(lsat_tc, training)
names(tc_rast) <- training$class

train_df <- plyr::ldply(tc_rast, rbind)
colnames(train_df) <- c("class", "brightness", "greenness", "wetness", "srtm")
train_df$class <- as.factor(train_df$class)

svm_model <- svm(class ~ ., data = train_df, cost = 100, gamma = 1)

#---------------------
# Run prediction and visualize

pred_rast <- raster::predict(lsat_tc, svm_model)
values(pred_rast) <- as.numeric(values(pred_rast))

plot(pred_rast)

#-----------------------------------
# Validation

polys_spdf <- as(polys, "Spatial")   #Need SPDF for validateMap()

validateMap(pred_rast, polys_spdf, responseCol = 'class_fct', 1000, mode = "classification")

#--------------------
# Use clumping function to remove speckling (likely aquaculture)

mat <- matrix(c(0, 1, 0, 1, 2, 1), ncol = 3, byrow = TRUE)
rc <- reclassify(pred_rast, mat)

rast_a <- rm_speckling(rc, speckle_size = 5, dir = 4)

plot(rast_a)

# reclassify water to 0 & rerun rm_speckling() to remove water speckles

mat <- matrix(c(-0.1, 0.1, 1, 0.9, 1.1, 0), ncol = 3, byrow = TRUE)
rc <- reclassify(rast_a, mat)

rast_b <- rm_speckling(rc, speckle_size = 10, dir = 4)

plot(rast_b)

#-----------------------------------
# Write map to file

writeRaster(water_mask, paste0(out_dir, site, "water_mask.tif"), format = "GTiff", overwrite = TRUE)

#-----------------------------------
# Code removed from "build_svm.R"


#-----------------------
# Build water mask
# Mask B flips land/water == 0 class to remove water speckles
# 
# mat_a <- matrix(c(0, 4, 0, 4.9, 5.1, 1), ncol = 3, byrow = TRUE)
# 
# lsat_mask_a <- reclassify(lsat_pred, mat_a) %>%
#   rm_speckling(speckle_size = 5, dir = 4)
# 
# mat_b <- matrix(c(-0.1, 0.1, 1, 0.9, 1.1, 0), ncol = 3, byrow = TRUE)
# 
# lsat_mask_b <- reclassify(lsat_mask_a, mat_b) %>%
#   rm_speckling(speckle_size = 10, dir = 4)
# 
# #----------------------
# # Apply water mask
# 
# mat_c <- matrix(c(-0.1, 0.1, NA, 0.9, 1.1, 1), ncol = 3, byrow = TRUE)
# lsat_water_mask <- reclassify(lsat_mask_b, mat_c)
# 
# mat_d <- matrix(c(0.9, 1.1, 1, 
#                   1.9, 2.1, 2,
#                   2.9, 3.1, 3,
#                   3.9, 4.1, 4,
#                   4.9, 5.1, 2,
#                   NA, NA, 5), ncol = 3, byrow = TRUE)
# 
# 
# masked_pred <- mask(lsat_pred, lsat_water_mask)
# final_lsat <- reclassify(masked_pred, mat_d)
# 
# #---------------------------
# # Visualize the two
# 
# par(mfrow = c(1, 2))
# plot(lsat_pred, main = "w/o mask")
# plot(final_lsat, main = "w/ mask")
