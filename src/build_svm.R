# Script for outputting SVM classification of tasselled cap transformation
# of Landsat data

library(e1071)
library(RStoolbox)
library(raster)
library(sf)
library(sp)
library(tidyverse)

#----------------------
# Load in helper functions

source("/home/jbukoski/research/scripts/thailand_stocks/src/helper_funcs.R")

#----------------------
# Specify in/out directories

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

#----------------------
# Adjust the below values to output the correct datasets

year <- "2015_"
site <- "krabi_"
split <- 0.6     # specify percentage to be used for training
seed <- round(runif(1, 1, 100000))

#----------------------
# Load in data

lsat <- brick(paste0(in_dir, year, site, "lsat.tif"))
names(lsat) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "srtm")

lsat_tc <- brick(paste0(in_dir, year, site, "lsat_tc.tif"))
names(lsat_tc) <- c("brightness", "greenness", "wetness", "srtm")

polys <- read_sf(paste0(in_dir, year, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat)) %>%
  st_as_sf()

#----------------------------------
# Split polygons into training and validation datasets

set.seed(seed)
train_idx <- sample.int(nrow(polys), nrow(polys)*split)

training <- polys[train_idx, ]
valid <- polys[-train_idx, ]

polys_spdf <- as(valid, "Spatial")   #Need SPDF for validateMap()

#------------------------------#
# Classification for LSAT data #
#------------------------------#

rast <- raster::extract(lsat, training)
names(rast) <- training$class

train_df <- plyr::ldply(rast, rbind)
colnames(train_df) <- c("class", "B1", "B2", "B3", "B4", "B5", "B6", 
                        "B7", "srtm")
train_df$class <- as.factor(train_df$class)

svm_lsat <- svm(class ~ ., data = train_df, cost = 100, gamma = 1)

#---------------------
# LSAT: Run prediction and visualize

lsat_pred <- raster::predict(lsat, svm_lsat)
values(lsat_pred) <- as.numeric(values(lsat_pred))

plot(lsat_pred)

#-----------------------
# Build water mask
# Mask B flips land/water == 0 class to remove water speckles

mat_a <- matrix(c(0, 4, 0, 4.9, 5.1, 1), ncol = 3, byrow = TRUE)

lsat_mask_a <- reclassify(lsat_pred, mat_a) %>%
  rm_speckling(speckle_size = 5, dir = 4)

mat_b <- matrix(c(-0.1, 0.1, 1, 0.9, 1.1, 0), ncol = 3, byrow = TRUE)

lsat_mask_b <- reclassify(lsat_mask_a, mat_b) %>%
  rm_speckling(speckle_size = 10, dir = 4)

#----------------------
# Apply water mask

mat_c <- matrix(c(-0.1, 0.1, NA, 0.9, 1.1, 1), ncol = 3, byrow = TRUE)
lsat_water_mask <- reclassify(lsat_mask_b, mat_c)

masked_pred <- mask(lsat_pred, lsat_water_mask)

mat_d <- matrix(c(0.9, 1.1, 1, 
                1.9, 2.1, 2,
                2.9, 3.1, 3,
                3.9, 4.1, 4,
                4.9, 5.1, 2,
                NA, NA, 5), ncol = 3, byrow = TRUE)

final_lsat <- reclassify(masked_pred, mat_d)

#-----------------------------------
# Validation & visualize

validateMap(final_lsat, polys_spdf, responseCol = 'class_fct', 
            1000, mode = "classification")

plot(final_lsat, main = "w/ water mask")

#-----------------------------------
# Write out LSAT raster

# writeRaster(final_lsat, paste0(out_dir, "prcssd_", year, site, "lsat.tif"), 
#             format = "GTiff", overwrite = TRUE)

#---------------------------------#
# Classification for LSAT_TC data #
#---------------------------------#

tc_rast <- raster::extract(lsat_tc, training)
names(tc_rast) <- training$class

train_df <- plyr::ldply(tc_rast, rbind)
colnames(train_df) <- c("class", "brightness", "greenness", "wetness", "srtm")
train_df$class <- as.factor(train_df$class)

svm_lsat_tc <- svm(class ~ ., data = train_df, cost = 100, gamma = 1)

#---------------------
# TC: Run prediction and visualize

lsat_tc_pred <- raster::predict(lsat_tc, svm_lsat_tc)
values(lsat_tc_pred) <- as.numeric(values(lsat_tc_pred))

plot(lsat_tc_pred, main = "lsat_tc_pred")

#-----------------------
# Build water mask; Uses matrices from lsat section
# Mask B flips land/water == 0 class to remove water speckles

lsat_tc_mask_a <- reclassify(lsat_tc_pred, mat_a) %>%
  rm_speckling(speckle_size = 5, dir = 4)

lsat_tc_mask_b <- reclassify(lsat_tc_mask_a, mat_b) %>%
  rm_speckling(speckle_size = 10, dir = 4)

#----------------------
# Apply water mask

lsat_tc_water_mask <- reclassify(lsat_tc_mask_b, mat_c)

masked_pred <- mask(lsat_tc_pred, lsat_tc_water_mask)

mat <- matrix(c(0.9, 1.1, 1, 
                1.9, 2.1, 2,
                2.9, 3.1, 3,
                3.9, 4.1, 4,
                4.9, 5.1, 2,
                NA, NA, 5), ncol = 3, byrow = TRUE)

final_lsat_tc <- reclassify(masked_pred, mat)

#---------------------------
# Visualize the two

par(mfrow = c(1, 2))
plot(lsat_tc_pred, main = "w/o mask")
plot(final_lsat_tc, main = "w/ mask")

#-----------------------------------
# Validation

validateMap(final_lsat_tc, polys_spdf, responseCol = 'class_fct', 
            1000, mode = "classification")

#-----------------------
# Write to file

# writeRaster(final_lsat_tc, paste0(out_dir, "prcssd", year, site, "svm_tc.tif"),
#             format = "GTiff", overwrite = TRUE)
