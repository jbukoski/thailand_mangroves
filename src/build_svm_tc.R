# Script for outputting SVM classification of tasselled cap transformation
# of Landsat data

# Input files are preprocessed landsat images (see lsat_preproceesing.R)

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

year_train <- "2017_"
year_valid <- "2007_"
site <- "krabi_"
split <- 1     # specify percentage to be used for training
seed <- round(runif(1, 1, 100000))

#----------------------
# Load in data
# 1. Tasseled cap transformed LSAT data (for training "year")
# 2. Tasseled cap transformed LSAT validation data (for validation "year")
# 3. Training polygons
# 4. Validation points

lsat_tc <- brick(paste0(in_dir, year_train, site, "lsat_tc.tif"))
names(lsat_tc) <- c("brightness", "greenness", "wetness", "srtm")

lsat_tc_valid <- brick(paste0(in_dir, year_valid, site, "lsat_tc.tif"))
names(lsat_tc_valid) <- c("brightness", "greenness", "wetness", "srtm")

polys <- read_sf(paste0(in_dir, year_train, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat)) %>%
  st_as_sf()

validation_pts <- read_sf(paste0(in_dir, site, "validation.shp")) %>%
  mutate(class = ifelse(class == "agricultur", "agriculture", 
                        ifelse(class == "aquacultur", "aquaculture", class)),
         class_fct = as.numeric(factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat))

#----------------------------------
# Split polygons into training and validation datasets

set.seed(seed)
train_idx <- sample.int(nrow(polys), nrow(polys)*split)

training <- polys[train_idx, ]
valid <- polys[-train_idx, ]

#polys_spdf <- as(valid, "Spatial")   #Need SPDF for validateMap()

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
# Build water mask
# Mask B flips land/water == 0 class to remove water speckles

mat_a <- matrix(c(0, 4, 0, 4.9, 5.1, 1), ncol = 3, byrow = TRUE)

lsat_mask_a <- reclassify(lsat_tc_pred, mat_a) %>%
  rm_speckling(speckle_size = 5, dir = 4)

mat_b <- matrix(c(-0.1, 0.1, 1, 0.9, 1.1, 0), ncol = 3, byrow = TRUE)

lsat_mask_b <- reclassify(lsat_mask_a, mat_b) %>%
  rm_speckling(speckle_size = 10, dir = 4)

#----------------------
# Apply water mask

mat_c <- matrix(c(-0.1, 0.1, NA, 0.9, 1.1, 1), ncol = 3, byrow = TRUE)
lsat_water_mask <- reclassify(lsat_mask_b, mat_c)

mat_d <- matrix(c(0.9, 1.1, 1, 
                  1.9, 2.1, 2,
                  2.9, 3.1, 3,
                  3.9, 4.1, 4,
                  4.9, 5.1, 2,
                  NA, NA, 5), ncol = 3, byrow = TRUE)


masked_pred <- mask(lsat_tc_pred, lsat_water_mask)
final_lsat_tc <- reclassify(masked_pred, mat_d)

#---------------------------
# Visualize the two

par(mfrow = c(1, 2))
plot(lsat_tc_pred, main = "w/o mask")
plot(final_lsat_tc, main = "w/ mask")

#-----------------------------------
# Validation

pred_tc_valid <- raster::predict(lsat_tc_valid, svm_lsat_tc) %>%
  mask(lsat_water_mask) %>%
  reclassify(mat_d)

plot(pred_tc_valid, main = "masked prediction on validation year")

validateMap(pred_tc_valid, validation_pts, responseCol = 'class_fct', 
            1000, mode = "classification")

plot(pred_tc_valid, main = paste0(year_valid, "validation"))
plot(validation_pts, add = T)

#-----------------------------------
# # Load in other TC images
# 
# lsat1987tc <- brick(paste0(in_dir, "1987_", site, "lsat_tc.tif"))
# names(lsat1987tc) <- c("brightness", "greenness", "wetness", "srtm")
# 
# lsat1997tc <- brick(paste0(in_dir, "1997_", site, "lsat_tc.tif"))
# names(lsat1997tc) <- c("brightness", "greenness", "wetness", "srtm")
# 
# 
# pred1987tc <- raster::predict(lsat1987tc, svm_lsat_tc)
# pred1997tc <- raster::predict(lsat1997tc, svm_lsat_tc)
# 
# par(mfrow = c(2,2))
# plot(pred1987tc, main = "1987")
# plot(pred1997tc, main = "1997")
# plot(pred_valid_tc, main = "2007")
# plot(final_lsat_tc, main = "2017")
# 
# dev.off()
# 
# 
# lsat_stack <- stack(lsat1987tc, lsat1997tc, lsat_tc_valid, lsat_tc)
# 
# par(mfrow = c(2,2))
# 
# for(i in seq(2, 16, 4)) {
#   hist(lsat_stack[[i]])
# }
# 
# #-----------------------
# # Write to file
# 
# writeRaster(final_lsat_tc, paste0(out_dir, "prcssd_", "2017_", site, "svm_tc.tif"),
#             format = "GTiff", overwrite = TRUE)
