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

lsat <- brick(paste0(in_dir, year_train, site, "lsat.tif"))
names(lsat) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

lsat_valid <- brick(paste0(in_dir, year_valid, site, "lsat.tif"))
names(lsat_valid) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

lsat_tc <- brick(paste0(in_dir, year_train, site, "lsat_tc.tif"))
names(lsat_tc) <- c("brightness", "greenness", "wetness", "srtm")

lsat_tc_valid <- brick(paste0(in_dir, year_valid, site, "lsat_tc.tif"))
names(lsat_tc_valid) <- c("brightness", "greenness", "wetness", "srtm")

polys <- read_sf(paste0(in_dir, year_train, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat)) %>%
  st_as_sf()

polys_valid <- read_sf(paste0(in_dir, year_valid, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat))

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
colnames(train_df) <- c("class", "B1", "B2", "B3", "B4", "B5", 
                        "B7", "srtm")
train_df$class <- as.factor(train_df$class)

svm_lsat <- svm(class ~ B2 + B3 + B4 + srtm, 
                data = train_df, cost = 100, gamma = 1)

#---------------------
# LSAT: Run prediction and visualize

lsat_pred <- raster::predict(lsat, svm_lsat)
values(lsat_pred) <- as.numeric(values(lsat_pred))

plot(lsat_pred, main = paste0(year_train, "training"))

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

plot(final_lsat, main = paste0(year_train, "training w/ water mask"))

#-----------------------------------
# Validate with validation year
# 1. Apply model prediction to 2007 data
# 2. Use 2007 training data to validate the prediction

lsat_pred_valid <- raster::predict(lsat_valid, svm_lsat)
values(lsat_pred_valid) <- as.numeric(values(lsat_pred_valid))

validateMap(lsat_pred_valid, polys_valid, responseCol = 'class_fct', 
            1000, mode = "classification")

plot(lsat_pred_valid, main = paste0(year_valid, "validation"))

#-----------------------------------
# Predict for all years

lsat1987 <- brick(paste0(in_dir, "1987_", site, "lsat.tif"))
names(lsat1987) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

lsat1997 <- brick(paste0(in_dir, "1997_", site, "lsat.tif"))
names(lsat1997) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

lsat2007 <- brick(paste0(in_dir, "2007_", site, "lsat.tif"))
names(lsat2007) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

lsat2017 <- brick(paste0(in_dir, "2017_", site, "lsat.tif"))
names(lsat2017) <- c("B1", "B2", "B3", "B4", "B5", "B7", "srtm")

pred1987 <- raster::predict(lsat1987, svm_lsat)
values(pred1987) <- as.numeric(values(pred1987))
pred1997 <- raster::predict(lsat1997, svm_lsat)
pred2007 <- raster::predict(lsat2007, svm_lsat)
pred2017 <- raster::predict(lsat2017, svm_lsat)


par(mfrow = c(2,2))
plot(pred1987, main = "1987")
plot(pred1997, main = "1997")
plot(pred2007, main = "2007")
plot(pred2017, main = "2017")

lsat_stack <- stack(lsat1987, lsat1997, lsat2007, lsat2017)

for(i in seq(7, 28, 7)) {
  hist(lsat_stack[[i]])
}



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

pred_valid_tc <- raster::predict(lsat_tc_valid, svm_lsat_tc)

par(mfrow = c(1, 2))
plot(final_lsat_tc, main = paste0(year_train, "training"))
plot(pred_valid_tc, main = paste0(year_valid, "validation"))

# Load in other TC images

lsat1987tc <- brick(paste0(in_dir, "1987_", site, "lsat_tc.tif"))
names(lsat1987tc) <- c("brightness", "greenness", "wetness", "srtm")

lsat1997tc <- brick(paste0(in_dir, "1997_", site, "lsat_tc.tif"))
names(lsat1997tc) <- c("brightness", "greenness", "wetness", "srtm")


pred1987tc <- raster::predict(lsat1987tc, svm_lsat_tc)
pred1997tc <- raster::predict(lsat1997tc, svm_lsat_tc)

par(mfrow = c(2,2))
plot(pred1987tc, main = "1987")
plot(pred1997tc, main = "1997")
plot(pred_valid_tc, main = "2007")
plot(final_lsat_tc, main = "2017")

dev.off()


lsat_stack <- stack(lsat1987tc, lsat1997tc, lsat_tc_valid, lsat_tc)

par(mfrow = c(2,2))

for(i in seq(2, 16, 4)) {
  hist(lsat_stack[[i]])
}

#-----------------------
# Write to file

# writeRaster(final_lsat_tc, paste0(out_dir, "prcssd", year, site, "svm_tc.tif"),
#             format = "GTiff", overwrite = TRUE)
