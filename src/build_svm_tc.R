# Script for outputting SVM classification of tasselled cap transformation
# of Landsat data

# Input files are preprocessed landsat images (see lsat_preproceesing.R)

library(cluster)
library(e1071)
library(randomForest)
library(raster)
library(rgdal)
library(rgeos)
library(RStoolbox)
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

band_names <- c("brightness", "greenness", "wetness", "srtm", "ndvi", "ndwi",
                "avg3", "avg5", "avg9", "var3", "var5", "var9")

lsat_tc <- brick(paste0(in_dir, year_train, site, "lsat_tc.tif"))
names(lsat_tc) <- band_names

lsat_tc_valid <- brick(paste0(in_dir, year_valid, site, "lsat_tc.tif"))
names(lsat_tc_valid) <- band_names

polys <- read_sf(paste0(in_dir, year_train, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat_tc)) %>%
  st_as_sf()

validation_pts <- read_sf(paste0(in_dir, site, "reg_vld_pts.shp")) %>%
  mutate(class = ifelse(class == "agricultur", "agriculture", 
                        ifelse(class == "aquacultur", "aquaculture", class)),
         class_fct = as.numeric(factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat_tc))

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
colnames(train_df) <- c("class", band_names)
train_df$class <- as.factor(train_df$class)

svm_lsat_tc <- svm(class ~ brightness + srtm + ndvi, 
                   data = train_df, cost = 100, gamma = 1)

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

pred_tc_valid <- raster::predict(lsat_tc_valid, svm_lsat_tc)
values(lsat_tc_pred) <- as.numeric(values(pred_tc_valid))

plot(pred_tc_valid, main = "masked prediction on validation year")


( accuracy <- validateMap(lsat_tc_pred, validation_pts, responseCol = 'class_fct', 
            1000, mode = "classification") )

acc_mat <- matrix(accuracy$performance$table, nrow = 5, ncol = 5)
users <- diag(acc_mat) / rowSums(acc_mat)
producers <- diag(acc_mat) / colSums(acc_mat)

users
producers

plot(pred_tc_valid, main = paste0(year_valid, "validation"))
plot(validation_pts, add = T)

#-----------------------------#
# Unsupervised classification #
#-----------------------------#
# Using kmeans clustering

v <- getValues(lsat_tc_valid)
i <- which(!is.na(v))
kmeans_rast <- kmeans(v[i], 12, iter.max = 100, nstart = 10)

kmeans_raster <- raster(lsat_tc)
kmeans_raster[i] <- kmeans_rast$cluster
plot(kmeans_raster)

kmeans_rc_mat <- matrix(c(0.9, 1.1, 4, 
                          1.9, 2.1, 2,
                          2.9, 3.1, 3,
                          3.9, 4.1, 1,
                          4.9, 5.1, 5), ncol = 3, byrow = TRUE)

kmeans_raster <- reclassify(kmeans_raster, kmeans_rc_mat)

validateMap(kmeans_raster, validation_pts, responseCol = 'class_fct', 
            1000, mode = "classification")

plot(lsat_tc_pred)

# Unsupervised classification using randomForest classification w/ kmeans

vx <- v[sample(nrow(v), 500), ]
rf <- randomForest(vx)
rf_prox <- randomForest(vx, ntree = 1000, proximity = TRUE)$proximity

E_rf <- kmeans(rf_prox, 5, iter.max = 100, nstart = 10)
rf <- randomForest(vx, as.factor(E_rf$cluster), ntree = 500)
rf_raster <- predict(lsat_tc_valid, rf)
plot(rf_raster)

#----------------------------------
# Load all layers and run predictions on them

tc_1987 <- brick(paste0(in_dir, "1987_", site, "lsat_tc.tif"))
names(tc_1987) <- band_names

tc_1997 <- brick(paste0(in_dir, "1997_", site, "lsat_tc.tif"))
names(tc_1997) <- band_names

tc_2007 <- brick(paste0(in_dir, "2007_", site, "lsat_tc.tif"))
names(tc_2007) <- band_names

tc_2017 <- brick(paste0(in_dir, "2017_", site, "lsat_tc.tif"))
names(tc_2017) <- band_names

pred_1987 <- raster::predict(tc_1987, svm_lsat_tc)
pred_1997 <- raster::predict(tc_1997, svm_lsat_tc)
pred_2007 <- raster::predict(tc_2007, svm_lsat_tc)
pred_2017 <- raster::predict(tc_2017, svm_lsat_tc)

par(mfrow = c(2,2))
plot(pred_1987, main = '1987')
plot(pred_1997, main = '1997')
plot(pred_2007, main = '2007')
plot(pred_2017, main = '2017')
