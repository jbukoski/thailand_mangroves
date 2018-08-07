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

out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"
in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"

#----------------------
# Adjust the below values to output the correct datasets

year <- "2015_"
site <- "krabi_"
split <- 0.6     # specify percentage to be used for training

#----------------------
# Load in data

rast <- brick(paste0(in_dir, year, site, "lsat.tif"))

srtm <- raster(paste0(in_dir, site, "srtm.tif")) %>%
  projectRaster(crs = crs(rast)) %>%
  resample(rast, method = "bilinear")

names(srtm) <- "srtm"

rast[[8]] <- srtm

polys <- read_sf(paste0(in_dir, year, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(rast)) %>%
  st_as_sf()

water_mask <- raster(paste0(out_dir, site, "water_mask.tif")) %>%
  round()

water_mask[water_mask == 1] <- NA

#----------------------------------
# Split polygons into training and validation datasets

set.seed(87654321)
train_idx <- sample.int(nrow(polys), nrow(polys)*0.6)

training <- polys[train_idx, ]
valid <- polys[-train_idx, ]

#-----------------------
# Tasseled cap transformation (B1: Brightness, B2: Greenness, B3: Wetness)

lsat_tc <- tasseledCap(rast[[c(1:5, 7)]], "Landsat5TM")
lsat_tc[[4]] <- srtm

#writeRaster(lsat_tc, "~/Desktop/thailand/lsat_tc.tif", format = "GTiff", overwrite = TRUE)

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

#-----------------------
# Build water mask

water_mat <- matrix(c(0, 4, 0, 4.9, 5.1, 1), ncol = 3, byrow = TRUE)
rc_a <- reclassify(pred_rast, water_mat)

rast_a <- rm_speckling(rc_a, speckle_size = 5, dir = 4)

plot(rast_a)

# reclassify water to 0 & rerun rm_speckling() to remove water speckles

mat <- matrix(c(-0.1, 0.1, 1, 0.9, 1.1, 0), ncol = 3, byrow = TRUE)
rc_b <- reclassify(rast_a, mat)

rast_b <- rm_speckling(rc_b, speckle_size = 10, dir = 4)

plot(rast_b)

#----------------------
# Apply water mask

mat <- matrix(c(-0.1, 0.1, NA, 0.9, 1.1, 1), ncol = 3, byrow = TRUE)
water_mask <- reclassify(rast_b, mat)

masked_pred <- mask(pred_rast, water_mask)

mat <- matrix(c(0.9, 1.1, 1, 
                1.9, 2.1, 2,
                2.9, 3.1, 3,
                3.9, 4.1, 4,
                4.9, 5.1, 2,
                NA, NA, 5), ncol = 3, byrow = TRUE)

rast <- reclassify(masked_pred, mat)

#---------------------------
# Visualize the two

par(mfrow = c(1, 2))
plot(pred_rast, main = "w/o mask")
plot(rast, main = "w/ mask")

#-----------------------------------
# Validation

polys_spdf <- as(valid, "Spatial")   #Need SPDF for validateMap()

validateMap(pred_rast, polys_spdf, responseCol = 'class_fct', 1000, mode = "classification")
validateMap(rast, polys_spdf, responseCol = 'class_fct', 1000, mode = "classification")

#-----------------------
# Write to file

#writeRaster(rast, paste0(out_dir, year, site, "svm_tc.tif"), format = "GTiff", overwrite = TRUE)

#-----------------------------------
# Unsupervised classification

# set.seed(040289)
# kmeans_model <- unsuperClass(lsat_tc, nSamples = 10000, nClasses = 5, nStarts = 25, nIter = 100)
# 
# plot(kmeans_model$map)

#writeRaster(kmeans_model$map, "~/Desktop/thailand/2015_kmeans.tif", format = "GTiff", overwrite = T)
