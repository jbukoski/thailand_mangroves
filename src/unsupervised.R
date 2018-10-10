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
site <- "nakorn_"
split <- 1     # specify percentage to be used for training
seed <- round(runif(1, 1, 100000))

#----------------------
# Load in data
# 1. Tasseled cap transformed LSAT data (for training "year")
# 2. Tasseled cap transformed LSAT validation data (for validation "year")
# 3. Training polygons
# 4. Validation points

band_names <- c("blue", "green", "red", "nir", "swir1", "swir2", "srtm",
                "brightness", "greenness", "wetness", "ndvi", "ndwi",
                "avg3", "avg5", "avg9", "var3", "var5", "var9")

lsat <- brick(paste0(in_dir, year_train, site, "lsat.tif"))
names(lsat) <- band_names

# PCA Transformation

set.seed(1)

# red_swir <- lsat[[3]] / lsat[[5]]
# swir_nir <- lsat[[5]] / lsat[[4]]
# 
# lsat <- stack(lsat[[3:5]], red_swir, swir_nir) 

sr <- sampleRandom(lsat[[7:11]], 50000)
plot(sr[, c(3,4)])

pca <- prcomp(sr, scale = TRUE)
plot(pca)
screeplot(pca)

pci <- predict(lsat, pca, index = 1:length(pca))

#-----------------------------#
# Unsupervised classification #
#-----------------------------#
# Using kmeans clustering

v <- getValues(pci)
i <- which(!is.na(v))
kmeans_rast <- kmeans(v[i], 20, iter.max = 1000, nstart = 25, algorithm = "MacQueen")

kmeans_raster <- raster(lsat)
kmeans_raster[i] <- kmeans_rast$cluster
plot(kmeans_raster)

writeRaster(kmeans_raster, "/home/jbukoski/Desktop/test.tif", format = "GTiff", overwrite = TRUE)

threes <- matrix(1, nrow = 5, ncol = 5)

rast <- focal(kmeans_raster, threes, fun = modal)

# Reclassification

# kmeans_rc_mat <- matrix(c(0.9, 1.1, 4, 
#                           1.9, 2.1, 2,
#                           2.9, 3.1, 3,
#                           3.9, 4.1, 1,
#                           4.9, 5.1, 5), ncol = 3, byrow = TRUE)
# 
# kmeans_raster <- reclassify(kmeans_raster, kmeans_rc_mat)

# Unsupervised classification using randomForest classification w/ kmeans

vx <- v[sample(nrow(v), 500), ]
rf <- randomForest(vx)
rf_prox <- randomForest(vx, ntree = 1000, proximity = TRUE)$proximity

E_rf <- kmeans(rf_prox, 5, iter.max = 100, nstart = 10)
rf <- randomForest(vx, as.factor(E_rf$cluster), ntree = 500)
rf_raster <- predict(lsat_valid, rf)
plot(rf_raster)