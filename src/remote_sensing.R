library(getlandsat)
library(RStoolbox)
library(raster)
library(sf)
library(rpart)
library(ggmap)
library(e1071)
library(tidyverse)
library(tidyr)

rast <- brick("~/Desktop/thailand/2015_landsat_image.tif")

srtm <- raster("~/Desktop/thailand/srtm_to_export.tif") %>%
  projectRaster(crs = crs(rast)) %>%
  resample(rast, method = "bilinear")

rast[[8]] <- srtm

polys <- read_sf("~/Desktop/thailand/training2015.shp")

#-----------------------

new_rast <- raster::extract(rast, polys)
names(new_rast) <- polys$class

training <- plyr::ldply(new_rast, rbind)
colnames(training) <- c("class", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8")
training$class <- as.factor(training$class)

classes <- c("agriculture", "aquaculture", "mangrove", "urban", "water")

index <- 1:nrow(training)
testindex <- sample(index, trunc(length(index)/3))
testset <- training[testindex,]
trainset <- training[-testindex,]

svm_model <- svm(class ~ ., data = trainset, cost = 100, gamma = 1)
svm_pred <- predict(svm_model, testset[, 2:9])

tab <- table(pred = svm_pred, true = testset[, 1])
colnames(tab) <- classes
rownames(tab) <- classes

tab

#--------------------------------------
# Run on raster

names(rast) <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8")

pred_rast <- raster::predict(rast, svm_model)
pred_rast <- round(pred_rast, 0)

plot(pred_rast)

writeRaster(pred_rast, filename = "~/Desktop/output_2.tif", format = "GTiff", overwrite = TRUE)
