# Script for outputting SVM classification of tasselled cap transformation
# of Landsat data

# Input files are preprocessed landsat images (see lsat_preproceesing.R)

library(cluster)
library(doParallel)
library(e1071)
library(foreach)
library(randomForest)
library(raster)
library(rgdal)
library(rgeos)
library(RStoolbox)
library(sf)
library(sp)
library(tidyverse)
library(ggthemes)

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

lsat_valid <- brick(paste0(in_dir, year_valid, site, "lsat.tif"))
names(lsat_valid) <- band_names

polys <- read_sf(paste0(in_dir, year_train, site, "training.shp")) %>% 
  mutate(class_fct = as.numeric(as.factor(class))) %>%
  as("Spatial") %>%
  spTransform(crs(lsat)) %>%
  st_as_sf()

validation_pts <- read_sf(paste0(in_dir, "pts_test.shp")) %>%
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

#------------------------------#
# Classification for LSAT data #
#------------------------------#

rast <- raster::extract(lsat, training)
names(rast) <- training$class

train_df <- plyr::ldply(rast, rbind)
colnames(train_df) <- c("class", band_names)
train_df$class <- as.factor(train_df$class)

svm_lsat <- svm(class ~ blue + green + red + nir + 
                  ndvi + ndwi + brightness + greenness + avg3, 
                data = train_df, cost = 100, gamma = 1)

#----------------------------------
# TC: Run prediction on training year and validation year & visualize

useCores <- detectCores() - 2
cl <- makeCluster(useCores)
registerDoParallel(cl)

rast_stack <- list(lsat, lsat_valid)

par_output <- foreach(i=1:length(rast_stack)) %dopar% {
  library(e1071)
  library(raster)
  r <- rast_stack[[i]]
  raster::predict(r, svm_lsat)
}

lsat_pred <- par_output[[1]]
valid_pred <- par_output[[2]]

stopCluster(cl)

#----------------------------------
# Visualize the classifications

par(mfrow = c(1,2))

plot(lsat_pred, main = "2017 predicted")
plot(valid_pred, main = "2007 predicted")

#-----------------------------------
# Accuracy assessment of prediction on the validation year

acc <- validateMap(valid_pred, validation_pts, responseCol = 'class_fct',
                          1000, mode = "classification")

acc_mat <- matrix(acc$performance$table, nrow = 5, ncol = 5)

usrs_acc <- diag(acc_mat) / rowSums(acc_mat)
prds_acc <- diag(acc_mat) / colSums(acc_mat)

plot(valid_pred, main = paste0(year_valid, "validation"))
plot(validation_pts, add = T)

#------------------------------------
# Run 3x3 focal smoothing using modal function

threes <- matrix(1, nrow = 3, ncol = 3)
fives <- matrix(1, nrow = 5, ncol = 5)
valid_smooth <- focal(valid_pred, w = threes, fun = modal)

plot(valid_smooth)

acc_smth <- validateMap(valid_smooth, validation_pts, responseCol = 'class_fct',
                        1000, mode = "classification")

acc_smth_mat <- matrix(acc_smth$performance$table, nrow = 5, ncol = 5)
usrs_smth_acc <- diag(acc_smth_mat) / rowSums(acc_smth_mat)
prds_smth_acc <- diag(acc_smth_mat) / colSums(acc_smth_mat)

#------------------------------------
# Assess validation points accuracy

valid_test <- raster::extract(valid_pred, validation_pts)

valid_test <- validation_pts %>%
  st_as_sf() %>%
  mutate(pred_fct = valid_test)

plot(valid_pred)

test <- valid_test %>%
  filter(!is.na(pred_fct)) %>%
  filter(class_fct != pred_fct) %>%
  #filter(pred_fct == 4) %>%
  as("Spatial") %>%
  plot(add = TRUE, col = .$pred_fct)

#----------------------------------
# Load all layers and run predictions on them

lsat1987 <- brick(paste0(in_dir, "1987_", site, "lsat.tif")) %>%
  setNames(band_names)
lsat1997 <- brick(paste0(in_dir, "1997_", site, "lsat.tif")) %>%
  setNames(band_names)
lsat2007 <- brick(paste0(in_dir, "2007_", site, "lsat.tif")) %>%
  setNames(band_names)
lsat2017 <- brick(paste0(in_dir, "2017_", site, "lsat.tif")) %>%
  setNames(band_names)


pred_1987 <- raster::predict(lsat1987, svm_lsat) %>%
  focal(w = threes, fun = modal)
pred_1997 <- raster::predict(lsat1997, svm_lsat) %>%
  focal(w = threes, fun = modal)
pred_2007 <- raster::predict(lsat2007, svm_lsat) %>%
  focal(w = threes, fun = modal)
pred_2017 <- raster::predict(lsat2017, svm_lsat) %>%
  focal(w = threes, fun = modal)


useCores <- detectCores() - 1
cl <- makeCluster(useCores)
registerDoParallel(cl)

rast_stack <- list(lsat1987, lsat1997, lsat2007, lsat2017)

par_output <- foreach(i=1:length(rast_stack)) %dopar% {
  
  library(e1071)
  library(raster)
  library(tidyverse)
  
  r <- rast_stack[[i]]
  raster::predict(r, svm_lsat) %>%
    focal(w = threes, fun = modal)
  
}

stopCluster(cl)

par(mfrow = c(1,1))
plot(par_output[[1]], main = "1987")
plot(par_output[[2]], main = "1997")
plot(par_output[[3]], main = "2007")
plot(par_output[[4]], main = "2017")

# Write to file

writeRaster(pred_1987, paste0(out_dir, site, "1987_", "svm.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(pred_1997, paste0(out_dir, site, "1997_", "svm.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(pred_2007, paste0(out_dir, site, "2007_", "svm.tif"), format = "GTiff", overwrite = TRUE)
writeRaster(pred_2017, paste0(out_dir, site, "2017_", "svm.tif"), format = "GTiff", overwrite = TRUE)


#--------------------------------------
# Table summary values

data <- rbind(
  c(tapply(pred_1987, pred_1987[], FUN = sum) * 0.09, 1987),
  c(tapply(pred_1997, pred_1997[], FUN = sum) * 0.09, 1997),
  c(tapply(pred_2007, pred_2007[], FUN = sum) * 0.09, 2007),
  c(tapply(pred_2017, pred_2017[], FUN = sum) * 0.09, 2017)
) %>% as.data.frame()

colnames(data) <- c("agri", "aqua", "mang", "urban", "water", "year")

test_data <- melt(data, id="year")

ggplot(data=test_data,
       aes(x=year, y=value, colour=variable)) +
  geom_line()


