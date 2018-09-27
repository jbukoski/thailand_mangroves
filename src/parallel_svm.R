# Parallel code for supervised classification of Landsat imagery

library(doParallel)
library(foreach)
library(tidyverse)
library(magrittr)

#----------------

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "nakorn_"

#-----------------------------------
# CODE TESTING 

useCores <- detectCores() - 1
cl <- makeCluster(useCores)
registerDoParallel(cl)

band_names <- c("blue", "green", "red", "nir", "swir1", "swir2", "srtm",
                "brightness", "greenness", "wetness", "ndvi", "ndwi",
                "avg3", "avg5", "avg9", "var3", "var5", "var9")

#years <- c("1987_")
years <- c("1987_", "1997_", "2007_", "2017_")

areas <- foreach(i=1:length(years)) %dopar% {
  
  library(e1071)
  library(raster)
  library(sf)
  library(tidyverse)
  
  lsat <- brick(paste0(in_dir, years[i], site, "lsat.tif"))
  names(lsat) <- band_names
  
  polys <- read_sf(paste0(in_dir, years[i], site, "training.shp")) %>% 
    mutate(class_fct = as.numeric(as.factor(class))) %>%
    as("Spatial") %>%
    spTransform(crs(lsat)) %>%
    st_as_sf()
  
  #----------------------------------
  # Split polygons into training and validation datasets
  
  seed <- round(runif(1, 1, 100000))
  set.seed(seed)
  split <- 1
  
  train_idx <- sample.int(nrow(polys), nrow(polys)*split)
  
  training <- polys[train_idx, ]
  valid <- polys[-train_idx, ]
  
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
                  data = train_df, cost = 1000, gamma = 1)
  
  threes <- matrix(1, nrow = 3, ncol = 3)
  
  asia_eqArea <- crs("+proj=aea +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 
                     +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  lsat_pred <- raster::predict(lsat, svm_lsat) %>%
    focal(w = threes, fun = modal) %>%
    projectRaster(crs = asia_eqArea, method = "ngb")
  
  write_raster(lsat_pred, paste0(out_dir, site, year[i], "classified.tif"), 
               format = "GTiff", overwrite = TRUE)
  
  aggregate(getValues(area(lsat_pred, weights = FALSE)),
            by = list(getValues(lsat_pred)), sum)
  
}

summary <- areas %>%
  reduce(left_join, by = "Group.1") %>%
  set_colnames(c("Class", gsub("_", "", years)))

               