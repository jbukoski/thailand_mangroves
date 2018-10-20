# Parallel code for supervised classification of Landsat imagery

library(doParallel)
library(foreach)
library(tidyverse)
library(magrittr)

#----------------

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "krabi_"

#-----------------------------------
# CODE TESTING 

useCores <- detectCores() - 1
cl <- makeCluster(useCores)
registerDoParallel(cl)

band_names <- c("blue", "green", "red", "nir", "swir1", "swir2", "srtm",
                "brightness", "greenness", "wetness", "ndvi", "ndwi",
                "avg3", "avg5", "avg9", "var3", "var5", "var9")

years <- c("1987_", "1997_", "2007_", "2017_")

# Clear sink

writeLines(c(""), "/home/jbukoski/Desktop/log.txt")

# Run parallel loop

areas <- foreach(i=1:length(years)) %dopar% {
  
  sink("/home/jbukoski/Desktop/log.txt", append=TRUE)
  
  library(e1071)
  library(raster)
  library(RStoolbox)
  library(sf)
  library(spatialEco)
  library(tidyverse)
  
  h2oMask <- raster(paste0(out_dir, site, years[i], "water_mask.tif"))
  
  lsat <- brick(paste0(in_dir, years[i], site, "lsat.tif")) %>%
    mask(h2oMask)
  names(lsat) <- band_names
  
  polys <- read_sf(paste0(in_dir, years[i], site, "training.shp")) %>% 
    mutate(class_fct = as.numeric(as.factor(class))) %>%
    as("Spatial") %>%
    spTransform(crs(lsat)) %>%
    st_as_sf("SpatialPolygons")
  
  pts <- spsample(as(polys, "Spatial"), n = 5000, "random") %>%
    point.in.poly(polys) %>%
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
  
  print(paste("extracting data from raster...", Sys.time(), years[i]))
  
  rast <- raster::extract(lsat, training)
  names(rast) <- training$class
  
  train_df <- plyr::ldply(rast, rbind)
  colnames(train_df) <- c("class", band_names)
  train_df$class <- as.factor(train_df$class)
  
  # rast_vals_df <- plyr::ldply(rast, rbind)
  # colnames(rast_vals_df) <- c("class", band_names)
  # rast_vals_df$class <- as.factor(rast_vals_df$class)
  # 
  # idx <- sample.int(nrow(rast_vals_df), nrow(rast_vals_df) * 0.6)
  # 
  # train_df <- rast_vals_df[idx, ]
  # valid_df <- rast_vals_df[-idx, ]
  
  print(paste("training model...", Sys.time(), years[i]))
  
  svm_lsat <- svm(class ~ blue + green + red + nir + swir1 + swir2 +
                    brightness + greenness + wetness + ndvi + ndwi, 
                  data = train_df, cost = 1000, gamma = 1)
  
  threes <- matrix(1, nrow = 3, ncol = 3)
  
  asia_eqArea <- crs("+proj=aea +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 
                     +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  print(paste("classifying imagery...", Sys.time(), years[i]))
  
  lsat_pred <- raster::predict(lsat, svm_lsat) %>%
    focal(w = threes, fun = modal) %>%
    projectRaster(crs = asia_eqArea, method = "ngb")
  
  rc_mat <- matrix(c(0.9, 1.1, 1,
                     1.9, 2.1, 2,
                     2.9, 3.1, 3,
                     3.9, 4.1, 1,
                     4.9, 5.1, 2), ncol = 3, byrow = TRUE)
  
  lsat_final <- reclassify(lsat_pred, rc_mat)
  
  print(paste("writing raster...", Sys.time(), years[i]))
  
  writeRaster(lsat_final, paste0(out_dir, site, years[i], "svm.tif"), 
              format = "GTiff", overwrite = TRUE)
  
  print(paste("aggregating values...", Sys.time(), years[i]))
  
  aggregate(getValues(area(lsat_final, weights = FALSE)),
            by = list(getValues(lsat_final)), sum)
  
}

#-----------------------------
# Aggregate statistics

summary <- areas %>%
  reduce(left_join, by = "Group.1") %>%
  set_colnames(c("Class", gsub("_", "", years))) %>%
  t %>%
  as.data.frame %>%
  mutate(year = rownames(.)) %>%
  filter(year != "Class") %>%
  gather(key = "class", value = "value", -year)
               

ggplot(data = summary, 
       aes(x = year, y = value, colour = class, group = class)) +
  geom_line()
