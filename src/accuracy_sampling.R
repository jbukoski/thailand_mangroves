## Accuracy assessment

library(raster)
library(rgdal)
library(sf)
library(tidyverse)

#---------------------------------

source("./src/helper_funcs.R")

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "nakorn_"
year <- "2017_"

#--------------------------------
# Define and set seed for reproducibility

seed <- 04021989

#---------------------------------
# Script for assessing the accuracy of the classified imagery

img <- brick(paste0(out_dir, site, year, "svm.tif")) %>%
  projectRaster(crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "), method = "ngb")

vals <- getValues(img) %>%
  tibble() %>%
  rename(value = ".") %>%
  mutate(rast_id = row_number(),
        class = ifelse(is.na(value), "none", 
                  ifelse(value == 1, "agri", 
                    ifelse(value == 2, "aqua", 
                      ifelse(value == 3, "mang", NA))))) %>%
  mutate(class = as.factor(class))

set.seed(seed)
valid <- vals %>%
  group_by(class) %>%
  sample_n(200, replace = FALSE)

pts <- valid %>%
  mutate(long = xyFromCell(img, rast_id)[, 1],
         lat = xyFromCell(img, rast_id)[, 2]) %>%
  filter(class != "none") %>%
  st_as_sf(coords = c("long", "lat"))

sp_pts <- as(pts, "Spatial")

plot(img)
plot(sp_pts, col = sp_pts$class, add = T)

#----------------------------
# Write out validation pts

writeOGR(sp_pts, 
         dsn = paste0(out_dir, site, year, "valid_pts"), 
         layer = "valid_pts",
         driver = "ESRI Shapefile",
         overwrite = TRUE)
