# Small script to convert the projection of the plot coordinates to
# WGS84 Projection

library(xlsx)
library(rgdal)
library(tidyverse)
library(magrittr)
library(readxl)

# Specify in/our directories

in_dir <- '/home/jbukoski/research/data/thailand_stocks/input/'
out_dir <- '/home/jbukoski/research/data/thailand_stocks/output/'

# Specify site for which you are working with

site <- 'thailand'

# Load data

dat <- read_excel(paste0(in_dir, site, '_data.xlsx'),
                  sheet = 'Metadata', col_names = TRUE) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub(" ", "_", colnames(.))) %>%
  set_colnames(gsub("[(),]", "", colnames(.))) %>%
  select(longitude, latitude) %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude))

# Convert data

utm_coords <- SpatialPoints(cbind(dat$latitude, dat$longitude), 
                            proj4string = CRS("+proj=utm +zone=47N")) 

wgs84_coords <- spTransform(utm_coords, CRS("+proj=longlat")) %>%
  as_data_frame() %>%
  set_colnames(c("longitude", "latitude"))

# Write data to file

write_csv(wgs84_coords, 
          paste0(out_dir, 'plot_coords.csv'))
