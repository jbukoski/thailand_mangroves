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


# Code to convert Rich's Trang data
# 
# library(rgdal)
# library(tidyverse)
# library(measurements)
# library(leaflet)
# 
# dat <- read_csv("/home/jbukoski/Desktop/plots.csv", col_names = TRUE) %>%
#   mutate(latitude = gsub("['\"^0]", "", latitude),
#          longitude = gsub("['\"]", "", longitude)) %>%
#   mutate(latitude = gsub("[?]", " ", latitude),
#          longitude = gsub("[?']", " ", longitude)) %>%
#   filter(grepl("^[0-9]+ [0-9]+.[0-9]+$", latitude)) %>%
#   mutate(latitude = as.numeric(measurements::conv_unit(latitude, from = "deg_dec_min", to = "dec_deg")),
#          longitude = as.numeric(measurements::conv_unit(longitude, from = "deg_dec_min", to = "dec_deg")))
# 
# xy <- dat[, 1:2]
# 
# write_csv(dat, "/home/jbukoski/Desktop/plots_adjusted.csv")
# 
# spdf <- SpatialPointsDataFrame(coords = xy, data = dat)
# 
# leaflet() %>%
#   addTiles() %>%
#   addMarkers(dat$longitude, dat$latitude)
