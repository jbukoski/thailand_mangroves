# Builds LUC plots

library(tidyverse)
library(ggthemes)
library(raster)
library(sf)
library(magrittr)
library(gridExtra)

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/site_map/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "nakorn_"

# Read in rasters

thailand <- read_sf(paste0(in_dir, "thailand_boundary.shp"))

svm1987 <- raster(paste0(out_dir, site, "1987_svm.tif"))
svm1997 <- raster(paste0(out_dir, site, "1997_svm.tif"))
svm2007 <- raster(paste0(out_dir, site, "2007_svm.tif"))
svm2017 <- raster(paste0(out_dir, site, "2017_svm.tif"))

# Produce data frames for plotting

svm1987_df <- svm1987 %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame() %>%
  set_colnames(c("class", "x", "y"))

svm1997_df <- svm1997 %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame() %>%
  set_colnames(c("class", "x", "y"))

svm2007_df <- svm2007 %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame() %>%
  set_colnames(c("class", "x", "y"))

svm2017_df <- svm2017 %>%
  as("SpatialPixelsDataFrame") %>%
  as.data.frame() %>%
  set_colnames(c("class", "x", "y"))

# Build plots

p1 <- ggplot() +
  geom_sf(data = thailand, fill = "#f9f0f9") +
  geom_tile(data = svm1987_df, aes(x=x, y=y, fill=class)) +
  theme_tufte() +
  scale_fill_viridis_c() +
  xlim(c(98.82, 99.05)) +
  ylim(c(7.825, 8.12)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(x = 99.02, y = 8.11, label= "1989"))

p2 <- ggplot() +
  geom_sf(data = thailand, fill = "#f9f0f9") +
  geom_tile(data = svm1997_df, aes(x=x, y=y, fill=class)) +
  theme_tufte() +
  scale_fill_viridis_c() +
  xlim(c(98.82, 99.05)) +
  ylim(c(7.825, 8.12)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(x = 99.02, y = 8.11, label= "1997"))

p3 <- ggplot() +
  geom_sf(data = thailand, fill = "#f9f0f9") +
  geom_tile(data = svm2007_df, aes(x=x, y=y, fill=class)) +
  theme_tufte() +
  scale_fill_viridis_c() +
  xlim(c(98.82, 99.05)) +
  ylim(c(7.825, 8.12)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(x = 99.02, y = 8.11, label= "2007"))

p4 <- ggplot() +
  geom_sf(data = thailand, fill = "#f9f0f9") +
  geom_tile(data = svm2017_df, aes(x=x, y=y, fill=class)) +
  theme_tufte() +
  scale_fill_viridis_c() +
  xlim(c(98.82, 99.05)) +
  ylim(c(7.825, 8.12)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(x = 99.02, y = 8.11, label= "2017"))

# Grid plots

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2, respect = TRUE, padding = 0)
