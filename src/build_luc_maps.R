# Builds LUC plots

library(tidyverse)
library(ggthemes)
library(raster)
library(sf)

in_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

pred_1987_4326 <- projectRaster(pred_1987, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
pred_1987_4326_df <- as(pred_1987_4326, "SpatialPixelsDataFrame") %>%
  as.data.frame() %>%
  mutate(layer = (round(layer, 0)))

p1 <- ggplot() +
  geom_sf(data = thailand, fill = "#f9f0f9") +
  geom_tile(data = pred_1987_4326_df, aes(x=x, y=y, fill=layer)) +
  theme_tufte() +
  scale_fill_viridis_c() +
  xlim(c(98.82, 99.05)) +
  ylim(c(7.825, 8.12)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p2 <- ggplot() +
  geom_sf(data = pred_1997) +
  theme_tufte() +
  theme(legend.position = "none")

p3 <- ggplot(pred_2007) +
  geom_sf(data = pred_2007) +
  theme_tufte() +
  theme() +
  theme(legend.position = "none")

p4 <- ggplot(pred_2007) +
  geom_sf(data = pred_2007) +
  theme_tufte() +
  theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2, respect = TRUE, padding = 0)