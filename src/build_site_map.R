# Script for creation of a site map for stocks study

library(cowplot)
library(ggmap)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(rgdal)
library(sf)



in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/site_map/"

se_asia <- read_sf(paste0(in_dir, "se_asia.shp"))
thailand <- read_sf(paste0(in_dir, "thailand_boundary.shp"))
mangroves <- read_sf(paste0(in_dir, "thailand_mangroves.shp"))
krabi <- read_sf(paste0(in_dir, "krabi_admin.shp"))
nakorn <- read_sf(paste0(in_dir, "nakorn_admin.shp"))
trang <- read_sf(paste0(in_dir, "trang_admin.shp"))
kre <- read_sf(paste0(in_dir, "kre_boundary.shp"))
ppm <- read_sf(paste0(in_dir, "ppm_boundary.shp"))
pre <- read_sf(paste0(in_dir, "pre_boundary.shp"))

p1 <- ggplot(se_asia) +
  geom_sf(data = se_asia, fill = "#f9f0f9") +
  geom_sf(data = thailand, fill = "#cacaca") +
  geom_sf(data = mangroves, fill = "dark green", col = "dark green") +
  geom_sf(data = krabi, col = "black", fill = NA, lwd = 0.25, linetype = "dashed") +
  geom_sf(data = nakorn, col = "black", fill = NA, lwd = 0.25, linetype = "dashed") +
  geom_sf(data = trang, col = "black", fill = NA, lwd = 0.25, linetype = "dashed") +
  #geom_sf(data = kre, fill = NA, col = "red", lwd = 1.1) +
  #geom_sf(data = ppm, fill = NA, col = "red", lwd = 1.1) +
  #geom_sf(data = pre, fill = NA, col = "red", lwd = 1.1) +
  theme_tufte() +
  xlim(c(95, 107)) + 
  ylim(c(5, 22)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))

p2 <- ggplot(thailand) +
  geom_sf(data = thailand, fill = "#cacaca") +
  geom_sf(data = krabi, col = "black", fill = NA, lwd = 0.25, linetype = "dashed") +
  geom_sf(data = nakorn, col = "black", fill = NA, lwd = 0.25, linetype = "dashed") +
  geom_sf(data = trang, col = "black", fill = NA, lwd = 0.25, linetype = "dashed") +
  geom_sf(data = mangroves, fill = "dark green", col = "dark green") +
  geom_sf(data = kre, fill = NA, col = "red", lwd = 1.1) +
  geom_sf(data = ppm, fill = NA, col = "red", lwd = 1.1) +
  geom_sf(data = pre, fill = NA, col = "red", lwd = 1.1) +
  geom_text(y = 8.2,  x = 99, label = "A", size = 4) +
  geom_text(y = 8.2,  x = 100.1, label = "B", size = 4) +
  geom_text(y = 7.5,  x = 99.5, label = "C", size = 4) +
  theme_tufte() +
  xlim(c(97.75, 100.75)) +
  ylim(c(7.0, 9.3)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))
  
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

grid.arrange(g1, g2, ncol = 2, respect = TRUE, widths = c(0.75, 1.29))

  