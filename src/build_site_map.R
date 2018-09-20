# Script for creation of a site map for stocks study

library(ggmap)
library(ggplot2)
library(rgdal)
library(sf)
library(ggthemes)
library(cowplot)

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/site_map/"

se_asia <- read_sf(paste0(in_dir, "se_asia.shp"))
thailand <- read_sf(paste0(in_dir, "thailand_boundary.shp"))
mangroves <- read_sf(paste0(in_dir, "thailand_mangroves.shp"))
krabi <- read_sf(paste0(in_dir, "krabi_admin.shp"))
nakorn <- read_sf(paste0(in_dir, "nakorn_admin.shp"))
kre <- read_sf(paste0(in_dir, "kre_boundary.shp"))
ppm <- read_sf(paste0(in_dir, "ppm_boundary.shp"))

p1 <- ggplot(se_asia) +
  geom_sf(data = se_asia, fill = "#f9f0f9") +
  geom_sf(data = thailand, fill = "#cacaca") +
  geom_sf(data = krabi, fill = "yellow", lwd = 0.1) +
  geom_sf(data = nakorn, fill = "yellow", lwd = 0.1) +
  geom_sf(data = mangroves, fill = "dark green", col = "dark green") +
  theme_tufte() +
  xlim(c(90, 110)) + 
  ylim(c(0, 25)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = unit(c(0, -1.5, 0, 0), "cm"))
  
p1

p2 <- ggplot(thailand) +
  geom_sf(data = thailand, fill = "#cacaca") +
  geom_sf(data = krabi, fill = "yellow", lwd = 0.1) +
  geom_sf(data = nakorn, fill = "yellow", lwd = 0.1) +
  geom_sf(data = mangroves, fill = "dark green", col = "dark green") +
  geom_sf(data = kre, fill = NA, col = "red", lwd = 1.1) +
  geom_sf(data = ppm, fill = NA, col = "red", lwd = 1.1) +
  theme_tufte() +
  xlim(c(97.75, 100.75)) +
  ylim(c(7.5, 9.25)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = unit(c(0, 0, 0, -1), "cm"))
  
p2

grid.arrange(p1, p2, nrow = 1, respect = TRUE, widths = c(1, 2), padding = 0)

plot_grid(p1, p2, align = "h", rel_widths = c(1, 2))

arrangeGrob(p1, p2)
