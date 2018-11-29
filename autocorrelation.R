library(lme4)
library(magrittr)
library(raster)
library(readxl)
library(sf)
library(spdep)
library(tidyverse)

#----------------------
# Load in helper functions & allometry functions

source("/home/jbukoski/research/scripts/thailand_stocks/src/helper_funcs.R")
source("/home/jbukoski/research/scripts/thailand_stocks/src/allometry.R")
source("/home/jbukoski/research/scripts/thailand_stocks/src/utm_to_latlong.R")

#----------------------
# Specify in/out directories

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

meta <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                   sheet="Metadata", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. ]", "_", colnames(.))) %>%
  select(-date, -crew_members, -data_checked_by, -data_entered_by,
         -accuracy, -dominant_species, -n_photo, -e_photo, -w_photo,
         -s_photo, -topography, -position, -disturbance, - other_notes,
         -`fringe,_transition,_interior`)

meta_clean <- read_sf(paste0(out_dir, "meta.shp"))


raw_trees <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                        sheet="Trees", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. ]", "_", colnames(.))) %>%
  select(-date, -recorder, -checked_by, -entered_by)

raw_saps <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                       sheet="Saplings", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. ]", "_", colnames(.))) %>%
  select(-date, -recorder, -checked_by, -entered_by)

raw_soil <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                       sheet="Soil", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. ]", "_", colnames(.))) %>%
  select(-date, -recorder, -checked_by, -entered_by,
         -actual_depth_start, -actual_depth_stop, -salinity_ppt,
         -soil_color, -texture, -depth_1, -depth_2, -depth_3)

#-----------------------------------------
# Specify necessary parameters

plot_size <- (7^2)*pi
subplot_size <- (2^2)*pi
transect_size <- 5 * plot_size

t_val <- qt(0.975, 7-1)

# Specify site areas in (what units?)

site_areas <- tibble(site = c("Krabi", "Nakorn"),
                     area = c(100000, 100000))

#------------------------------------------------------------------------------

#------------------------#
# Processing of raw data #
#------------------------#

# Create species code and calculate basal area

trees <- raw_trees %>% 
  id_taxon(raw_trees$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2)),
         basal_area = pi * (dbh_cm/100/2)^2) %>%
  dplyr::select(-genus, -species)

# Calculate above-ground biomass using species-specific allometric equations. 
# Where species-specific equations are not available, used Komiyama et al 2005 
# general equation with species specific wood densities

trees <- trees %>%
  left_join(allom_lookup, by = c("sps_code" = "sps_code")) %>%
  mutate(params = map2(dbh_cm, density, list)) %>%
  mutate(agb = invoke_map_dbl(ag_form, params)) %>%
  mutate(bgb = invoke_map_dbl(bg_form, params))

# Adjust AGB variable based on 'status' variable
# Calculate cone if base_cm measurement exists, otherwise assume a cylinder

trees <- trees %>%
  mutate(top_diam = ifelse(status == 3 & !is.na(base_cm) & height_m < 1.37, 
                           base_cm - (100 * height_m * ((base_cm - dbh_cm) / (100 * height_m))),
                           base_cm - (100 * height_m * ((base_cm - dbh_cm) / 137)))) %>%
  mutate(stump.vol = ifelse(status == 3 & !is.na(base_cm) & height_m < 1.37, 
                            (pi * 100 * height_m) / 12 * (base_cm^2 + top_diam^2 + (base_cm * top_diam)),
                            ifelse(status == 3 & !is.na(base_cm) & height_m >= 1.37, 
                                   (pi * 100 * height_m) / 12 * (base_cm^2 + top_diam^2 + (base_cm * top_diam)), 
                                   NA))) %>%
  mutate(adj_agb = ifelse(is.na(status), agb,
                          ifelse(status == 1, 0.95 * agb, 
                                 ifelse(status == 2, 0.8 * agb, 
                                        ifelse(status == 3 & !is.na(base_cm), 
                                               density * stump.vol / 1000, 
                                               density * pi * dbh_cm * height_m * 100 / 1000)))))

trees <- trees %>%
  mutate(agb = adj_agb,
         biomass = adj_agb + bgb) %>%
  dplyr::select(-ag_form, -bg_form, -ag_ref, -bg_ref, 
                -params, -top_diam, -stump.vol, - adj_agb)

#-------------------------------------------------------------------------------
# Processing for saplings
# Follows same steps as for trees

saps <- raw_saps %>% 
  id_taxon(raw_saps$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2))) %>%
  dplyr::select(-genus, -species)

saps <- saps %>%
  left_join(allom_lookup, by = c("sps_code" = "sps_code")) %>%
  mutate(params = map2(dbh_cm, density, list)) %>%
  mutate(agb = invoke_map_dbl(ag_form, params)) %>%
  mutate(bgb = invoke_map_dbl(bg_form, params))

saps <- saps %>%
  mutate(adj_agb = ifelse(is.na(status), agb,
                          ifelse(status == 1, 0.95 * agb, 0.8 * agb)))

saps <- saps %>%
  mutate(agb = adj_agb,
         biomass = agb + bgb) %>%
  dplyr::select(-ag_form, -bg_form, -ag_ref, -bg_ref, -params, - adj_agb)

#-------------------------------------------------------------
# Biomass: Join trees and saplings

biomass <- bind_rows(mutate(trees, stage = "tree"), mutate(saps, stage = "sapling"))

#-------------------------------------------------------------------------------
# Compute tree and sapling based forest structure variables.

trees_structure <- trees %>%
  select(site, plot, subplot, dbh_cm, status, sps_code, basal_area) %>%
  group_by(site, plot, subplot) %>%
  mutate(subplot_dbh = mean(dbh_cm),
         subplot_ba = 10000 * sum(basal_area) / plot_size,
         subplot_n = 10000 * n() / plot_size) %>%
  select(site, plot, subplot, subplot_dbh, subplot_ba, subplot_n) %>%
  distinct %>%
  group_by(site, plot) %>%
  mutate(plot_dbh = mean(subplot_dbh),
         plot_dbh_se = sqrt(var(subplot_dbh) / n()),
         plot_ba = mean(subplot_ba),
         plot_ba_se = sqrt(var(subplot_ba) / n()),
         plot_n = mean(subplot_n),
         plot_n_se = sqrt(var(subplot_n) / n())) %>%
  select(site, plot, subplot, subplot_dbh, subplot_ba, subplot_n,
         plot_dbh, plot_dbh_se, plot_ba, plot_ba_se,
         plot_n, plot_n_se) %>%
  distinct


plot_biomass <- biomass %>%
  select(site, plot, subplot, biomass, agb, bgb) %>%
  left_join(site_areas, by = "site") %>%
  group_by(site, plot, subplot) %>%
  summarise(ttl_tau = mean(area) * ( sum(biomass) / plot_size),
            agb_tau = mean(area) * ( sum(agb) / plot_size),
            bgb_tau = mean(area)* ( sum(bgb) / plot_size)) %>%
  group_by(site, plot, subplot) %>%
  left_join(site_areas, by = "site") %>%
  mutate(avg_plot_ttl = mean(ttl_tau),
         avg_plot_agb = mean(agb_tau),
         avg_plot_bgb = mean(bgb_tau),
         avg_plot_ttl_se = sqrt(var(ttl_tau) / n()),
         avg_plot_agb_se = sqrt(var(agb_tau) / n()),
         avg_plot_bgb_se = sqrt(var(bgb_tau) / n()),
         plot_ttl_ha = 10 * avg_plot_ttl / area,
         plot_agb_ha = 10 * avg_plot_agb / area,
         plot_bgb_ha = 10 * avg_plot_bgb / area,
         plot_ttl_ha_se = 10 * avg_plot_ttl_se / area,
         plot_agb_ha_se = 10 * avg_plot_agb_se / area,
         plot_bgb_ha_se = 10 * avg_plot_bgb_se / area) %>%
  select(site, plot, subplot, plot_ttl_ha, plot_ttl_ha_se,
         plot_agb_ha, plot_agb_ha_se,
         plot_bgb_ha, plot_bgb_ha_se) %>%
  distinct

#--------------------------------------------------
#Species classification

sps_rii_init <- trees %>%
  group_by(site, plot, subplot) %>%
  mutate(site_n = n(), 
         site_ba = sum(basal_area)) %>%
  group_by(site, plot, subplot, sps_code) %>%
  mutate(sps_ba = sum(basal_area),
         sps_n = n(),
         rel_sps_ba = mean(sps_ba)/mean(site_ba) * 100,
         rel_sps_n = mean(sps_n)/mean(site_n) * 100,
         rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  select(site, plot, subplot, sps_code, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(site, plot, subplot, -rel_imp) %>%
  distinct

# Compute relative frequency

subplot_sps_freq <- trees %>%
  group_by(site, plot, subplot, sps_code) %>%
  add_count(sps_code) %>%
  rename(sps_n = n) %>%
  select(site, plot, subplot, sps_code, sps_n) %>%
  distinct

subplot_freq <- trees %>% group_by(site, plot, subplot) %>%
  add_count() %>%
  select(site, plot, subplot, n) %>%
  distinct

freq <- subplot_sps_freq %>%
  left_join(subplot_freq, by = c("site", "plot", "subplot")) %>%
  mutate(rel_freq = sps_n / n * 100)

# Join species importance indexes

sps_rii <- sps_rii_init %>%
  left_join(freq, by = c("site", "plot", "subplot", "sps_code")) %>%
  select(site, plot, subplot, sps_code, genus, species, rel_sps_n, rel_freq, rel_sps_ba) %>%
  distinct %>%
  rowwise %>%
  mutate(rel_imp = sum(rel_sps_n, rel_freq, rel_sps_ba) / 3,
         subplot_id = paste0(substring(site, 1, 2), plot, subplot)) %>%
  arrange(subplot_id, -rel_imp)

# Determine subplot sps groups

idx <- aggregate(rel_imp ~ subplot_id, data = sps_rii, max)

dom_sps <- sps_rii %>%
  filter(rel_imp %in% idx$rel_imp) %>%
  mutate(dom_sps = ifelse(rel_imp > 80, paste0(sps_code), paste0(sps_code, "_mix")))

#--------------------------------------------------
# Plot classification by genus

genus_rii_init <- trees %>%
  group_by(site, plot, subplot) %>%
  mutate(site_n = n(), 
         site_ba = sum(basal_area)) %>%
  group_by(site, plot, subplot, genus) %>%
  mutate(sps_ba = sum(basal_area),
         sps_n = n(),
         rel_sps_ba = mean(sps_ba)/mean(site_ba) * 100,
         rel_sps_n = mean(sps_n)/mean(site_n) * 100,
         rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  select(site, plot, subplot, genus, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(site, plot, subplot, -rel_imp) %>%
  distinct

# Compute relative frequency

subplot_genus_freq <- trees %>%
  group_by(site, plot, subplot, genus) %>%
  add_count(genus) %>%
  rename(sps_n = n) %>%
  select(site, plot, subplot, genus, sps_n) %>%
  distinct

subplot_freq <- trees %>% group_by(site, plot, subplot) %>%
  add_count() %>%
  select(site, plot, subplot, n) %>%
  distinct

freq <- subplot_genus_freq %>%
  left_join(subplot_freq, by = c("site", "plot", "subplot")) %>%
  mutate(rel_freq = sps_n / n * 100)

# Join species importance indexes

genus_rii <- genus_rii_init %>%
  left_join(freq, by = c("site", "plot", "subplot", "genus")) %>%
  select(site, plot, subplot, genus, rel_sps_n, rel_freq, rel_sps_ba) %>%
  distinct %>%
  rowwise %>%
  mutate(rel_imp = sum(rel_sps_n, rel_freq, rel_sps_ba) / 3,
         subplot_id = paste0(substring(site, 1, 2), plot, subplot)) %>%
  arrange(subplot_id, -rel_imp)

# Determine subplot sps groups

idx <- aggregate(rel_imp ~ subplot_id, data = genus_rii, max)

dom_genus <- genus_rii %>%
  filter(rel_imp %in% idx$rel_imp) %>%
  mutate(dom_genus = ifelse(rel_imp > 80, paste0(substring(genus, 1, 2)), paste0(substring(genus, 1, 2), "_mix")))


#--------------------------------------------------
# Processing soil

soil <- raw_soil %>%
  mutate(c_dens = bulk_density * (percent_c/100),
         int_volume = ifelse(interval == 5, 
                             ((avg_depth/100) - (100/100)) * 10000, 
                             ((int_b/100) - (int_a/100)) * 10000),
         soc_per_ha = int_volume * c_dens)

soil_summary <- soil %>%
  group_by(site, plot, subplot) %>%
  mutate(subplot_soc = sum(soc_per_ha, na.rm = TRUE),
         subplot_bd = mean(bulk_density),
         subplot_c_dens = mean(c_dens),
         subplot_poc = mean(percent_c)) %>%
  select(site, plot, subplot, subplot_bd,
         subplot_poc, subplot_c_dens, subplot_soc) %>%
  distinct()


#-------------------------------------------------------------------------------
# Join to spatial

full_data <- meta_clean %>% 
  left_join(trees_structure, by = c("site", "plot", "subplot")) %>%
  left_join(plot_biomass, by = c("site", "plot", "subplot")) %>%
  left_join(soil_summary, by = c("site", "plot", "subplot")) %>%
  select(-plot_ttl_ha_se, -plot_agb_ha_se, -plot_bgb_ha_se,
         -plot_dbh, -plot_dbh_se, -plot_ba, -plot_ba_se, -plot_n, -plot_n_se) %>%
  rename(subplot_ttl_ha = plot_ttl_ha, subplot_agb_ha = plot_agb_ha, subplot_bgb_ha = plot_bgb_ha) %>%
  mutate(subplot_id = paste0(substring(site, 1, 2), plot, subplot),
         plot_id = paste0(substring(site, 1,2), plot)) %>%
  left_join(select(dom_sps, rel_imp, site, plot, subplot, dom_sps), by = c("site", "plot", "subplot")) %>%
  rename(sps_rii = rel_imp) %>%
  left_join(select(dom_genus, rel_imp, site, plot, subplot, dom_genus), by = c("site", "plot", "subplot")) %>%
  rename(genus_rii = rel_imp) 


#writeOGR(as(full_data, "Spatial"), paste0(out_dir), layer = "processed_subplot_data.shp", driver = "ESRI Shapefile")
#write_csv(full_data, paste0(out_dir, "processed_subplot_data.csv"))

#-------------------------------------------------------------------------------

coords <- full_data %>%
  filter(site == "Nakorn", plot != 3) %>%
  drop_na() %>%
  as("Spatial")
  
nb <- knn2nb(knearneigh(coords, k = round(nrow(coords)*0.3, 0)))

summary(nb)

nbm <- nb2mat(nb, style = "B")
nbw <- nb2listw(nb, style='B')
nbw

mc <- moran.mc(coords$subplot_ttl_ha, nbw, nsim=99)
mc

plot(mc)

full_data %>% 
  #filter(site == "Krabi") %>% 
  ggplot(aes(subplot_ttl_ha, subplot_poc, group = plot)) + 
  geom_point()

#------------------------------------------
# Mixed model

null <- lme4::lmer(subplot_soc ~ subplot_ttl_ha + (1|plot), 
                   data = drop_na(krabi), REML = FALSE)

lmm <- lme4::lmer(subplot_soc ~ subplot_ttl_ha + dom_sps + (1|plot), 
                  data = drop_na(krabi), REML = FALSE)

anova(null, lmm)
