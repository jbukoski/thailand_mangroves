#######################################################
## Site analysis for Krabi and Pak Panang field data ##
#######################################################
#
# Imports plotwise field data for Krabi River Estuary
# and the Pak Panang Mangrove and generates site wide
# estimates of forest structure, species composition
# and ecosystem carbon stocks
#
# Inputs:
#  1. Excel file containing field data (publish on Harvard DataVerse later?)
#
# Outputs
#  2. CSV files of processed data
#
#
#--------------------------------------
# Load libraries and begin script

library("ggplot2")
library("lmfor")
library("magrittr")
library("tidyverse")
library("readxl")

#----------------------
# Load in helper functions & allometry functions

source("/home/jbukoski/research/scripts/thailand_stocks/src/helper_funcs.R")
source("/home/jbukoski/research/scripts/thailand_stocks/src/allometry.R")

#----------------------
# Specify in/out directories

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

#---------------------------------------
## Load in data

meta <- read_excel(paste0(in_dir, "thailand_data.xlsx"), sheet="Metadata", col_names = T)

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

raw_seedlings <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                        sheet="Seedlings", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. ]", "_", colnames(.))) %>%
  select(-date, -recorder, -checked_by, -entered_by)

raw_cwd <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                  sheet="CWD", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[ .]", "_", colnames(.))) %>%
  select(-date, -recorder, -checked_by, -entered_by)

raw_soil <- read_excel(paste0(in_dir, "thailand_data.xlsx"), 
                   sheet="Soil", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. ]", "_", colnames(.))) %>%
  select(-date, -recorder, -checked_by, -entered_by,
         -actual_depth_start, -actual_depth_stop, -salinity_ppt,
         -soil_color, -texture, -depth_1, -depth_2, -depth_3)

raw_aqua <- read_excel(paste0(in_dir, "aquaculture_data.xlsx"),
                       sheet="Soil", col_names = T) %>%
  set_colnames(tolower(colnames(.))) %>%
  set_colnames(gsub("[. /]", "_", colnames(.))) %>%
  select(site, site_category, plot = angie_transect_id,
         interval = d_interval, soil_depth = mean_soil__depth_cm,
         bulk_density = soil_bulk_density_g_cm3, percent_c = soil_perc_c, c_dens = carbon_density)

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

#-----------------------------------------------------------------------------
# Processing for coarse woody debris (CWD)

# Convert CWD measurements to mass C in Mg/ha following Kauffman & Donato 2012
# Mean specific gravities (g/cm^3) of CWD classes taken from K&D, 2012

cwd_params <- tibble(size = c("fine", "small", "medium", "large"),
                     density = c(0.48, 0.64, 0.71, 0.69),
                     avg_diam = c(0.43, 1.47, 4.52, NA))

cwd <- raw_cwd %>%
  dplyr::select(-remarks) %>%
  gather(size, n, -site, -plot, -subplot, -transect) %>%
  separate(size, c("size", "status")) %>%
  left_join(cwd_params, by = "size") %>%
  mutate(trnsct_lngth = ifelse(size == "fine", 2,
                               ifelse(size == "small", 3,
                                      ifelse(size == "medium", 5, 12)))) %>%
  mutate(volume = ifelse(size != "large", (pi^2) * ((n * (avg_diam^2)) / (8 * trnsct_lngth)),
                         (pi^2) * (n^2) / (8 * trnsct_lngth)),
         mass = volume * density,
         adj_mass = ifelse(is.na(status), mass,
                           ifelse(status == "rotten", mass * 0.5, mass))) %>%
  group_by(site, plot, subplot, size, status) %>%
  summarise(total = mean(n), mass = mean(mass)) %>%
  group_by(site, plot, subplot) %>%
  mutate(subplot_mass = sum(mass)) %>%
  group_by(site, plot) %>%
  mutate(plot_mass = mean(subplot_mass))

#------------------------------------------------------
# Processing for soil

soil <- raw_soil %>%
  mutate(c_dens = bulk_density * (percent_c/100),
         int_volume = ifelse(interval == 5, 
                             ((avg_depth/100) - (100/100)) * 10000, 
                             ((int_b/100) - (int_a/100)) * 10000),
         soc_per_ha = int_volume * c_dens)

#-------------------------------------------------------------------------------

#----------------------------------------#
# Forest structure & species composition #
#----------------------------------------#

# Compute tree and sapling based forest structure variables.

trees_structure <- trees %>%
  select(site, plot, subplot, dbh_cm, status, sps_code, basal_area) %>%
  group_by(site, plot, subplot) %>%
  mutate(subplot_dbh = mean(dbh_cm),
         subplot_ba = sum(basal_area),
         subplot_n = n()) %>%
  select(site, plot, subplot, subplot_dbh, subplot_ba, subplot_n) %>%
  distinct %>%
  group_by(site, plot) %>%
  mutate(plot_dbh = mean(subplot_dbh),
         plot_dbh_se = sqrt(var(subplot_dbh) / n()),
         plot_ba = 10000 * mean(subplot_ba) / plot_size,
         plot_ba_se = 10000 * sqrt(var(subplot_ba) / n()) / plot_size,
         plot_n = 10000 * mean(subplot_n) / plot_size,
         plot_n_se = 10000 * sqrt(var(subplot_n) / n()) / plot_size) %>%
  select(site, plot, plot_dbh, plot_dbh_se, plot_ba, plot_ba_se,
         plot_n, plot_n_se) %>%
  distinct %>%
  group_by(site) %>%
  mutate(site_dbh = mean(plot_dbh),
         site_dbh_se = sqrt(var(plot_dbh) / n()),
         site_ba = mean(plot_ba),
         site_ba_se = sqrt(var(plot_ba) / n()),
         site_n = mean(plot_n),
         site_n_se = sqrt(var(plot_n) / n())) %>%
  select(1, 9:14) %>%
  distinct

saps_structure <- saps %>%
  select(site, plot, subplot, dbh_cm, status, sps_code) %>%
  group_by(site, plot, subplot) %>%
  mutate(subplot_dbh = mean(dbh_cm),
         subplot_n = n()) %>%
  select(site, plot, subplot, subplot_dbh, subplot_n) %>%
  distinct %>%
  group_by(site, plot) %>%
  mutate(plot_dbh = mean(subplot_dbh),
         plot_dbh_se = sqrt(var(subplot_dbh) / n()),
         plot_n = 10000 * mean(subplot_n) / plot_size,
         plot_n_se = 10000 * sqrt(var(subplot_n) / n()) / plot_size) %>%
  select(site, plot, plot_dbh, plot_dbh_se, plot_n, plot_n_se) %>%
  distinct

# Compute relative importance index for each species at the site level
# Compute rel density and rel dominance first, rel frequency later

sps_rii <- trees %>%
  group_by(site) %>%
  mutate(site_n = n(), 
         site_ba = sum(basal_area)) %>%
  group_by(site, sps_code) %>%
  mutate(sps_ba = sum(basal_area),
         sps_n = n(),
         rel_sps_ba = mean(sps_ba)/mean(site_ba) * 100,
         rel_sps_n = mean(sps_n)/mean(site_n) * 100,
         rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  select(site, sps_code, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(site, -rel_imp)

# Compute relative frequency

freq <- trees %>%
  select(site, plot, sps_code) %>%
  distinct %>%
  mutate(presence = 1)

freq_frame <- trees %>%
  distinct(site, plot) %>%
  left_join(expand(trees, plot, sps_code), by = "plot") %>%
  left_join(freq, by = c("site", "plot", "sps_code")) %>%
  mutate(presence = ifelse(is.na(presence), 0, presence))

rel_freq <- freq_frame %>%
  group_by(site, sps_code) %>%
  summarise(sps_n = sum(presence),
         ttl_n = n(),
         rel_freq = 100 * sps_n / ttl_n) %>%
  select(site, sps_code, rel_freq) %>%
  filter(rel_freq > 0) %>%
  arrange(site, -rel_freq)

# Join species importance indexes

sps_rii <- sps_rii %>%
  left_join(rel_freq, by = c("site", "sps_code")) %>%
  select(site, sps_code, genus, species, rel_sps_n, rel_freq, rel_sps_ba) %>%
  distinct %>%
  rowwise %>%
  mutate(rel_imp = sum(rel_sps_n, rel_freq, rel_sps_ba) / 3) %>%
  arrange(site, -rel_imp)

#-------------------------------------------------------------------------------

#----------------------------#
# Carbon & biomass estimates #
#----------------------------#

# Obtain estimate of biomass per hectare based on all plots
# Outputs estimates of biomass in Mg/ha at the subplot, plot, and site aggregation
# Any of the trees, saplings, or biomass tibbles can be run through
# Calculation of mean and variance follows Gregoire and Valentine 20XX?

plot_biomass <- biomass %>%
  select(site, plot, subplot, biomass, agb, bgb) %>%
  left_join(site_areas, by = "site") %>%
  group_by(site, plot, subplot) %>%
  summarise(ttl_tau = mean(area) * ( sum(biomass) / plot_size),
            agb_tau = mean(area) * ( sum(agb) / plot_size),
            bgb_tau = mean(area)* ( sum(bgb) / plot_size)) %>%
  group_by(site, plot) %>%
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
  select(site, plot, plot_ttl_ha, plot_ttl_ha_se,
         plot_agb_ha, plot_agb_ha_se,
         plot_bgb_ha, plot_bgb_ha_se) %>%
  distinct

site_biomass <- plot_biomass %>%
  group_by(site) %>%
  select(site, plot_ttl_ha, plot_agb_ha, plot_bgb_ha) %>%
  distinct %>%
  mutate(site_ttl_ha = mean(plot_ttl_ha),
         site_agb_ha = mean(plot_agb_ha),
         site_bgb_ha = mean(plot_bgb_ha),
         site_ttl_ha_se = sqrt(var(plot_ttl_ha) / n()),
         site_agb_ha_se = sqrt(var(plot_agb_ha) / n()),
         site_bgb_ha_se = sqrt(var(plot_bgb_ha) / n())) %>%
  select(site, site_ttl_ha, site_ttl_ha_se, site_agb_ha, 
         site_agb_ha_se, site_bgb_ha, site_bgb_ha_se) %>%
  distinct

# Print summary table of total biomass & aboveground vs belowground carbon
# Note different units (biomass vs C) of reported values.

biomass_c_summary <- plot_biomass %>%
  select(plot, site, plot_ttl_ha, plot_ttl_ha_se, 
         plot_agb_ha, plot_agb_ha_se, 
         plot_bgb_ha, plot_bgb_ha_se) %>%
  mutate(plot_agb_c_ha = plot_agb_ha * 0.47,
         plot_agb_c_ha_se = plot_agb_ha_se * 0.47,
         plot_bgb_c_ha = plot_bgb_ha * 0.39,
         plot_bgb_c_ha_se = plot_bgb_ha_se * 0.39) %>%
  select(-plot_agb_ha, -plot_agb_ha_se, -plot_bgb_ha, -plot_bgb_ha_se) %>%
  distinct()

# Print summary table of CWD carbon estimates (mg / ha) using processed CWD table

cwd_summary <- cwd %>%
  select(site, plot, subplot, subplot_mass, plot_mass) %>%
  distinct() %>%
  group_by(site, plot) %>%
  mutate(plot_mass_var = var(subplot_mass),
         plot_mass_se = sqrt(var(subplot_mass) / n()),
         n = n()) %>%
  ungroup() %>%
  select(site, plot, plot_mass, plot_mass_se) %>%
  distinct() %>%
  mutate(plot_mass_c = plot_mass * 0.5,
         plot_mass_c_se = plot_mass_se * 0.5)


#------------------------------------------------------------------------------

#-----------------------------#
# Analysis of soil properties #
#-----------------------------#

soil_summary <- soil %>%
  group_by(site, plot, subplot) %>%
  mutate(subplot_c = sum(soc_per_ha)) %>%
  select(site, plot, subplot, subplot_c) %>%
  distinct() %>%
  group_by(site, plot) %>%
  mutate(plot_c = mean(subplot_c),
         plot_c_se = sqrt(var(subplot_c) / n())) %>%
  select(site, plot, plot_c, plot_c_se) %>%
  distinct

#-----------------------
# Join all carbon pools

c_summary <- bind_cols(site = biomass_c_summary$site,
                       plot = biomass_c_summary$plot,
                       agc = biomass_c_summary$plot_agb_c_ha,
                       bgc = biomass_c_summary$plot_bgb_c_ha,
                       cwd = cwd_summary$plot_mass_c,
                       soc = soil_summary$plot_c,
                       agc_se = biomass_c_summary$plot_agb_c_ha_se,
                       bgc_se = biomass_c_summary$plot_bgb_c_ha_se,
                       cwd_se = cwd_summary$plot_mass_c_se,
                       soc_se = soil_summary$plot_c_se) %>%
  mutate(total = sum(agc, bgc, cwd, soc),
         total_se = sqrt(agc_se^2 + bgc_se^2 + cwd_se^2 + soc_se^2)) %>%
  group_by(site) %>%
  mutate(agc_avg = mean(agc),
         agc_se = sqrt(var(agc) / n()),
         bgc_avg = mean(bgc),
         bgc_se = sqrt(var(bgc) / n()),
         cwd_avg = mean(cwd),
         cwd_se = sqrt(var(cwd) / n()),
         soc_avg = mean(soc),
         soc_se = sqrt(var(soc) / n())) %>%
  select(site, agc_avg, agc_se, bgc_avg, bgc_se, 
         cwd_avg, cwd_se, soc_avg, soc_se) %>%
  distinct
  
c_summary %>%
  mutate(total = agc_avg + bgc_avg + cwd_avg + soc_avg,
         total_se = sqrt(sum(agc_se^2 + bgc_se^2 + cwd_se^2 + soc_se^2))) %>%
  View



#------------------------------------------------------------------------------

#---------------------------#
# Summarize soil properties #
#---------------------------#

# Summarize by depth interval
# Can insert the aquaculture data frame or the soil data frame

soil %>%
  select(site, plot, interval, bulk_density, percent_c, c_dens, avg_depth) %>%
  group_by(site, plot) %>%
  mutate(plot_depth = mean(avg_depth),
         plot_bd = mean(bulk_density),
         plot_poc = mean(percent_c),
         plot_c_dens = mean(c_dens),
         plot_depth_se = sqrt(var(avg_depth) / n()),
         plot_bd_se = sqrt(var(bulk_density) / n()),
         plot_poc_se = sqrt(var(percent_c) / n()),
         plot_c_dens_se = sqrt(var(c_dens) / n())) %>%
  arrange(site, plot) %>%
  select(site, plot, plot_bd, plot_bd_se,
         plot_poc, plot_poc_se,
         plot_c_dens, plot_c_dens_se,
         plot_depth, plot_depth_se) %>%
  distinct

