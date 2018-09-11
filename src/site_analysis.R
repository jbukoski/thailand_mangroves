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
trees <- read_excel(paste0(in_dir, "thailand_data.xlsx"), sheet="Trees", col_names = T)
saps <- read_excel(paste0(in_dir, "thailand_data.xlsx"), sheet="Saplings", col_names = T)
seedlings <- read_excel(paste0(in_dir, "thailand_data.xlsx"), sheet="Seedlings", col_names = T)
cwd <- read_excel(paste0(in_dir, "thailand_data.xlsx"), sheet="CWD", col_names = T)
soil <- read_excel(paste0(in_dir, "thailand_data.xlsx"), sheet="Soil", col_names = T)

#-----------------------------------------
# Specify necessary parameters

plot_size <- 5*(7^2)*pi
subplot_size <- 5*(2^2)*pi
sampled_area <- plot_size*7 
sub_sampled_area <- subplot_size*7

site_areas <- tibble(site = c("Krabi", "Nakorn"),
                     area = c(100000, 100000))

#------------------------------------------------------------------------------
# Clean up the dataset

colnames(trees) <- tolower(gsub("[. ]", "_", colnames(trees)))

trees <- trees %>% 
  id_taxon(trees$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2)),
         basal_area = pi * (dbh_cm/100/2)^2) %>%
  dplyr::select(-genus, -species)

#------------------------------------------------------------------------------
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
  mutate(biomass = adj_agb + bgb) %>%
  dplyr::select(-ag_form, -bg_form, -ag_ref, -bg_ref, -params, -top_diam, -stump.vol)

#-------------------------------------------------------------------------------
# Compute sapling biomass

colnames(saps) <- tolower(gsub("[. ]", "_", colnames(saps)))

saps <- saps %>% 
  id_taxon(saps$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2))) %>%
  dplyr::select(-genus, -species)

saps <- saps %>%
  left_join(allom_lookup, by = c("sps_code" = "sps_code")) %>%
  mutate(params = map2(dbh_cm, density, list)) %>%
  mutate(agb = invoke_map_dbl(ag_form, params)) %>%
  mutate(bgb = invoke_map_dbl(bg_form, params))

# Adjust by status, no status == 3 conditions included for saplings

saps <- saps %>%
  mutate(adj_agb = ifelse(is.na(status), agb,
                          ifelse(status == 1, 0.95 * agb, 0.8 * agb)))

# Compute total biomass as function of bgb & adj_agb

saps <- saps %>%
  mutate(biomass = adj_agb + bgb) %>%
  dplyr::select(-ag_form, -bg_form, -ag_ref, -bg_ref, -params)

#-------------------------------------------------------------------------------
# Calculate biomass for the coarse-woody debris pool based on default mean 
# diameters and densities given in the Kauffman & Donato 2012

colnames(cwd) <- tolower(gsub("[. ]", "_", colnames(cwd)))

# Define the mean specific gravities (g/cm^3) of the wood classes; taken from K&D, 2012

cwd_params <- tibble(size = c("fine", "small", "medium", "large"),
                     density = c(0.48, 0.64, 0.71, 0.69),
                     avg_diam = c(0.43, 1.47, 4.52, NA))

#Calculating the mass and volume per plot for each of the size classes, beginning with fine CWD

new_cwd <- cwd %>%
  dplyr::select(-remarks) %>%
  gather(size, n, -site, -date, -recorder, -plot, 
         -subplot, -transect, -data_checked_by, -data_entered_by) %>%
  separate(size, c("size", "status")) %>%
  left_join(cwd_params, by = "size") %>%
  mutate(trnsct_lngth = ifelse(size == "fine", 2,
                               ifelse(size == "small", 3,
                                      ifelse(size == "medium", 5, 12)))) %>%
  mutate(volume = ifelse(size != "large", (pi^2) * ((n * (avg_diam^2)) / (8 * trnsct_lngth)),
                         (pi^2) * (n^2) / (8 * trnsct_lngth)),
         mass = volume * density,
         adj_mass = ifelse(is.na(status), mass,
                           ifelse(status == "rotten", mass * 0.5, mass)))

# Create a summary table

summary_cwd <- new_cwd %>%
  group_by(site, plot, size, status) %>%
  summarise(total = mean(n), mass = mean(mass)) %>%
  group_by(plot) %>%
  mutate(plot_mass = sum(mass))

site_cwd <- summary_cwd %>%
  group_by(site) %>%
  summarize(cwd = mean(plot_mass),
            cwd_se = mean(sd(unique(plot_mass) / sqrt(n_distinct(plot_mass)))))

#-------------------------------------------------------------------------------
#Obtain estimate of biomass per hectare based on all plots

trees_ci <- trees %>%
  left_join(site_areas, by = "site") %>%
  mutate(N = area) %>%
  group_by(site) %>%
  dplyr::summarise(t_val = qt(0.975, sampled_area-1),
                   tau_cap = mean(N) / sampled_area * sum(biomass),
                   v_cap = mean(N)^2 * (1/sampled_area - 1/mean(N)) * var(biomass),
                   se_tau_cap = sqrt(v_cap),
                   mean_biomass = tau_cap / (mean(N)/10000) / 1000,
                   lower = (tau_cap - (se_tau_cap * t_val)) / (mean(N)/10000)/1000,
                   upper = (tau_cap + (se_tau_cap * t_val)) / (mean(N)/10000)/1000,
                   moe = 100*((se_tau_cap * t_val)/(mean(N) / 10000))/mean_biomass/1000,
                   ag_tau = mean(N)/sampled_area * sum(agb),
                   ag_var = mean(N)^2 * (1/sampled_area - 1/mean(N)) * var(agb),
                   ag_se = sqrt(ag_var) / (mean(N)/10000) / 1000,
                   ag_mean = ag_tau / (mean(N)/10000) / 1000,
                   bg_tau = mean(N)/sampled_area * sum(bgb),
                   bg_var = mean(N)^2 * (1/sampled_area - 1/mean(N)) * var(bgb),
                   bg_se = sqrt(bg_var) / (mean(N)/10000) / 1000,
                   bg_mean = bg_tau / (mean(N)/10000) / 1000
  )

#-------------------------------------------------------------------------------
#Obtain estimate of biomass per hectare based on all plots

saps_ci <- saps %>%
  left_join(site_areas, by = "site") %>%
  mutate(N = area) %>%
  group_by(site) %>%
  dplyr::summarise(tau_cap = mean(N)/sampled_area * sum(biomass),
                   mean_biomass = tau_cap / (mean(N)/10000) / 1000,
                   v_cap = mean(N)^2 * (1/sampled_area - 1/mean(N)) * var(biomass),
                   se_tau_cap = sqrt(v_cap),
                   t_val = qt(0.975, sampled_area-1),
                   lower = (tau_cap - (se_tau_cap * t_val))/(mean(N)/10000)/1000,
                   upper = (tau_cap + (se_tau_cap * t_val))/(mean(N)/10000)/1000,
                   moe = 100*((se_tau_cap * t_val)/(mean(N) / 10000))/mean_biomass/1000
  )

#-------------------------------------------------------------------------------
# Classify each plot by species
# Species metrics at the plot level

plot_sps <- trees %>%
  group_by(site) %>%
  dplyr::mutate(plot_n = n(), 
                plot_ba = sum(basal_area)) %>%
  group_by(site, genus, species) %>%
  dplyr::summarize(sps_ba = sum(basal_area),
                   sps_n = n(),
                   rel_sps_ba = mean(sps_ba)/mean(plot_ba) * 100,
                   rel_sps_n = mean(sps_n)/mean(plot_n) * 100,
                   rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  dplyr::select(site, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(-rel_imp)

plot_sps %>%
  dplyr::select(-species) %>%
  group_by(site, genus) %>%
  filter(rel_imp > 1) %>%
  summarise(a_rel_sps_n = sum(rel_sps_n),
            b_rel_sps_ba = sum(rel_sps_ba),
            c_rel_imp = sum(rel_imp)) %>%
  gather(key, value, -site, -genus) %>%
  mutate(key = factor(key, labels = c("Rel. Dens", "Rel. BA", "Rel. Imp."))) %>%
  arrange(site, genus) %>%
  ggplot(aes(x = key, y = value, fill = genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_bar(stat = "identity", position = "stack", color="black", 
           show.legend=FALSE, width = 0.7) +
  facet_grid(site ~ .) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  labs(y = NULL,
       x = NULL) +
  scale_fill_discrete(name="Genera", 
                      labels = c("Avicennia", 
                                 "Bruguiera",
                                 "Excoecaria",
                                 "Rhizophora",
                                 "Xylocarpus"))

# Calculate for subplots

subplot_sps <- trees %>%
  group_by(site, plot, subplot) %>%
  dplyr::mutate(agg_n = n(),
                agg_ba = sum(basal_area)) %>%
  group_by(site, plot, subplot, genus, species) %>%
  dplyr::summarize(sps_ba = sum(basal_area),
                   sps_n = n(),
                   rel_sps_ba = mean(sps_ba)/mean(agg_ba)*100,
                   rel_sps_n = mean(sps_n)/mean(agg_n)*100,
                   rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  dplyr::select(site, plot, subplot, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(plot, subplot, -rel_imp)

#----------------------------------------------------------------------------------------------
#Summarize forest structure

plot_structure <- trees %>%
  group_by(site, plot) %>%
  dplyr::summarise(n = n(), 
                   ba = sum(basal_area),
                   dbh = mean(dbh_cm),
                   min_dbh = min(dbh_cm),
                   max_dbh = max(dbh_cm)) %>%
  mutate(plot_size = plot_size,
         ba_per_ha = 10000*ba/plot_size,
         dbh_per_ha = dbh,
         n_per_ha = 10000*n/plot_size)

site_structure <- plot_structure %>%
  group_by(site) %>%
  dplyr::summarise(
    n = mean(n_per_ha),
    n_se = sd(n_per_ha) / sqrt(n_distinct(n)),
    ba = mean(ba_per_ha),
    ba_se = sd(ba_per_ha) / sqrt(n_distinct(ba_per_ha)),
    dbh = mean(dbh_per_ha),
    dbh_se = sd(dbh_per_ha) / sqrt(n_distinct(dbh_per_ha))
  )
  
trees %>%
  group_by(site) %>%
  ggplot() +
  geom_histogram(aes(x = dbh_cm, fill = site), stat = "bin", binwidth = 1, colour = "black") +
  facet_grid(site ~ .) +
  theme_bw() +
  labs(x = "DBH (cm)", y = "Frequency") +
  theme(text = element_text(size = 22)) +
  scale_fill_discrete(name = "Site")

plot_structure %>%
  select(site, plot, ba_per_ha, dbh_per_ha, n_per_ha) %>%
  mutate(n_per_ha = n_per_ha/100) %>%
  gather(key, value, -site, -plot) %>%
  mutate(key = factor(key, labels = c("BA (sq.m/ha)", "Avg. DBH (cm)", "Trees/ha (100s)"))) %>%
  arrange(site, key) %>%
  ggplot(aes(x = site, group = site, y = value, fill = site)) +
  stat_boxplot(geom='errorbar', width = 1) +
  geom_boxplot() +
  facet_grid(key ~ .) +
  theme_bw() +
  coord_flip() +
  scale_fill_discrete(name = NULL) +
  theme(legend.position="bottom") +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 24))
  
#-----------------------------------------------------------------------------------------------
#Soil analysis

names(soil) <- tolower(names(soil))
names(soil) <- gsub('[.]', '_', names(soil))

soil <- soil %>%
  mutate(c_dens = bulk_density * (percent_c/100),
         int_volume = ifelse(interval == 5, 
                             ((avg_depth/100) - (100/100)) * 10000, 
                             ((int_b/100) - (int_a/100)) * 10000),
         soc_per_ha = int_volume * c_dens)

soil_summary <- soil %>%
  group_by(site, plot, subplot) %>%
  dplyr::mutate(plot_c = sum(soc_per_ha)) %>%
  group_by(site, plot, subplot) %>%
  dplyr::summarize(soc = mean(plot_c))

soil_plot <- soil %>%
  group_by(site, plot, subplot) %>%
  dplyr::mutate(subplot_c = sum(soc_per_ha)) %>%
  group_by(site, plot) %>%
  summarise(plot_soc = mean(subplot_c),
            min_soc = min(subplot_c),
            max_soc = max(subplot_c))
            

soil_site <- soil_plot %>%
  group_by(site, plot) %>%
  dplyr::mutate(se_plot_c = sd(unique(plot_soc)) / sqrt(n_distinct(plot_soc))) %>%
  group_by(site) %>%
  dplyr::summarize(soc = mean(plot_soc),
                   soc_se = mean(unique(se_plot_c)))

#-----------------------------------------------------------------------------
# Generate plots

soil_plot %>%
  filter(plot <= 7) %>%
  ggplot(aes(x = plot, y = plot_soc, col = site)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = min_soc, 
                    ymax = max_soc), width = 0.15,
                size = 1.1) + 
  labs(y = "Soil Organic Carbon (Mg C/ha)",
       x = "Plot") +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  ylim(0, 1500) +
  scale_colour_discrete(name="Site", 
                      labels = c("Krabi", "Nakorn"))


soil_site %>% 
  ggplot(aes(x = site, y = soc)) + 
  geom_errorbar(aes(ymin = soc - soc_se, 
                    ymax = soc + soc_se), width = 0.05) + 
  geom_point() + 
  ggtitle("Soil Carbon Plot Means +/- 1 St. Error") +
  ylab("Soil Organic Carbon (Mg/ha)") +
  xlab("Site") +
  theme_bw() +
  ylim(0, 1000)

#------------------------------------------------------------------------------
# Add other relevant parameters to characterize the plot (i.e. mean DBH and Tons ag.biomass)

error_summary <- trees_ci %>%
  left_join(soil_site, by = c("site" = "site")) %>%
  left_join(site_cwd, by = c("site" = "site")) %>%
  dplyr::select(site, ag_se, bg_se, soc_se, cwd_se) %>%
  mutate(ag_se = ag_se,
         bg_se = bg_se,
         soc_se = soc_se,
         cwd_se = cwd_se) %>%
  gather(pool, error, -site) %>%
  arrange(site) %>%
  mutate(pool = rep(c("agc", "bgc", "soc", "cwd"), 2))

summary <- trees_ci %>%
  dplyr::select(site, mean_biomass, ag_mean, bg_mean) %>%
  left_join(soil_site, by = c("site" = "site")) %>%
  left_join(site_cwd, by = c("site" = "site")) %>%
  mutate(bgc = bg_mean *.46 * -1,
         agc = ag_mean *.46,
         soc = soc * -1,
         cwd = cwd * 0.5) %>%
  dplyr::select(site, agc, bgc, soc, cwd) %>%
  gather(pool, value, -site) %>%
  left_join(error_summary, by = c("site", "pool")) %>%
  mutate(error2 = error)

##---------------------------------------------------------------------------------
## Visualizing the data

# Ordering of stacked bars is determined by levels of the factor

test_dat <- tibble(
  site = c(rep("Krabi", 4), rep("Nakorn", 4)),
  pool = rep(c(rep("Above", 2), rep("Below", 2)), 2),
  min_err = c(72.3, NA, -888.9 - 53.41, NA, 108 + 8.10, NA, -39.9 - 251 - 18 - 8.17, NA),
  max_err = c(72.3  + 18.88, NA, -888.9, NA, 108 + 8.10 + 23.2 + 2.85, NA, -39.9 - 251, NA)
)

p1 <- summary %>% 
  mutate(pool = factor(pool, levels = c("agc", "cwd", "soc", "bgc"))) %>%
  ggplot(aes(x = site, y = value, fill = pool)) +
  geom_bar(stat = "identity", width = 0.3) +
  geom_bar(stat = "identity", width = 0.3, color="black", show.legend=FALSE)+
  theme_bw() +
  xlab("Site") +
  ylab("Carbon storage (Mg C/ha)") +
  ylim(-1000, 200) +
  scale_fill_discrete(name="Ecosystem C Pools", 
                      breaks = c("agc", "cwd", "bgc", "soc"),
                      labels = c("Aboveground Biomass", 
                                 "Coarse Woody Debris",
                                 "Belowground Biomass",
                                 "Soil Organic Carbon")) +
  theme(text = element_text(size = 22))

p1 +  
  geom_linerange(aes(x = test_dat$site,
                     ymin = test_dat$min_err, 
                     ymax = test_dat$max_err)) +
  geom_segment(aes(x = 0.95, xend = 1.05, y = -888.9 - 53.41, yend = -888.9 - 53.41)) +
  geom_segment(aes(x = 0.95, xend = 1.05, y = 72.3  + 18.88, yend = 72.3  + 18.88)) +
  geom_segment(aes(x = 1.95, xend = 2.05, 
                   y = -39.9 - 251 - 18 - 8.17, 
                   yend = -39.9 - 251 - 18 - 8.17)) +
  geom_segment(aes(x = 1.95, xend = 2.05, 
                   y = 108 + 8.10 + 23.2 + 2.85, 
                   yend = 108 + 8.10 + 23.2 + 2.85))

#-----------------------------------------------------------------------------------
# Build a summary table for visualizations of the data at the subplot level 

sp_summary <- trees %>%
  group_by(site, plot, subplot) %>%
  left_join(site_areas, by = "site") %>%
  mutate(agb_tot = area*(sum(agb)/(pi*(7^2))),
         agb_ha = (10*agb_tot)/area,
         bgb_tot = area*(sum(bgb)/(pi*(7^2))),
         bgb_ha = (10*bgb_tot)/area,
         n_sps = nlevels(factor(sps_code))) %>%
  left_join(meta, by = c("site" = "Site", "plot" = "Plot", "subplot" = "Subplot")) %>%
  left_join(soil_summary, by = c("site", "plot", "subplot")) %>%
  mutate(tot_c = soc + agb_ha + bgb_ha,
         distance = `Distance from shoreline`) %>%
  dplyr::select(site, plot, subplot, agb_tot, agb_ha, bgb_tot, bgb_ha, soc, tot_c, distance, n_sps) %>%
  distinct()

trees %>%
  filter(plot <= 7) %>%
  group_by(site, plot) %>%
  summarize(n = n(),
            agb_se = sd(agb)/sqrt(n) * 0.43,
            agb = 10000*sum(agb)/plot_size/1000 * 0.43,
            bgb_se = sd(bgb)/sqrt(n) * 0.43,
            bgb = 10000*sum(bgb)/plot_size/1000 * 0.43) %>%
  ggplot(aes(x = plot, y = agb + bgb, col = site)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = agb + bgb - agb_se - bgb_se, 
                    ymax = agb + bgb + agb_se + bgb_se), width = 0.15,
                size = 1.1) + 
  labs(y = "Biomass Carbon (Mg C/ha)",
       x = "Plot") +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  scale_colour_discrete(name="Site", 
                        labels = c("Krabi", "Nakorn"))


            