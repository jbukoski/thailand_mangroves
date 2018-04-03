# This script analyzes the sampling data of mangrove carbon stocks in Southeast Asia

cat("Site analysis          ", format(Sys.time(), "%A, %d %B %Y"),
    "\nMangrove sampling: May-August 2015          ", "Jacob J. Bukoski",
    "\n================================================================\n")

rm(list=ls())

setwd("~/Dropbox/mangrove-work/data")

library("ggplot2")
library("lmfor")
library("tidyverse")
library("readxl")

#------------------------------------------------------------------------------
## Load in data and set necessary site level values

meta <- read_excel("thailand_data.xlsx", sheet="Metadata", col_names = T)
trees <- read_excel("thailand_data.xlsx", sheet="Trees", col_names = T)
saps <- read_excel("thailand_data.xlsx", sheet="Saplings", col_names = T)
seedlings <- read_excel("thailand_data.xlsx", sheet="Seedlings", col_names = T)
cwd <- read_excel("thailand_data.xlsx", sheet="CWD", col_names = T)
soil <- read_excel("thailand_data.xlsx", sheet="Soil", col_names = T)

plot_size <- 5*(7^2)*pi
subplot_size <- 5*(2^2)*pi
sampled_area <- plot_size*7 
sub_sampled_area <- subplot_size*7
krabi_ha <- 102120000
nakorn_ha <- 56800000

site_areas <- tibble(site = c("Krabi", "Nakorn"),
                     area = c(102120000, 56800000))

#------------------------------------------------
# Source helper functions and allometry table

source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/helper_funcs.R")
source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/allometry.R")
source("/home/jbukoski/Dropbox/mangrove-work/model-files/Function-SummarySE.R")

#------------------------------------------------------------------------------
# Clean up the dataset.

colnames(trees) <- tolower(gsub("[.]", "_", colnames(trees)))

trees <- trees %>% 
  id_taxon(trees$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2)),
         basal_area = 0.00007854 * dbh_cm^2) %>%
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

# Adjust ag.biomass variable based on Status variable
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
                          ifelse(status == 1, 0.95*agb, 0.8*agb)))

# Compute total biomass as function of bgb & adj_agb

saps <- saps %>%
  mutate(biomass = adj_agb + bgb) %>%
  dplyr::select(-ag_form, - bg_form, - ag_ref, 
                - bg_ref, -params)

#-------------------------------------------------------------------------------
#Calculate biomass for the coarse-woody debris pool based on default mean diameters
#and densities given in the Kauffman and Donato protocols

colnames(cwd) <- tolower(gsub("[. ]", "_", colnames(cwd)))

#define the mean specific gravities (g/cm^3) of the wood classes; taken from K&D, 2012

cwd_params <- tibble(size = c("fine", "small", "medium", "large"),
                     density = c(0.48, 0.64, 0.71, 0.69),
                     avg_diam = c(0.43, 1.47, 4.52, NA))
  
#Calculating the mass and volume per plot for each of the size classes, beginning with fine CWD
  
new_cwd <- cwd %>%
  dplyr::select(-remarks) %>%
  gather(size, n, -site, -date, -recorder, -plot, -subplot, -transect, -data_checked_by, -data_entered_by) %>%
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
  dplyr::group_by(site, plot, size, status) %>%
  dplyr::summarise(total = mean(n), mass = mean(mass)) %>%
  dplyr::group_by(plot) %>%
  dplyr::mutate(plot_mass = sum(mass))

site_cwd <- summary_cwd %>%
  group_by(site) %>%
  dplyr::summarize(cwd = mean(plot_mass),
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
                   lower = (tau_cap - (se_tau_cap * t_val))/(mean(N)/10000)/1000,
                   upper = (tau_cap + (se_tau_cap * t_val))/(mean(N)/10000)/1000,
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
  group_by(site, plot) %>%
  dplyr::mutate(plot_n = n(), 
                plot_ba = sum(basal_area)) %>%
  group_by(site, plot, genus, species) %>%
  dplyr::summarize(sps_ba = sum(basal_area),
                   sps_n = n(),
                   rel_sps_ba = mean(sps_ba)/mean(plot_ba) * 100,
                   rel_sps_n = mean(sps_n)/mean(plot_n) * 100,
                   rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  select(site, plot, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(plot, -rel_imp)

ggplot(plot_sps, aes(species, rel_imp)) +
  geom_boxplot()
  
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
  select(site, plot, subplot, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(plot, subplot, -rel_imp)

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

soil_site <- soil %>%
  group_by(site, plot, subplot) %>%
  dplyr::mutate(plot_c = sum(soc_per_ha)) %>%
  group_by(site, plot) %>%
  dplyr::mutate(se_plot_c = sd(unique(plot_c)) / sqrt(n_distinct(plot_c))) %>%
  group_by(site) %>%
  dplyr::summarize(soil_c = mean(plot_c),
                   soc_se = mean(unique(se_plot_c)))

soil_summary %>% 
  ggplot(aes(x = plot, y = mean_plot_c)) + 
  geom_errorbar(aes(ymin = mean_plot_c - se_plot_c, 
                    ymax = mean_plot_c + se_plot_c), width = 0.1) + 
  geom_point() + 
  ggtitle("Soil Carbon Plot Means +/- 1 St. Error")

#------------------------------------------------------------------------------
# Add other relevant parameters to characterize the plot (i.e. mean DBH and Tons ag.biomass)

error_summary <- trees_ci %>%
  left_join(soil_site, by = c("site" = "site")) %>%
  left_join(site_cwd, by = c("site" = "site")) %>%
  select(site, ag_se, bg_se, soc_se, cwd_se) %>%
  mutate(ag_se = ag_se,
         bg_se = bg_se,
         soc_se = soc_se,
         cwd_se = cwd_se) %>%
  gather(pool, error, -site) %>%
  arrange(site) %>%
  mutate(pool = rep(c("agc", "bgc", "soc", "cwd"), 2))

summary <- trees_ci %>%
  select(site, mean_biomass, ag_mean, bg_mean) %>%
  left_join(soil_site, by = c("site" = "site")) %>%
  left_join(site_cwd, by = c("site" = "site")) %>%
  mutate(bgc = bg_mean *.46 * -1,
         agc = ag_mean *.46,
         soc = soil_c * -1,
         cwd = cwd * 0.5) %>%
  select(site, agc, bgc, soc, cwd) %>%
  gather(pool, value, -site) %>%
  left_join(error_summary, by = c("site", "pool")) %>%
  mutate(error2 = error)
  
##---------------------------------------------------------------------------------
## Visualizing the data

# Ordering of stacked bars is determined by levels of the factor

summary %>%
  mutate(pool = factor(pool, levels = c("agc", "cwd", "soc", "bgc"))) %>%
  ggplot(aes(x = site, y = value, fill = pool)) +
  geom_bar(stat = "identity", width = 0.3) +
  geom_bar(stat = "identity", width = 0.3, color="black", show.legend=FALSE)+
  theme_bw() +
  xlab("Site") +
  ylab("Carbon storage (Mg C / ha)") +
  ylim(-1000, 200)

ggplot(summary) +
  geom_bar(aes(x=site, y = value), stat = "identity") +
  geom_errorbar(aes(ymin = error, ymax = error2), stat = "identity")

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
  select(site, plot, subplot, agb_tot, agb_ha, bgb_tot, bgb_ha, soc, tot_c, distance, n_sps) %>%
  distinct()
  
sp_summary %>%
  filter(site == "Nakorn") %>%
  ggplot(aes(distance, soc)) +
  geom_point(aes(col = site)) +
  facet_grid(. ~ site)

lm_dat <- sp_summary %>%
  filter(site == "Nakorn", !is.na(soc))
  
lm_soc <- lm(soc ~ distance, data = lm_dat)
