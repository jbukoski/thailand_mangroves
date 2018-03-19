# This script analyzes the sampling data of mangrove carbon stocks in Southeast Asia

cat("Site analysis          ", format(Sys.time(), "%A, %d %B %Y"),
    "\nMangrove sampling: May-August 2015          ", "Jacob J. Bukoski",
    "\n================================================================\n")

rm(list=ls())
rm(list=ls()[!(ls() %in% c("read", "total.summary", "foo.mix", "foo.nakorn", "test"))])
##rm(list=setdiff(ls(), read))  #Mechanism by which to clear list except select items

setwd("~/Dropbox/mangrove-work/data")
library("xlsx")
library("ggplot2")
library("lmfor")
library("tidyverse")

#------------------------------------------------------------------------------
#Interactive prompt to import the data for each site; if this works, can adjust later 
#to fill in parameters such as plot.size (n), total hectares (N), sampled.area, etc.

#NOTE: NEED TO ADJUST TOTAL HA FOR VIETNAMESE SITES - PLACE FILLED WITH 1,000,000 sq m

read <- function() {
  ANSWER = readline("Which site would you like to analyze? \n
                    Select one of: krabi, nakorn, nam_dinh, hai_phong, quang_ninh: ")
  if(ANSWER == "krabi")
  {cat("Okay, loading the Krabi (Thailand) data...\n")
    all.trees <<- read.xlsx("krabi_data.xlsx", sheetName="Trees", header=T);
    saplings <<- read.xlsx("krabi_data.xlsx", sheetName="Saplings", header=T);
    seedlings <<- read.xlsx("krabi_data.xlsx", sheetName="Seedlings", header=T);
    soil <<- read.xlsx("krabi_data.xlsx", sheetName="Soil", header=T);
    cwd <<- read.xlsx("krabi_data.xlsx", sheetName="CWD", header=T)
    plot.size<<-5*(7^2)*pi; subplot.size<<-5*(2^2)*pi;sampled.area<<-plot.size*7; 
    sub.sampled.area<<-subplot.size*7; N<<-102120000; site<<-"Krabi";
    cat("...loaded!")} else
      if (ANSWER == "nakorn")
      {cat("Okay, loading the Nakorn Si Thammarat (Thailand) data...\n");
        all.trees <<- read.xlsx("nakorn_data.xlsx", sheetName="Trees", header=T);
        saplings <<- read.xlsx("nakorn_data.xlsx", sheetName="Saplings", header=T);
        seedlings <<- read.xlsx("nakorn_data.xlsx", sheetName="Seedlings", header=T);
        soil <<- read.xlsx("nakorn_data.xlsx", sheetName="Soil", header=T);
        cwd <<- read.xlsx("nakorn_data.xlsx", sheetName="CWD", header=T)
        plot.size<<-5*(7^2)*pi; subplot.size<<-5*(2^2)*pi; sampled.area<<-plot.size*10;
        sub.sampled.area<<-subplot.size*10; N<<-56800000; site<<-"Nakorn"
        cat("...loaded!")} else
          if (ANSWER == "nam_dinh")
          {cat("Okay, loading the Nam Dinh (Vietnam) data...\n");
            all.trees <<- read.xlsx("nam_dinh_data.xlsx", sheetName="Trees", header=T);
            seedlings <<- read.xlsx("nam_dinh_data.xlsx", sheetName="Seedlings", header=T);
            soil <<- read.xlsx("nam_dinh_data.xlsx", sheetName="Soil", header=T);
            cwd <<- read.xlsx("nam_dinh_data.xlsx", sheetName="CWD", header=T)
            plot.size<<-18; sampled.area<<-plot.size*4; N<<-23530000; site<<-"Nam Dinh"
            cat("...loaded!")} else
              if (ANSWER == "hai_phong")
              {cat("Okay, loading the Hai Phong (Vietnam) data...\n");
                all.trees <<- read.xlsx("hai_phong_data.xlsx", sheetName="Trees", header=T);
                saplings <<- read.xlsx("hai_phong_data.xlsx", sheetName="Saplings", header=T);
                seedlings <<- read.xlsx("hai_phong_data.xlsx", sheetName="Seedlings", header=T);
                soil <<- read.xlsx("hai_phong_data.xlsx", sheetName="Soil", header=T);
                cwd <<- read.xlsx("hai_phong_data.xlsx", sheetName="CWD", header=T)
                plot.size<<-6*(7^2)*pi; subplot.size<<-6*(2^2)*pi; 
                sampled.area<<-plot.size*3; sub.sampled.area<<-subplot.size*3; N<<-24330000; site<<-"Hai Phong" 
                cat("...loaded!")} else
                  if (ANSWER == "quang_ninh")
                  {cat("Okay, loading the Quang Ninh (Vietnam) data...\n");
                    all.trees <<- read.xlsx("quang_ninh_data.xlsx", sheetName="Trees", header=T);
                    saplings <<- read.xlsx("quang_ninh_data.xlsx", sheetName="Saplings", header=T);
                    seedlings <<- read.xlsx("quang_ninh_data.xlsx", sheetName="Seedlings", header=T);
                    soil <<- read.xlsx("quang_ninh_data.xlsx", sheetName="Soil", header=T)
                    cwd <<- read.xlsx("quang_ninh_data.xlsx", sheetName="CWD", header=T)
                    plot.size<<-6*(7^2)*pi; subplot.size<<-6*(2^2)*pi; 
                    sampled.area<<-plot.size*3; sub.sampled.area<<-subplot.size*3; N<<-162680000; site<<-"Quang Ninh"
                    cat("...loaded!")}
}

read()

#------------------------------------------------
# Source helper functions and allometry table

source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/helper_funcs.R")
source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/allometry.R")

#------------------------------------------------------------------------------
# Clean up the dataset.

trees <- all.trees
colnames(trees) <- tolower(colnames(trees))

trees <- trees %>% 
  id_taxon(trees$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2)),
         basal_area = 0.00007854*dbh.cm^2) %>%
  dplyr::select(-genus, -species, everything())

#------------------------------------------------------------------------------
# Calculate above-ground biomass using species-specific allometric equations. 
# Where species-specific equations are not available, used Komiyama et al 2005 
# general equation with species specific wood densities

trees <- trees %>%
  left_join(allom_lookup, by = c("sps_code" = "sps_code")) %>%
  mutate(params = map2(dbh.cm, density, list)) %>%
  mutate(agb = invoke_map_dbl(ag_form, params)) %>%
  mutate(bgb = invoke_map_dbl(bg_form, params))

# Adjust ag.biomass variable based on Status variable
# Calculate cone if base.cm measurement exists, otherwise assume a cylinder

trees <- trees %>%
  mutate(top.diam = ifelse(status == 3 & !is.na(base.cm) & height.m < 1.37, 
                           base.cm-(100*height.m*((base.cm-dbh.cm)/(100*height.m))),
                           base.cm-(100*height.m*((base.cm-dbh.cm)/137)))) %>%
  mutate(stump.vol = ifelse(status == 3 & !is.na(base.cm) & height.m < 1.37, (pi*100*height.m)/12*(base.cm^2+top.diam^2+(base.cm*top.diam)),
                            ifelse(status == 3 & !is.na(base.cm) & height.m >= 1.37, (pi*100*height.m)/12*(base.cm^2+top.diam^2+(base.cm*top.diam)), NA))) %>%
  mutate(adj_agb = ifelse(is.na(status), agb,
                          ifelse(status == 1, 0.95*agb, 
                                 ifelse(status == 2, 0.8*agb, 
                                        ifelse(status == 3 & !is.na(base.cm), density*stump.vol/1000, density*pi*dbh.cm*height.m*100/1000)))))

# Compute total biomass

trees <- trees %>%
  mutate(biomass = adj_agb + bgb) %>%
  dplyr::select(-ag_form, - bg_form, - ag_ref, 
                - bg_ref, -params, - top.diam, - stump.vol)
  
#-------------------------------------------------------------------------------
# Compute sapling biomass

saps = saplings

colnames(saps) <- tolower(colnames(saps))

saps <- saps %>% 
  id_taxon(saps$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2))) %>%
  dplyr::select(-genus, -species)

saps <- saps %>%
  left_join(allom_lookup, by = c("sps_code" = "sps_code")) %>%
  mutate(params = map2(dbh.cm, density, list)) %>%
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

head(saps)
  
#-------------------------------------------------------------------------------
#Calculate biomass for the coarse-woody debris pool based on default mean diameters
#and densities given in the Kauffman and Donato protocols

colnames(cwd) <- tolower(colnames(cwd))

#define the mean specific gravities (g/cm^3) of the wood classes; taken from K&D, 2012

cwd_params <- tibble(size = c("fine", "small", "medium", "large"),
                     density = c(0.48, 0.64, 0.71, 0.69),
                     avg_diam = c(0.43, 1.47, 4.52, NA))
  
#Calculating the mass and volume per plot for each of the size classes, beginning with fine CWD
  
new_cwd <- cwd %>%
  dplyr::select(-remarks) %>%
  gather(size, n, -site, -date, -recorder, -data.checked.by, -data.entered.by, 
         -plot, -subplot, -transect) %>%
  separate(size, c("size", "status")) %>%
  left_join(cwd_params, by = "size") %>%
  mutate(trnsct_lngth = ifelse(size == "fine", 2,
                               ifelse(size == "small", 3,
                                      ifelse(size == "medium", 5, 12)))) %>%
  mutate(volume = ifelse(size != "large", (pi^2) * ((n*(avg_diam^2))/(8*trnsct_lngth)),
                         (pi^2) * (n^2) / (8*trnsct_lngth)),
         mass = volume * density,
         adj_mass = ifelse(is.na(status), mass,
                           ifelse(status == "rotten", mass * 0.5, mass)))

# Create a summary table

summary_cwd <- new_cwd %>%
  dplyr::group_by(plot, size, status) %>%
  dplyr::summarise(total = mean(n), mass = mean(mass)) %>%
  dplyr::group_by(plot) %>%
  dplyr::mutate(plot_mass = sum(mass))

head(summary_cwd)
  
#-------------------------------------------------------------------------------
#Obtain estimate of biomass per hectare based on all plots

trees_ci <- trees %>%
  group_by(site) %>%
  dplyr::summarise(tau_cap = N/sampled.area * sum(biomass),
                   mean_biomass = tau_cap / (N/10000) / 1000,
                   v_cap = N^2 * (1/sampled.area - 1/N) * var(biomass),
                   se_tau_cap = sqrt(v_cap),
                   t_val = qt(0.975, sampled.area-1),
                   lower = (tau_cap - (se_tau_cap * t_val))/(N/10000)/1000,
                   upper = (tau_cap + (se_tau_cap * t_val))/(N/10000)/1000,
                   moe = 100*((se_tau_cap * t_val)/(N / 10000))/mean_biomass/1000
                   )

trees_ci
#-------------------------------------------------------------------------------
#Obtain estimate of biomass per hectare based on all plots

saps_ci <- saps %>%
  group_by(site) %>%
  dplyr::summarise(tau_cap = N/sampled.area * sum(biomass),
                   mean_biomass = tau_cap / (N/10000) / 1000,
                   v_cap = N^2 * (1/sampled.area - 1/N) * var(biomass),
                   se_tau_cap = sqrt(v_cap),
                   t_val = qt(0.975, sampled.area-1),
                   lower = (tau_cap - (se_tau_cap * t_val))/(N/10000)/1000,
                   upper = (tau_cap + (se_tau_cap * t_val))/(N/10000)/1000,
                   moe = 100*((se_tau_cap * t_val)/(N / 10000))/mean_biomass/1000
  )

saps_ci

#-------------------------------------------------------------------------------
# Classify each plot by species
# At plot level, incorporate frequency metric
# At subplot level, do not incorporate frequency metric

# Species metrics at the plot level, does not include frequency

plot_sps <- trees %>%
  group_by(plot) %>%
  dplyr::mutate(plot_n = n(), 
                plot_ba = sum(basal_area)) %>%
  group_by(plot, genus, species) %>%
  dplyr::summarize(sps_ba = sum(basal_area),
                   sps_n = n(),
                   rel_sps_ba = mean(sps_ba)/mean(plot_ba) * 100,
                   rel_sps_n = mean(sps_n)/mean(plot_n) * 100,
                   rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  select(plot, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(plot, -rel_imp)

ggplot(plot_sps, aes(species, rel_imp)) +
  geom_boxplot()
  
# Calculate for subplots

subplot_sps <- trees %>%
  group_by(plot, subplot) %>%
  dplyr::mutate(agg_n = n(),
                agg_ba = sum(basal_area)) %>%
  group_by(plot, subplot, genus, species) %>%
  dplyr::summarize(sps_ba = sum(basal_area),
                   sps_n = n(),
                   rel_sps_ba = mean(sps_ba)/mean(agg_ba)*100,
                   rel_sps_n = mean(sps_n)/mean(agg_n)*100,
                   rel_imp = (rel_sps_ba + rel_sps_n) / 2) %>%
  select(plot, subplot, genus, species, rel_sps_ba, rel_sps_n, rel_imp) %>%
  arrange(plot, subplot, -rel_imp)



test <- subplot_sps %>%
  filter(plot == 2, subplot == 1) %>%
  ggplot(aes(rel_imp ~ species)) +
  geom_bar()

#------------------------------------------------------------------------------
#Calculate relative frequency for each species within each plot
#Frequency is defined as presence of a species within a plot

f = aggregate(Species ~ Subplot + Plot, data = trees, 
              FUN=function(x) paste(unique(x), collapse=', '))

n.plot = length(levels(as.factor(trees$Plot)))
n.subplot = length(levels(as.factor(trees$Subplot)))
n.sps = length(levels(as.factor(trees$Species)))

a=matrix(0, nrow(f), n.sps)
colnames(a) = levels(trees$Species)
freq = as.data.frame(cbind(f, a))

counts = freq[,4:ncol(freq)]

Species = as.factor(trees$Species)

freq[-(1:3)] = sapply(colnames(freq[-(1:3)]), grepl, x=freq$Species  ) + 0

colnames(a) = as.character(levels(trees$Species))
sums = a

Plot = levels(as.factor(trees$Plot))
plot.sums = rowsum(freq[,c(4:ncol(freq))], group = freq$Plot)
rel.sums = (100*(plot.sums/as.numeric(n.subplot)))
plot.sums = cbind(Plot,rel.sums)

sps.freqs = sps.stems[,1:2]

Rel.Freq = rep(0,nrow(sps.freqs)); 
sps.freqs = cbind(sps.freqs,Rel.Freq)

for(j in 1:nrow(plot.sums))
  for(k in 1:ncol(plot.sums))
    for(l in 1:nrow(sps.freqs))
      if(plot.sums$Plot[j] == sps.freqs$Plot[l] &&
           colnames(plot.sums[k]) == sps.freqs$Species[l])
      {sps.freqs$Rel.Freq[l] = plot.sums[j,k]}

head(sps.freqs)
#------------------------------------------------------------------------------
#Combine relative basal area, relative stem density, and relative 
#frequency to compute the relative importance value

stems.species = as.character(sps.stems$Species)

rel.imp = sps.stems[,1:2]
rel.imp = cbind(rel.imp, round(sps.BAs$rel.basal.area,1), 
                round(sps.stems$rel.stems,1), round(sps.freqs$Rel.Freq,1))

attach(rel.imp)

importance = (rel.imp[,3]  + rel.imp[,4] + rel.imp[,5])/3

rel.imp = cbind(rel.imp, round(importance,1))

colnames(rel.imp) = c("Species", "Plot", "Rel.BA", 
                      "Rel.Stems", "Rel.Freq", "Rel.Importance")

rel.imp

##NOTE: Figure out how to sort the rel.imp data frame by rel.imp. 
##End goal is to print both first AND second most important species

#------------------------------------------------------------------------------
#Append the species name to the importance value pulled for each of the plots

attach(rel.imp)

imp.values = aggregate(Rel.Importance ~ Plot, data=rel.imp, max)
Species = rep(0,nrow(imp.values))
imp.values = cbind(imp.values,Species)

for(k in 1:nrow(imp.values)) {
    for(j in 1:nrow(rel.imp)) {
  if(imp.values$Rel.Importance[k] == rel.imp$Rel.Importance[j])
  imp.values$Species[k] = as.character(rel.imp$Species[j])}}

#-----------------------------------------------------------------------------------------------
#Soil analysis

attach(soil)
soil.c = as.data.frame(cbind(Plot, Subplot, Interval, Int.A, Int.B,
                             round(Avg.depth, 1), round(Bulk.density,2), Percent.C))
soil.c$Plot = as.numeric(soil.c$Plot)
soil.c$C.density = as.numeric(Bulk.density*(Percent.C/100))

soil.c$int.volume = ((Int.B/100)-(Int.A/100))*10000

# Calculate the depth of the bottom-most interval and apply the fifth c.density value to it

for(k in 1:nrow(soil.c))
    if(Interval[k] == 4)
    {soil.c$int.volume[k] = ((0.5)*10000)} else
      if(Interval[k] == 5) 
      {soil.c$int.volume[k] = ((Avg.depth[k]/100) - (100/100))*10000} else
      {soil.c$int.volume[k] = ((Int.B[k]/100)-(Int.A[k]/100))*10000}

soil.c$C.per.ha = soil.c$int.volume*soil.c$C.density
head(soil.c)

subplot.soil.c = aggregate(C.per.ha ~ Subplot + Plot, data = soil.c, FUN="sum")
plot.soil.c = aggregate(C.per.ha ~ Plot, data=subplot.soil.c, FUN="mean")

head(plot.soil.c)

## Source the SummarySE function (borrowed from: 
##    http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_%28ggplot2%29/)

source("~/Dropbox/mangrove-work/model-files/Function-SummarySE.R")

plot.soil.summary = summarySE(subplot.soil.c, measurevar="C.per.ha", groupvars=c("Plot"))

ggplot(plot.soil.summary, aes(x=Plot, y=C.per.ha)) + 
  geom_errorbar(aes(ymin=C.per.ha-se, ymax=C.per.ha+se), width=.1) + geom_point() + 
  ggtitle("Soil Carbon Plot Means +/- 1 St. Error")

## Repeat summary stats for all plots combined, i.e. no grouping variable (to cut down on error)

site.soil.summary = summarySE(subplot.soil.c, measurevar="C.per.ha")

#------------------------------------------------------------------------------
#Add other relevant parameters to characterize 
#the plot (i.e. mean DBH and Tons ag.biomass)

Avg.DBH =round(aggregate(DBH.cm ~ Plot, data=trees, FUN="mean"),1)
Mg.ag.biomass = round(ag.biomass.per.ha/1000,1)
Mg.bg.biomass = round(bg.biomass.per.ha/1000,1)
Mg.biomass = round(biomass.per.ha/1000,1)

summary = cbind.data.frame(site, imp.values$Plot, imp.values$Species, imp.values$Rel.Importance, 
                Avg.DBH$DBH.cm, BAs.per.ha, Mg.ag.biomass, Mg.bg.biomass, Mg.biomass, cwd.mass$Total,
                plot.soil.c$C.per.ha)

colnames(summary) = c("Site", "Plot","Dom.Sps","Imp. Value (%)","Avg.DBH (cm)", "BA (m^2/ha)",
                      "AGB (Mg/ha)", "BGB (Mg/ha)", "Total Biomass (Mg/ha)", 
                      "CWD (Mg/ha)", "SOC (Mg/ha)")

summary = as.data.frame(summary)

print(summary)

if(exists("total.summary") == F) {total.summary = summary} else
  {total.summary = rbind(total.summary,summary)}

#Summarize the data by dominant genera

sps.summary <- aggregate(. ~ Dom.Sps, data=summary, FUN="mean")

#-----------------------------------------------------------------------------------------------

