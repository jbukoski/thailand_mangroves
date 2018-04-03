## Generate a look up table for the allometric equations

library(tidyverse)

# Initialize the table

allom_lookup <- tibble(genus = character(), 
                       species = character(), 
                       sps_code = character(),
                       density = double(),
                       ag_form = list(), 
                       bg_form = list(), 
                       ag_ref = character(),
                       bg_ref = character()
                       )

# Load the allometric equations

source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/allometry_equations.R")


# Build general Komiyama et al 2005 functions

allom_lookup <- allom_lookup %>%
  add_row(genus = "aegiceras", species = "corniculatum", sps_code = "aeco", density = NA,
          ag_form = c(NA), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "ren2010", bg_ref = "komiyama2005") %>%
  add_row(genus = "avicennia", species = "alba", sps_code = "aval", density = 0.587, 
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "avicennia", species = "marina", sps_code = "avma", density = 0.65,
          ag_form = c(avma_ag_comley2005), 
          bg_form = c(avma_bg_comley2005), 
          ag_ref = "comley2005", bg_ref = "comley2005") %>%
  add_row(genus = "avicennia", species = "officinalis", sps_code = "avof", density = 0.605,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "cylindrica", sps_code = "brcy", density = 0.720,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "gymnorrhiza", sps_code = "brgy", density = 0.710,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "parviflora", sps_code = "brpa", density = 0.760,
          ag_form = c(function(dbh, dens){10^(-0.7045+2.5336*(log10(dbh)))}), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = NA, bg_ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "sexangula", sps_code = "brse", density = 0.740,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "excoecaria", species = "agallocha", sps_code = "exag", density = 0.416,
          ag_form = c(function(dbh, dens){exp(1)^(1.0996*(log(dbh^2))-0.8572)}), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "hossain2015", bg_ref = "komiyama2005") %>%
  add_row(genus = "kandelia", species = "obovata", sps_code = "kaob", density = 0.525,
          ag_form = c(function(dbh, dens){0.251*0.525*(dbh^2.46)}), 
          bg_form = c(general_bg_komiyama2005),
          ag_ref = "khan????", bg_ref = "komiyama2005") %>%
  add_row(genus = "lumnitzera", species = "racemosa", sps_code = "lura", density = 0.710,
          ag_form = c(function(dbh, dens){1.184+(dbh)^2.384}), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "kangkuso2015", bg_ref = "komiyama2005") %>%
  add_row(genus = "rhizophora", species = "apiculata", sps_code = "rhap", density = 0.850,
          ag_form = c(rhap_ag_ong2004), 
          bg_form = c(rhap_bg_ong2004), 
          ag_ref = "ong2004", bg_ref = "ong2004") %>%
  add_row(genus = "rhizophora", species = "mucronata", sps_code = "rhmu", density = 0.821,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "rhizophora", species = "stylosa", sps_code = "rhst", density = 0.840,
          ag_form = c(rhst_ag_comley2005), 
          bg_form = c(rhst_bg_comley2005), 
          ag_ref = "comley2005", bg_ref = "comley2005") %>%
  add_row(genus = "sonneratia", species = "caseolaris", sps_code = "soca", density = 0.389,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "sonneratia", species = "ovata", sps_code = "soov", density = 0.850,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "scyphiphora", species = "hydrophyllacea", sps_code = "schy", density = 0.90,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "xylocarpus", species = "granatum", sps_code = "xygr", density = 0.567,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005") %>%
  add_row(genus = "xylocarpus", species = "moluccensis", sps_code = "xymo", density = 0.611,
          ag_form = c(general_ag_komiyama2005), 
          bg_form = c(general_bg_komiyama2005), 
          ag_ref = "komiyama2005", bg_ref = "komiyama2005")


# ## Generate test data to examine the lookup table
# 
# test_data <- tibble(species = c("aval", "brcy", "brgy", "avma", "aval"),
#                diameter = c(5.5, 1.3, 19.1, 8.6, 7.7))
# 
# 
# ## Use tidyverse functions to map allometry equations to tree measurements
# 
# biom <- test_data %>%
#   left_join(allom_lookup, by=c("species" = "sps_code")) %>%
#   mutate(params = map2(diameter, density, list)) %>%
#   mutate(agb = invoke_map_dbl(ag_form, params))


# #manually enter species wood densities here, **none available for KO at the moment
# sps = rbind("AA", "AM", "AO", "BC", "BG", "BP", "BS", "EA", "KO", "LR", "RA", "RM", "RS", "SH", "SC", "SO", "XG", "XM")
# densities = c((0.56+0.67+0.53)/3, 0.650, (0.59+0.62)/2, 0.720, (0.66+0.76)/2, (0.74+0.78)/2, 0.740, 
#               (0.39+0.48+0.379)/3, (0.512+0.460+0.557+0.57)/4, 0.710, 0.85, (0.74+0.82+0.904)/3, 
#               0.84,0.9,(0.387+0.390)/2, 0.850, (0.557+0.525+0.62)/3, 0.611)
# sps.densities= as.data.frame(cbind(sps,densities))
# colnames(sps.densities)=c("Species", "Density")
# sps.densities$Density = as.numeric(as.character(sps.densities$Density))

