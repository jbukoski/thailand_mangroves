
## Generate a look up table for the allometric equations

# Initialize the table

allom_lookup <- tibble(genus = character(), 
                       species = character(), 
                       sps_code = character(),
                       density = double(),
                       ag_form = list(), 
                       bg_form = list(), 
                       ref = character()
                       )

# Build general Komiyama et al 2005 functions

komi2005ag <- function(dbh, dens){
  0.251*dens*(dbh^2.46)
}

komi2005bg <- function(dbh, dens){
  0.199*(dens^0.899)*(dbh^2.22)
}


allom_lookup <- allom_lookup %>%
  add_row(genus = "aegiceras", species = "corniculatum", sps_code = "aeco", density = NA,
          ag_form = c(NA), bg_form = c(NA), ref = "ren2010") %>%
  add_row(genus = "avicennia", species = "alba", sps_code = "aval", density = 0.506, 
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "avicennia", species = "marina", sps_code = "avma", density = NA,
          ag_form = c(function(dbh, dens){0.1848*(dbh^2.3524)}), bg_form = NA, ref = "clough1997") %>%
  add_row(genus = "avicennia", species = "officinalis", sps_code = "avof", density = 0.670,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "cylindrica", sps_code = "brcy", density = 0.749,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "gymnorrhiza", sps_code = "brgy", density = 0.699,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "bruguiera", species = "parviflora", sps_code = "brpa", density = NA,
          ag_form = c(function(dbh, dens){10^(-0.7045+2.5336*(log10(DBH.cm[k])))}), bg_form = c(NA), ref = NA) %>%
  add_row(genus = "bruguiera", species = "sexangula", sps_code = "brse", density = 0.808,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "excoecaria", species = "agallocha", sps_code = "exag", density = NA,
          ag_form = c(function(dbh, dens){exp(1)^(1.0996*(log(dbh^2))-0.8572)}), bg_form = NA, ref = "hossain2015") %>%
  add_row(genus = "kandelia", species = "obovata", sps_code = "kaob", density = NA,
          ag_form = c(function(dbh, dens){0.251*0.525*(dbh^2.46)}), bg_form = NA, ref = "khan????") %>%
  add_row(genus = "lumnitzera", species = "racemosa", sps_code = "lura", density = NA,
          ag_form = c(function(dbh, dens){1.184+(DBH.cm[k])^2.384}), bg_form = NA, ref = "kangkuso2015") %>%
  add_row(genus = "rhizophora", species = "apiculata", sps_code = "rhap", density = NA,
          ag_form = c(function(dbh, dens){exp(1)^((2.318*log(dbh)) - 1.671)}), bg_form = NA, ref = "ong2004") %>%
  add_row(genus = "rhizophora", species = "mucronata", sps_code = "rhmu", density = 0.701,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "rhizophora", species = "stylosa", sps_code = "rhst", density = NA,
          ag_form = c(function(dbh, dens){10^(-0.6564+2.4292*(log10(dbh)))}), bg_form = NA, ref = "clough1997") %>%
  add_row(genus = "sonneratia", species = "caseolaris", sps_code = "soca", density = 0.340,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "scyphiphora", species = "hydrophyllacea", sps_code = "schy", density = 0.90,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "xylocarpus", species = "granatum", sps_code = "xygr", density = 0.528,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005") %>%
  add_row(genus = "xylocarpus", species = "moluccensis", sps_code = "xymo", density = 0.531,
          ag_form = c(komi2005ag), bg_form = c(komi2005bg), ref = "komiyama2005")


## Generate test data to examine the lookup table

test_data <- tibble(species = c("aval", "brcy", "brgy", "avma", "aval"),
               diameter = c(5.5, 1.3, 19.1, 8.6, 7.7))


## Use tidyverse functions to map allometry equations to tree measurements

biom <- data %>%
  left_join(select(allom_lookup, sps_code, density, ag_form, bg_form), 
            by=c("species" = "sps_code")) %>%
  mutate(params = map2(diameter, density, list)) %>%
  mutate(agb = invoke_map_dbl(ag_form, params))
