# helper functions

## A helper function to identify genus and species from the species code

id_taxon <- function(df, species_code) {
  
  df %>%
    mutate(genus = ifelse(species_code %in% c("AA", "AM", "AO"), "avicennia",
                          ifelse(species_code %in% c("BC", "BG", "BP"), "bruguiera", 
                                 ifelse(species_code %in% c("EA"), "excoecaria",
                                        ifelse(species_code %in% c("RA", "RM"), "rhizophora",
                                               ifelse(species_code %in% c("SH"), "scyphiphora", 
                                                      ifelse(species_code %in% c("XG", "XM"), "xylocarpus", "unknown"))))))) %>%
    mutate(species = ifelse(species_code == "AA", "alba",
                            ifelse(species_code == "AM", "marina",
                                   ifelse(species_code == "AO", "officinalis",
                                          ifelse(species_code == "BC", "cylindrica",
                                                 ifelse(species_code == "BG", "gymnorrhiza",
                                                        ifelse(species_code == "BP", "parviflora",
                                                               ifelse(species_code == "EA", "agallocha",
                                                                      ifelse(species_code == "RA", "apiculata",
                                                                             ifelse(species_code == "RM", "mucronata",
                                                                                    ifelse(species_code == "SH", "hydrophyllacea",
                                                                                           ifelse(species_code == "XG", "granatum",
                                                                                                  ifelse(species_code == "XM", "moluccensis", "unknown")))))))))))))

}


test_df <- tibble(sps = c("AA", "AM", "AA", "AA", "RM", "AO"))

id_taxon(test_df, test_df$sps)

