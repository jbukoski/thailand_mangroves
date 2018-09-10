# helper functions

## A helper function to identify genus and species from the species code

id_taxon <- function(df, species_code) {
  
  df %>%
    mutate(genus = ifelse(species_code %in% c("AA", "AM", "AO"), "avicennia",
                          ifelse(species_code %in% c("BC", "BG", "BP", "BS"), "bruguiera", 
                          ifelse(species_code %in% c("EA"), "excoecaria",
                          ifelse(species_code %in% c("LR"), "lumnitzera",
                          ifelse(species_code %in% c("RA", "RM"), "rhizophora",
                          ifelse(species_code %in% c("SH"), "scyphiphora",
                          ifelse(species_code %in% c("SO"), "sonneratia",
                          ifelse(species_code %in% c("XG", "XM"), "xylocarpus", "unknown"))))))))) %>%
    mutate(species = ifelse(species_code == "AA", "alba",
                            ifelse(species_code == "AM", "marina",
                            ifelse(species_code == "AO", "officinalis",
                            ifelse(species_code == "BC", "cylindrica",
                            ifelse(species_code == "BG", "gymnorrhiza",
                            ifelse(species_code == "BP", "parviflora",
                            ifelse(species_code == "BS", "sexangula",
                            ifelse(species_code == "EA", "agallocha",
                            ifelse(species_code == "LR", "racemosa",
                            ifelse(species_code == "RA", "apiculata",
                            ifelse(species_code == "RM", "mucronata",
                            ifelse(species_code == "SH", "hydrophyllacea",
                            ifelse(species_code == "SO", "ovata",
                            ifelse(species_code == "XG", "granatum",
                            ifelse(species_code == "XM", "moluccensis", "unknown"))))))))))))))))

}

# A helper function to remove speckling

rm_speckling <- function(raster, speckle_size, dir) {
  rclump <- clump(raster, directions = dir)
  clumpFreq <- as.data.frame(freq(rclump))
  drop_idx <- clumpFreq$value[which(clumpFreq$count < speckle_size)]
  
  rclumpSieve <- rclump
  rclumpSieve[rclump %in% drop_idx] <- NA
  
  mat <- matrix(c(NA, NA, 0, -Inf, Inf, 1), ncol = 3, byrow = TRUE)
  rast <- reclassify(rclumpSieve, mat)
  
  return(rast)
}

# A function for a summary that provides standard deviation, standard error, and 95% confidence interval

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

