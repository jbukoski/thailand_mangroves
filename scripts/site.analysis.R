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

#------------------------------------------------------------------------------
# Clean up the dataset.

trees <- all.trees

source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/helper_funcs.R")

colnames(trees) <- tolower(colnames(trees))

trees <- trees %>% 
  id_taxon(trees$species) %>%
  mutate(sps_code = paste0(substr(genus, 1, 2), substr(species, 1, 2)))


#------------------------------------------------------------------------------
#Calculate above-ground biomass using species-specific allometric equations; where species-specific
#equations are not available, used Komiyama et al 2005 general equation with species specific wood
#densities

source("/home/jbukoski/manuscripts/thailand_stocks/thailand_mangroves/scripts/allometry.R")

test_trees <- trees %>%
  left_join(allom_lookup, by = c("sps_code" = "sps_code")) %>%
  mutate(params = map2(dbh.cm, density, list)) %>%
  mutate(agb = invoke_map_dbl(ag_form, params))

#-------
#old code

trees <- all.trees;
  
  ag.biomass = rep(0,nrow(trees))
  leaf.biomass = rep(0,nrow(trees))
  prop.biomass = rep(0,nrow(trees))
  trees = trees[,1:13]
  trees = cbind(trees,ag.biomass,leaf.biomass,prop.biomass)
  head(trees)
  
  attach(trees)
  
  for(k in 1:nrow(trees)) #have to adjust this to remove rows with NA (excel read in extra lines)
    if(Species[k] == "AA") #Avicennia alba, general equation w/ AA density
    {trees$ag.biomass[k] = 0.251*0.506*((DBH.cm[k])^2.46)} else
      if(Species[k] == "AC") #Aegiceras corniculatum, Tam et al "total aerial" equation
      {trees$ag.biomass[k] = 0.280*((DBH.cm[k]^2*Height.m[k])^0.693)} else # Equation from Ren et al. 2010
        if(Species[k] == "AM") #Avicennia marina, 2 equations!
        {trees$ag.biomass[k] = 0.1848*((DBH.cm[k])^2.3524)} else
        # {trees$ag.biomass[k] = 10^(-0.7506+2.2990*(log(DBH.cm[k], base=10)))}  #Clough et al. 1997 (for VN sites)
          if(Species[k] == "AO") #Avicennia officinalis, general equation
          {trees$ag.biomass[k] = 0.251*0.670*((DBH.cm[k])^2.46)} else
            if(Species[k] == "BC") #Bruguiera cylindrica, general equation w/ BC density
            {trees$ag.biomass[k] = 0.251*0.749*((DBH.cm[k])^2.46)} else
              if(Species[k] == "BG") #Bruguiera gymnorrhiza, gen. equation w/ BG density
              {trees$ag.biomass[k] = 0.251*0.699*((DBH.cm[k])^2.46)} else
                if(Species[k] == "BP") #Bruguiera parviflora, where is this from?
                {trees$ag.biomass[k] = 10^(-0.7045+2.5336*(log10(DBH.cm[k])))} else
                  if(Species[k] == "BS") #Bruguiera sexangula, gen. equation w/ BG density
                  {trees$ag.biomass[k] = 0.251*0.808*((DBH.cm[k])^2.46)} else
                    if(Species[k] == "EA") #Excoecaria agallocha; Note: eq. from Hossain et al, 2015
                    {trees$ag.biomass[k] = exp(1)^(1.0996*(log((DBH.cm[k])^2))-0.8572)} else
                      if(Species[k] == "KO") #Kandelia obovata, Khan equation using height and diameter
                      {trees$ag.biomass[k] = 0.251*0.525*((DBH.cm[k])^2.46)} else
                        #{trees$ag.biomass[k] = 3.203*(10^-2)*((DBH.cm[k]^2*Height.m[k])^1.058)} else
                        if(Species[k] == "LR") #Lumnitzera racemosa, general equation using diam
                        # {trees$ag.biomass[k] = 1.788+(2.529*(log(DBH.cm[k], base=10)))} else
                        {trees$ag.biomass[k] = 1.184+(DBH.cm[k])^2.384} else  #Kangkuso et al. 2015
                          if(Species[k] == "RA") #Rhizophora apiculata (Ong 2004, total AGB)
                          {trees$ag.biomass[k] = exp(1)^((2.318*log(DBH.cm[k])) - 1.671)} else
                          # trees$leaf.biomass[k] = 10^(-1.8571+(2.1072*(log(DBH.cm[k],base=10))))} else
                            if(Species[k] == "RM") #Rhizophora mucronata, gen equation w/ RM density
                            {trees$ag.biomass[k] = 0.251*0.701*((DBH.cm[k])^2.46)} else
                              if(Species[k] == "RS") #Rhizophora stylosa, Clough et al 1997; Australia
                              {trees$ag.biomass[k] = 10^(-0.6564+2.4292*(log(DBH.cm[k], base=10)))} else
                                if(Species[k] == "SC") #Sonneratia caseolaris, Komiyama et al 2005 generic
                                {trees$ag.biomass[k] = 0.251*0.340*((DBH.cm[k])^2.46)} else
                                  if(Species[k] == "SH") #Scyphiphora hydrophyllacea; general equation
                                  {trees$ag.biomass[k] = 0.251*0.9*((DBH.cm[k])^2.46)} else
                                    if(Species[k] == "SO") #Shorea obtusa, Komiyama et al 2005 general
                                    {trees$ag.biomass[k] = 0.251*0.850*((DBH.cm[k])^2.46)} else
                                      if(Species[k] == "XG") #Xylocarpus granatum, gen equation w/ XG density
                                      {trees$ag.biomass[k] = 0.251*0.528*((DBH.cm[k])^2.46)} else 
                                        if(Species[k] == "XM") #Xylocarpus moluccensis, gen eq. w/ XM density
                                        {trees$ag.biomass[k] = 0.251*0.531*((DBH.cm[k])^2.46)}
  
  for(k in 1:nrow(trees))
    if(Species[k] == "RA" & DBH.cm[k] <= 5)
    {trees$prop.biomass[k] = 0.101*trees$ag.biomass[k]} else
      if(Species[k] == "RA" & DBH.cm[k] > 5 & DBH.cm[k] <= 10)
      {trees$prop.biomass[k] = 0.204*trees$ag.biomass[k]} else
        if(Species[k] == "RA" & DBH.cm[k] > 10 & DBH.cm[k] <= 15)
        {trees$prop.biomass[k] = 0.356*trees$ag.biomass[k]} else
          if(Species[k] == "RA" & DBH.cm[k] > 15 & DBH.cm[k] <= 20)
          {trees$prop.biomass[k] = 0.273*trees$ag.biomass[k]} else
            if(Species[k] == "RA" & DBH.cm[k] > 20)
            {trees$prop.biomass[k] = 0.210*trees$ag.biomass[k]}
  
  for(k in 1:nrow(trees))
    if(Species[k] == "RA")
    {trees$ag.biomass[k] = trees$ag.biomass[k] + trees$leaf.biomass[k] + trees$prop.biomass[k]}
              
  rm(ag.biomass, leaf.biomass, prop.biomass)  #Clean things up, rm empty vectors
  head(trees)  #Check to see if ag.biomass values are calculated
  
  #------------------------------------------------------------------------------------
  # Adjust ag.biomass variable based on Status variable
  
  trees$top.diam = rep(0, nrow(trees))
  trees$stump.vol = rep(0, nrow(trees))
  
  for(k in 1:nrow(trees))
    if(is.na(Status[k]) == "TRUE")
    {trees$ag.biomass.new[k] = trees$ag.biomass[k]} else
      if(Status[k] == "1")
      {trees$ag.biomass.new[k] = 0.95*trees$ag.biomass[k]} else #Reduce biomass by 5% for loss of leaves
        if(Status[k] == "2")
        {trees$ag.biomass.new[k] = 0.8*trees$ag.biomass[k]} else  #Reduce biomass by 20% for loss of branches
          if(Status[k] == "3" & is.na(Base.cm[k])=="FALSE" & Height.m[k] >= 1.37)
            {trees$top.diam[k] = Base.cm[k]-((100*Height.m[k])*((Base.cm[k]-DBH.cm[k])/137))
            trees$stump.vol[k] = ((pi*100*Height.m[k])/12)*
              (Base.cm[k]^2+top.diam[k]^2+(Base.cm[k]*top.diam[k]))
            trees$ag.biomass.new[k] = (0.69*trees$stump.vol[k])/1000} else
            if(Status[k] == "3" & is.na(Base.cm[k])=="FALSE" & Height.m[k] < 1.37)
              {trees$top.diam[k] = Base.cm[k]-((100*Height.m[k])*((Base.cm[k]-DBH.cm[k])/(100*Height.m[k])))
              trees$stump.vol[k] = ((pi*100*Height.m[k])/12)*
                  (Base.cm[k]^2+top.diam[k]^2+(Base.cm[k]*top.diam[k]))
              trees$ag.biomass.new[k] = (0.69*trees$stump.vol[k])/1000} #Using wood density of large CWD from K&D 2012
  
  trees$ag.biomass = trees$ag.biomass.new
  trees = trees[,1:14]
  head(trees)
  
  #-------------------------------------------------------------------------------
  # Calculate belowground biomass using Komiyama's general equation 
  # and species specific wood densities; *NOTE: C content of roots lower than AG parts, should
  # multiply by 39% to convert biomass to C content; densities taken from the Global Wood Density 
  # Database; multiple values are averaged from SE Asia; KO not available - averaged across KC
  
  #manually enter species wood densities here, **none available for KO at the moment
  sps = rbind("AA", "AM", "AO", "BC", "BG", "BP", "BS", "EA", "KO", "LR", "RA", "RM", "RS", "SH", "SC", "SO", "XG", "XM")
  densities = c((0.56+0.67+0.53)/3, 0.650, (0.59+0.62)/2, 0.720, (0.66+0.76)/2, (0.74+0.78)/2, 0.740, 
                (0.39+0.48+0.379)/3, (0.512+0.460+0.557+0.57)/4, 0.710, 0.85, (0.74+0.82+0.904)/3, 
                0.84,0.9,(0.387+0.390)/2, 0.850, (0.557+0.525+0.62)/3, 0.611)
  sps.densities= as.data.frame(cbind(sps,densities))
  colnames(sps.densities)=c("Species", "Density")
  sps.densities$Density = as.numeric(as.character(sps.densities$Density))
  
  #-------------------------------------------------------------------------------
  #Create a vector to be filled with the belowground biomass data; uses Komiyama et al general equation 
  #with species densities
  
  bg.biomass <- rep(0, nrow(trees))   
  trees <- cbind(trees, bg.biomass)
  
  #Adjust levels of trees data frame to match levels of species densities; avoids error below
  trees$Species = factor(trees$Species, 
                         levels=levels(sps.densities$Species))
  
  trees <- trees %>%
    left_join(sps.densities, by = "Species") %>%
    mutate(bg.biomass = ifelse(Species == "KO", 0.745*(ag.biomass^0.810), 
                              ifelse(Species == "RA", 0.00698*(DBH.cm^2.61), 
                                    ifelse(Species == "RS", 10^(-0.583*(log10(DBH.cm)^1.86)), 0.199*(Density^0.899)*(DBH.cm^2.22))))) %>%
    mutate(ag.biomass = round(ag.biomass, 2), bg.biomass = round(bg.biomass, 2))
  
  rm(bg.biomass)
  
  #-------------------------------------------------------------------------------
  #Create a vector for total biomass; sum of ag.biomass and bg.biomass
  
  trees <- trees %>%
    mutate(biomass = ag.biomass + bg.biomass)
  
  head(trees)
  
  #-------------------------------------------------------------------------------
  #Calculate the biomass of saplings
  saps = saplings
  
  sap.AGB = rep(0,nrow(saps))
  sap.leaf = rep(0,nrow(saps))
  sap.prop = rep(0,nrow(saps))
  saps = saps[,1:10]
  saps = cbind.data.frame(saps,sap.AGB,sap.leaf,sap.prop)
  head(saps)
  
  attach(saps)
  
  for(k in 1:nrow(saps)) #have to adjust this to remove rows with NA (excel read in extra lines)
    if(Species[k] == "AA") #Avicennia alba, general equation w/ AA density
    {saps$sap.AGB[k] = 0.251*0.506*((DBH.cm[k])^2.46)} else
      if(Species[k] == "AC") #Aegiceras corniculatum, Tam et al "total aerial" equation
      {saps$sap.AGB[k] = 0.280*((DBH.cm[k]^2*Height.m[k])^0.693)} else # Equation from Ren et al. 2010
        if(Species[k] == "AM") #Avicennia marina, 2 equations!
        {saps$sap.AGB[k] = 0.1848*((DBH.cm[k])^2.3524)} else
          # {saps$sap.AGB[k] = 10^(-0.7506+2.2990*(log(DBH.cm[k], base=10)))}  #Clough et al. 1997 (for VN sites)
          if(Species[k] == "AO") #Avicennia officinalis, general equation
          {saps$sap.AGB[k] = 0.251*0.670*((DBH.cm[k])^2.46)} else
            if(Species[k] == "BC") #Bruguiera cylindrica, general equation w/ BC density
            {saps$sap.AGB[k] = 0.251*0.749*((DBH.cm[k])^2.46)} else
              if(Species[k] == "BG") #Bruguiera gymnorrhiza, gen. equation w/ BG density
              {saps$sap.AGB[k] = 0.251*0.699*((DBH.cm[k])^2.46)} else
                if(Species[k] == "BP") #Bruguiera parviflora, where is this from?
                {saps$sap.AGB[k] = 10^(-0.7045+2.5336*(log(DBH.cm[k], base=10)))} else
                  if(Species[k] == "BS") #Bruguiera sexangula, gen. equation w/ BG density
                  {saps$sap.AGB[k] = 0.251*0.808*((DBH.cm[k])^2.46)} else
                    if(Species[k] == "EA") #Excoecaria agallocha; Note: eq. from Hossain et al, 2015
                    {saps$sap.AGB[k] = exp(1)^(1.0996*(log((DBH.cm[k])^2))-0.8572)} else
                      if(Species[k] == "KO") #Kandelia obovata, general equation using height and diameter
                      {saps$sap.AGB[k] = 0.251*0.525*((DBH.cm[k])^2.46)} else
                      #{saps$sap.AGB[k] = 3.203*(10^-2)*((DBH.cm[k]^2*Height.m[k])^1.058)} else
                        if(Species[k] == "LR") #Lumnitzera racemosa, general equation using height and diameter (in )
                        # {saps$sap.AGB[k] = 1.788+(2.529*(log(DBH.cm[k], base=10)))} else
                        {saps$sap.AGB[k] = 1.184+(DBH.cm[k])^2.384} else  #Kangkuso et al. 2015
                          if(Species[k] == "RA") #Rhizophora apiculata (Ong 2004, total AGB)
                          {saps$sap.AGB[k] = exp(1)^((2.318*log(DBH.cm[k])) - 1.671)} else
                            # saps$sap.leaf[k] = 10^(-1.8571+(2.1072*(log(DBH.cm[k],base=10))))} else
                            if(Species[k] == "RM") #Rhizophora mucronata, gen equation w/ RM density
                            {saps$sap.AGB[k] = 0.251*0.701*((DBH.cm[k])^2.46)} else
                              if(Species[k] == "RS") #Rhizophora stylosa, Clough et al 1997; Australia
                              {saps$sap.AGB[k] = 10^(-0.6564+2.4292*(log(DBH.cm[k], base=10)))} else
                                if(Species[k] == "SC") #Sonneratia caseolaris, Komiyama et al 2005 generic
                                {saps$sap.AGB[k] = 0.251*0.340*((DBH.cm[k])^2.46)} else
                                  if(Species[k] == "SH") #Scyphiphora hydrophyllacea; general equation
                                  {trees$ag.biomass[k] = 0.251*0.9*((DBH.cm[k])^2.46)} else
                                    if(Species[k] == "SO") #Shorea obtusa, Komiyama et al 2005 general
                                    {saps$sap.AGB[k] = 0.251*0.850*((DBH.cm[k])^2.46)} else
                                      if(Species[k] == "XG") #Xylocarpus granatum, gen equation w/ XG density
                                      {saps$sap.AGB[k] = 0.251*0.528*((DBH.cm[k])^2.46)} else 
                                        if(Species[k] == "XM") #Xylocarpus moluccensis, gen eq. w/ XM density
                                        {saps$sap.AGB[k] = 0.251*0.531*((DBH.cm[k])^2.46)}
  
  for(k in 1:nrow(saps))
    if(Species[k] == "RA" & DBH.cm[k] <= 5)
    {saps$sap.prop[k] = 0.101*saps$sap.AGB[k]} else
      if(Species[k] == "RA" & DBH.cm[k] > 5 & DBH.cm[k] <= 10)
      {saps$sap.prop[k] = 0.204*saps$sap.AGB[k]} else
        if(Species[k] == "RA" & DBH.cm[k] > 10 & DBH.cm[k] <= 15)
        {saps$sap.prop[k] = 0.356*saps$sap.AGB[k]} else
          if(Species[k] == "RA" & DBH.cm[k] > 15 & DBH.cm[k] <= 20)
          {saps$sap.prop[k] = 0.273*saps$sap.AGB[k]} else
            if(Species[k] == "RA" & DBH.cm[k] > 20)
            {saps$sap.prop[k] = 0.210*saps$sap.AGB[k]}
  
  for(k in 1:nrow(saps))
    if(Species[k] == "RA")
    {saps$sap.AGB[k] = saps$sap.AGB[k] + saps$sap.leaf[k] + saps$sap.prop[k]}
  
  rm(sap.AGB,sap.leaf,sap.prop)  #Clean things up, rm empty vector
  head(saps)  #Check to see if sap.AGB values are calculated
  
  #------------------------------------------------------------------------------------
  # Adjust sap.AGB variable based on Status variable
  
  saps$top.diam=rep(0,nrow(saps))
  saps$stump.vol=rep(0,nrow(saps))
  
  attach(saps)
  
  for(k in 1:nrow(saps))
    if(is.na(Status[k]) == "TRUE")
    {saps$sap.AGB.new[k] = saps$sap.AGB[k]} else
      if(Status[k] == "1")
      {saps$sap.AGB.new[k] = 0.95*saps$sap.AGB[k]} else #Reduce biomass by 5% for loss of leaves
        if(Status[k] == "2")
        {saps$sap.AGB.new[k] = 0.8*saps$sap.AGB[k]} else  #Reduce biomass by 20% for loss of branches
          if(Status[k] == "3" & is.na(Base.cm[k])=="FALSE" & Height.m[k] >= 1.37)
          {saps$top.diam[k] = Base.cm[k]-((100*Height.m[k])*((Base.cm[k]-DBH.cm[k])/137))
          saps$stump.vol[k] = ((pi*100*Height.m[k])/12)*
            (Base.cm[k]^2+top.diam[k]^2+(Base.cm[k]*top.diam[k]))
          saps$sap.AGB.new[k] = (0.69*saps$stump.vol[k])/1000} else
            if(Status[k] == "3" & is.na(Base.cm[k])=="FALSE" & Height.m[k] < 1.37)
            {saps$top.diam[k] = Base.cm[k]-((100*Height.m[k])*((Base.cm[k]-DBH.cm[k])/(100*Height.m[k])))
            saps$stump.vol[k] = ((pi*100*Height.m[k])/12)*
              (Base.cm[k]^2+top.diam[k]^2+(Base.cm[k]*top.diam[k]))
            saps$sap.AGB.new[k] = (0.69*saps$stump.vol[k])/1000} #Using wood density of large CWD from K&D 2012
  
  saps$sap.AGB = saps$sap.AGB.new
  saps = saps[,1:14]
  head(saps)
  
  #-------------------------------------------------------------------------------
  #Calculate sapling belowground biomass; Create a vector to be filled with the 
  #sapling belowground biomass data; uses Komiyama et al general equation with 
  #species densities
  
  saps.BGB = rep(0,nrow(saps))   
  saps = cbind(saps,saps.BGB)
  
  #Adjust levels of saps data frame to match levels of species densities; avoids error below
  saps$Species = factor(saps$Species, levels=levels(sps.densities$Species))
  
  for(k in 1:nrow(saps))
    for(j in 1:nrow(sps.densities))
        if(saps$Species[k] == "RA") # From Ong, 2004
        {saps$saps.BGB[k] = 0.00698*(DBH.cm[k]^2.61)} else
          if(saps$Species[k] == "RS")
          {saps$saps.BGB[k] = 10^(-0.583*(log(DBH.cm[k],base=10)^1.86))} else
            if(saps$Species[k] == sps.densities$Species[j] & saps$Species[k] != c("RA","RS"))
            {saps$saps.BGB[k] = 0.199*(sps.densities$Density[j]^0.899)*
              ((saps$DBH.cm)[k]^2.22)}
  
  rm(saps.BGB)
  
  attach(saps)
  saps$sap.AGB = round(sap.AGB,2)
  saps$saps.BGB = round(saps.BGB,2)
  
  #-------------------------------------------------------------------------------
  #Create a vector for total biomass; sum of ag.biomass and bg.biomass
  
  for(k in 1:nrow(saps))
  {saps$biomass[k] = round(sap.AGB[k] + saps.BGB[k], 2)}
  
  head(saps)
  
  #-------------------------------------------------------------------------------
  #Calculate biomass for the coarse-woody debris pool based on default mean diameters 
  #and densities given in the Kauffman and Donato protocols
  
  #define the mean specific gravities (g/cm^3) of the wood classes; taken from K&D, 2012
  cwd.gravs = as.data.frame(cbind(0.48, 0.64, 0.71, 0.69))
  cwd.diams = as.data.frame(cbind(0.43, 1.47, 4.52, NA))
  cwd.params = rbind(cwd.gravs, cwd.diams)
  
  colnames(cwd.params)=c("Fine", "Small", "Medium", "Large")
  rownames(cwd.params)=c("Density", "Mean.Diams")
  
  #Calculating the mass and volume per plot for each of the size classes, beginning with fine CWD
  
  fine.cts = aggregate(Fine ~ Plot, data=cwd, sum)
  fine.volume = rep(0,nrow(fine.cts)); fine.mass = rep(0, nrow(fine.cts))
  fine.cts=cbind(fine.cts,fine.volume,fine.mass)
  
  for(k in 1:nrow(fine.cts))
    {fine.cts$fine.volume[k]=pi^2*((fine.cts$Fine[k]*(cwd.params$Fine[2]^2))/(8*2*20))
      fine.cts$fine.mass[k]=fine.cts$fine.volume[k]*cwd.params$Fine[1]}
  
  #Repeating here for small CWD
  
  small.cts = aggregate(Small ~ Plot, data=cwd, sum)
  small.volume = rep(0,nrow(small.cts)); small.mass = rep(0, nrow(small.cts))
  small.cts=cbind(small.cts,small.volume,small.mass)
  
  for(k in 1:nrow(small.cts))
    {small.cts$small.volume[k]=pi^2*((small.cts$Small[k]*(cwd.params$Small[2]^2))/(8*3*20))
      small.cts$small.mass[k]=small.cts$small.volume[k]*cwd.params$Small[1]}
  
  #Repeating here for medium CWD
  
  medium.cts = aggregate(Medium ~ Plot, data=cwd, sum)
  medium.volume = rep(0,nrow(medium.cts)); medium.mass = rep(0, nrow(medium.cts))
  medium.cts=as.data.frame(cbind(medium.cts,medium.volume,medium.mass))
  
  for(k in 1:nrow(medium.cts))
  {medium.cts$medium.volume[k]=pi^2*((medium.cts$Medium[k]*(cwd.params$Medium[2]^2))/(8*5*20))
  medium.cts$medium.mass[k]=medium.cts$medium.volume[k]*cwd.params$Medium[1]}
  
  #Repeating here for large, rotten CWD, though slightly different formula
  
  large.rot.volume = rep(0,nrow(cwd))
  
  for(k in 1:nrow(cwd))
  {large.rot.volume[k]=pi^2*(sum(cwd$Large.rotten[k]^2)/(8*12*20))}
  
  large = cbind(cwd,large.rot.volume)
  large.rot = aggregate(large.rot.volume ~ Plot, data=large, sum)
  large.mass = large.rot$large.rot.volume*cwd.params$Large[1]
  large.rot = cbind(large.rot, large.mass*0.5)   #Multiple here by 50% to account for rot
  colnames(large.rot) = c("Plot", "Volume", "Mass")
  
  #Repeating here for non-rotten large CWD
  
  large.snd.volume = rep(0,nrow(cwd))
  
  for(k in 1:nrow(cwd))
  {large.snd.volume[k]=pi^2*(sum(cwd$Large.sound[k]^2)/(8*12*20))}
  
  large = cbind(cwd,large.snd.volume)
  large.snd = aggregate(large.snd.volume ~ Plot, data=large, sum)
  large.snd.mass = large.snd$large.snd.volume*cwd.params$Large[1]
  large.snd = cbind(large.snd, large.snd.mass)
  colnames(large.snd) = c("Plot", "Volume", "Mass")
   
  #Combining here into a summary matrix of CWD mass by size class
  
  cwd.mass = cbind(fine.cts[1], round(fine.cts$fine.mass,2), round(small.cts$small.mass,2),
          round(medium.cts$medium.mass,2), round(large.rot$Mass,2), round(large.snd$Mass,2))
  cwd.mass = cbind(cwd.mass, rep(0,nrow(cwd.mass)))
  
  for(k in 1:nrow(cwd.mass))
  cwd.mass[k,7]=sum(cwd.mass[k,1:6])
  
  colnames(cwd.mass) = c("Plot", "Fine", "Small", "Medium", "Large.rotten", "Large.sound", "Total")
  
  rm(fine.mass,fine.volume,small.volume,medium.mass,medium.volume,large.rot.volume,large.snd.volume)
  
  cwd.mass = round(cwd.mass,1)
  
  #-------------------------------------------------------------------------------
  #Obtain estimate of biomass per hectare based on all plots
  
  attach(trees)
  
  tau.cap = N/sampled.area * sum(biomass)
  mean = tau.cap/(N/10000) #divide by number of hectares
  
  v.cap = N^2 * (1/sampled.area - 1/N) * var(biomass)
  se.tau.cap = sqrt(v.cap)
  
  t.val = qt(0.975, sampled.area-1)
  
  tree.CI = cbind(mean, (tau.cap - (se.tau.cap * t.val))/(N/10000), (tau.cap + (se.tau.cap * t.val))/(N/10000), 
                100*((se.tau.cap * t.val)/(N/10000))/mean)
  tree.CI[,1:3] = tree.CI[,1:3]/1000
  colnames(tree.CI) = c("Mean Biomass (Mg)", "Lower bound", "Upper bound", "MoE(%)")
  rownames(tree.CI) = c(as.character(levels(Site[1])))
  
  tree.CI
  #-------------------------------------------------------------------------------
  #Obtain estimate of biomass per hectare based on all plots
  
  attach(saps)
  
  tau.cap = N/sub.sampled.area * sum(biomass)
  mean = tau.cap/(N/10000) #divide by number of hectares
  
  v.cap = N^2 * (1/sub.sampled.area - 1/N) * var(biomass)
  se.tau.cap = sqrt(v.cap)
  
  t.val = qt(0.975, sub.sampled.area-1)
  
  sap.CI = cbind(mean, (tau.cap - (se.tau.cap * t.val))/(N/10000), 
                 (tau.cap + (se.tau.cap * t.val))/(N/10000), 100*((se.tau.cap * t.val)/(N/10000))/mean)
  sap.CI[,1:3] = sap.CI[,1:3]/1000
  colnames(sap.CI) = c("Mean Biomass (Mg)", "Lower bound", "Upper bound", "MoE(%)")
  rownames(sap.CI) = c(as.character(levels(Site[1])))
  
  sap.CI
  
  #-------------------------------------------------------------------------------
  #Classify each plot by species
  
  basal.area = rep(0, nrow(trees))
  
  attach(trees)
  
  for(k in 1:nrow(trees)) 
    basal.area[k] = 0.00007854*DBH.cm[k]^2 
  
  trees.2 = cbind(trees, basal.area)
  head(trees.2)
  
  #Cts = aggregate(basal.area ~ Species + Plot, data=trees, NROW)
  nt.per.plot = aggregate(DBH.cm ~ Plot, data=trees, NROW)
  nt.per.ha = nt.per.plot[,2] * (10000/(plot.size))
  
  plot.BAs = aggregate(basal.area ~ Plot, data=trees, sum)
  BAs.per.ha = plot.BAs[,2] * (10000/(plot.size))
  BAs.per.ha = round(BAs.per.ha,1)
  
  plot.ag.biomass = aggregate(ag.biomass ~ Plot, data=trees, sum)
  plot.bg.biomass = aggregate(bg.biomass ~ Plot, data=trees, sum)
  plot.biomass = aggregate(biomass ~ Plot, data=trees, sum)
  
  ag.biomass.per.ha = plot.ag.biomass[,2] * (10000/(plot.size))
  bg.biomass.per.ha = plot.bg.biomass[,2] * (10000/plot.size)
  biomass.per.ha = plot.biomass[,2] * (10000/plot.size)
  
  max.dbh = aggregate(DBH.cm ~ Species + Plot, data=trees, max)
  
  #-------------------------------------------------------------------------------
  #Calculate relative basal area for each of the plots
  
  sps.BAs = aggregate(basal.area ~ Species + Plot, data=trees, sum)
  attach(sps.BAs)
  
  rel.basal.area = rep(0, nrow(sps.BAs))
  plot.ba = rep(0, nrow(sps.BAs))
  sps.BAs = cbind(sps.BAs, plot.ba, rel.basal.area)
  
  for(k in 1:nrow(sps.BAs))
    for(j in 1:nrow(plot.BAs))
      if(sps.BAs$Plot[k] == plot.BAs$Plot[j])
      {sps.BAs$plot.ba[k] = plot.BAs$basal.area[j]
       sps.BAs$rel.basal.area[k] = 
         100*sps.BAs$basal.area[k]/sps.BAs$plot.ba[k]}
  
  #------------------------------------------------------------------------------
  #Calculate relative stem density
  
  Plot.stems = aggregate(basal.area ~ Plot, data=trees, NROW) #capital P here to denote different from plot.stems w/in sps.stems
  
  sps.stems = aggregate(basal.area ~ Species + Plot, data=trees, NROW)
  attach(sps.stems)
  
  rel.stems = rep(0, nrow(sps.BAs))
  plot.stems = rep(0, nrow(sps.BAs))
  sps.stems = cbind(sps.stems, plot.stems, rel.stems)
  
  colnames(sps.stems) = c("Species", "Plot", "stems", "plot.stems", "rel.stems")
  
  for(k in 1:nrow(sps.stems))
    for(j in 1:nrow(Plot.stems))
      if(sps.stems$Plot[k] == Plot.stems$Plot[j])
      {sps.stems$plot.stems[k] = Plot.stems$basal.area[j]
       sps.stems$rel.stems[k] = 
         100*sps.stems$stems[k]/sps.stems$plot.stems[k]}
  
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
  
