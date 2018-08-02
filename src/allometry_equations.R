# Allometric equation functions, labelled by <species_code>_<aboveground/belowground>_<authoryear>
# Example: total aboveground biomass equation for Rhizophora apiculata by Ong et al., 2004:
#   - rhap_ag_ong2004


# General functions, Komiyama et al., 2005, function of diameter at breast height and wood density

general_ag_komiyama2005 <- function(dbh, dens){
  
  0.251*dens*(dbh^2.46)
  
}

general_bg_komiyama2005 <- function(dbh, dens){
  
  0.199*(dens^0.899)*(dbh^2.22)
  
}

# Avicennia marina, Comley et al., 2005 (Australia)

avma_ag_comley2005 <- function(dbh, dens) {
  
  10^(-0.511*2.113*log10(dbh))
  
}

avma_bg_comley2005 <- function(dbh, dens) {
  
  10^(0.106*1.171*log10(dbh))
  
}


# Rhizophora apiculata, Ong et al., 2004 (Malaysia)

rhap_ag_ong2004 <- function(dbh, dens) {
  
  if(dbh >= 15) {
    exp(1)^(2.420*log(dbh) - 1.832)
  } else {
    exp(1)^(2.318*log(dbh) - 1.671)
  }
  
}

rhap_bg_ong2004 <- function(dbh, dens) {
  ## multiple equations for total belowground biomass dependent upon DBH
  
  if(dbh >= 15) {
    exp(1)^(2.611*log(dbh) - 3.454)
  } else {
    exp(1)^(1.522*log(dbh) - 1.707)
  }

}

rhap_prop_ong2004 <- function(dbh, dens) {

  if(dbh >= 15) {
    exp(1)^(2.546*log(dbh) - 2.94)
  }

}


# Rhizophora stylosa, Comley et al., 2005 (Australia)

rhst_ag_comley2005 <- function(dbh, dens) {
  
  10^(-0.696*2.465*log10(dbh))
  
}

rhst_bg_comley2005 <- function(dbh, dens) {
  
  10^(-0.583*1.860*log10(dbh))
  
}

