# Accuracy assessment

library(janitor)
library(raster)
library(rgdal)
library(sf)
library(tidyverse)

#---------------------------------

source("./src/helper_funcs.R")

in_dir <- "/home/jbukoski/research/data/thailand_stocks/input/"
out_dir <- "/home/jbukoski/research/data/thailand_stocks/output/"

site <- "nakorn_"
year <- "2017_"

#---------------------------------

acc <- readOGR(dsn = paste0(out_dir, site, year, "valid_pts/"), layer = "valid_pts")
rast <- raster(paste0(out_dir, site, year, "svm.tif"))

mat <- table(acc$class, acc$ref_class) %>%
  matrix(ncol = 3) %>%
  as_tibble

areas <- aggregate(getValues(area(rast, weights = FALSE)), 
                   by = list(getValues(rast)), sum)

weights <- areas %>%
  mutate(ha = x / 10000,
         total = sum(ha),
         wgt = ha / total) %>%
  rbind(colSums(.))

N <- weights$ha[nrow(weights)]

new_mat <- mat %>% mutate(Total = rowSums(.)) %>%
  rbind(colSums(.)) %>%
  cbind(weights$ha, round(weights$wgt, 4)) %>%
  rename(ha = "weights$ha",
         wgt = "round(weights$wgt, 4)")

my_func <- function(val, sum, wgt) {
  z <- wgt * (val / sum)
  return(z)
}

test <- new_mat %>%
  mutate(v1 = my_func(V1, 200, wgt),
         v2 = my_func(V2, 200, wgt),
         v3 = my_func(V3, 200, wgt)) %>%
  select(v1, v2, v3) %>%
  slice(-nrow(.)) %>%
  cbind(rowSums(.)) %>%
  rbind(colSums(.)) %>%
  rename(total = "rowSums(.)")

se_func <- function(val, sum, wgt) {
  z <- wgt^2 * (val/sum) * (1 - (val/sum)) / (sum - 1)
  
  return(z)
}

agri_se <- sqrt(sum(se_func(new_mat$V1, new_mat$Total, new_mat$wgt)[1:3]))
aqua_se <- sqrt(sum(se_func(new_mat$V2, new_mat$Total, new_mat$wgt)[1:3]))
mang_se <- sqrt(sum(se_func(new_mat$V3, new_mat$Total, new_mat$wgt)[1:3]))

#-----------------------------
# Calculate User's, Producer's and Overall accuracy
# Formulas from page 301-302 in Warner et al, 2009 (Sage Handbook of RS)

diag_p <- diag(as.matrix(test[1:3, 1:3]))
diag_n <- diag(as.matrix(new_mat[1:3, 1:3]))

overall_var <- function(N_i, N, U_i, n_i) {
  z <- (N_i^2 / N^2) * (U_i * (1-U_i) / (n_i - 1))
  return(z)
}

marg_ttl <- function(nij, Ni, ni) {
  z = nij * ( Ni / ni)
  return(z)
}

calc_prod <- function(i) {
  val <- (diag_n[i] * areas$x[i] / new_mat$Total[i]) / sum(marg_ttl(new_mat[1:3, i], areas$x[1:3], new_mat$Total[1:3]))
  return(val)
}


users_var <- function(Ui, ni) {
  z <- Ui * (1 - Ui) / (ni - 1)
  return(z)
}

prods_var <- function(Nj, Pj, Uj, nj, nij) {
  term_a <- 1 / (marg_ttl(nij, Nj, nj))
  
}

calc_prods_se <- function(j) {
  
  c <- c(1, 2, 3)
  idx <- j != c
  i <- c[idx]
  
  for(j in j) {
    
    term_a <- areas$x[j]^2 * ((1 - prods[j])^2) * users[j] * (1 - users[j]) / (new_mat$Total[j] - 1)
    term_b <- c()
    
    for(i in i) {
      term_b_add <- areas$x[i]^2 * (new_mat[i, j] / 200) * (1 - new_mat[i, j]/200) / (200-1)
      term_b <- c(term_b, term_b_add)
    }
    
    numerator <- term_a + prods[j]^2 * (sum(term_b))
    z <- marg_ttl(new_mat[1:3, j], areas$x[1:3], new_mat$Total[1:3])
    
  }
  
  denom <- sum(z)
  se <- sqrt( (1 / denom^2) * numerator)
  
  return(se)
}

#--------------------------
# Apply functions to calculate accuracy and accuracy uncertainties

overall <- sum(diag_p)
users <- round(diag_p / test$total[1:3], 4)
prods <- round(c(calc_prod(1), calc_prod(2), calc_prod(3)), 4)

overall_se <- sqrt(sum(overall_var(new_mat$ha[1:3], N, users, 200)))
users_se <- round(sqrt((users_var(users, new_mat$Total[1:3]))), 4)
prods_se <- round(c(calc_prods_se(1), calc_prods_se(2), calc_prods_se(3)), 4)
