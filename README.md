Repository for code and data analysis associated with Thailand carbon stock estimation in Krabi Mangrove Forest and the Pak Phanang Mangrove Forest

Author: Jacob J. Bukoski

Research questions:

  1. Is dominant species a significant predictor of carbon stocks in mangroves?

Goals:

  1. Quantify differences in carbon stocks, species composition and forest structure between Pak Phanang and Krabi River Estuary mangrove forests
  2. Test for significant effect of species composition on aboveground and belowground carbon stocks
  3. Test for control of distance on carbon stocks in mangroves

Data analysis plan:

  1. Quantify component specific carbon stocks in Pak Phanang and Krabi River Estuary (site level)
  2. Subplot-level measures of forest structure, biomass and carbon stocks -- test for significant effects of species composition, forest structure, and distance from water

Ancillary goals:

  1. R package for mangrove allometry
  2. Publish tree-level dataset in Data Dryad

#### Overview of scripts in /src/

Remote sensing analysis:

  1. GEE Code for KRE and PPM -> exports relevant Landsat imagery & SRTM data clipped to each site.
    - [Krabi](https://code.earthengine.google.com/0e9bdf9b5ab1f84b83c2ab5ac6588ebc) (21 Sep. 2018)
    - [Nakorn](https://code.earthengine.google.com/1ce04f7129e3f9f29518b232e24d7c54) (21 Sep. 2018)
  2. `lsat_hist_matching.R` - histogram matches the three preceding decadal time dates to the 2017 imagery.
  3. `lsat_preprocessing.R` - computes Landsat transformations, stacks as a raster, and exports as a finalized "processed" Landsat.
  4. `build_svm.R` - trains, validates, and outputs the SVM classifications of the Landsat imagery for each site.
  
Field data analysis:

  5. `site_analaysis.R` - processes field data from the mangrove forests and abandoned aquaculture ponds

Visualization:

  6. `build_site_map.R` - builds Figure 1, the site map for KRE and PPM
  7. `build_field_plots.R` - builds figures corresponding to visualization of forest structure and carbon stock field data
  8. `build_luc_maps.R` - builds figures corresponding to remote sensing analysis of land use change
  
Other:

  9. `allometry.R` - helper script with allometry functions
  10. `allometry_equations.R` - helper script with actually allometric equations for mangroves
  ~11. `build_water_mask.R` - helper script to build a water mask (may not actually use?)~
  12. `helper_funcs.R`- generic helper functions for the analysis
  13. `unsupervised.R` - quick examinition of unsupervised classification of the RS imagery
  14. `utm_to_latlong.R` - helper function to convert plot lat/long values to decimal degree
