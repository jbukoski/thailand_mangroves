# R file to generate plots

#-----------------------------------------------------------------------------
# Generate plots

soil_plot %>%
  filter(plot <= 7) %>%
  ggplot(aes(x = plot, y = plot_soc, col = site)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = min_soc, 
                    ymax = max_soc), width = 0.15,
                size = 1.1) + 
  labs(y = "Soil Organic Carbon (Mg C/ha)",
       x = "Plot") +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  ylim(0, 1500) +
  scale_colour_discrete(name="Site", 
                        labels = c("Krabi", "Nakorn"))


soil_site %>% 
  ggplot(aes(x = site, y = soc)) + 
  geom_errorbar(aes(ymin = soc - soc_se, 
                    ymax = soc + soc_se), width = 0.05) + 
  geom_point() + 
  ggtitle("Soil Carbon Plot Means +/- 1 St. Error") +
  ylab("Soil Organic Carbon (Mg/ha)") +
  xlab("Site") +
  theme_bw() +
  ylim(0, 1000)

#------------------------------------------------------------------------------
# Add other relevant parameters to characterize the plot (i.e. mean DBH and Tons ag.biomass)

error_summary <- trees_ci %>%
  left_join(soil_site, by = c("site" = "site")) %>%
  left_join(site_cwd, by = c("site" = "site")) %>%
  dplyr::select(site, ag_se, bg_se, soc_se, cwd_se) %>%
  mutate(ag_se = ag_se,
         bg_se = bg_se,
         soc_se = soc_se,
         cwd_se = cwd_se) %>%
  gather(pool, error, -site) %>%
  arrange(site) %>%
  mutate(pool = rep(c("agc", "bgc", "soc", "cwd"), 2))

summary <- trees_ci %>%
  dplyr::select(site, mean_biomass, ag_mean, bg_mean) %>%
  left_join(soil_site, by = c("site" = "site")) %>%
  left_join(site_cwd, by = c("site" = "site")) %>%
  mutate(bgc = bg_mean *.46 * -1,
         agc = ag_mean *.46,
         soc = soc * -1,
         cwd = cwd * 0.5) %>%
  dplyr::select(site, agc, bgc, soc, cwd) %>%
  gather(pool, value, -site) %>%
  left_join(error_summary, by = c("site", "pool")) %>%
  mutate(error2 = error)

##---------------------------------------------------------------------------------
## Visualizing the data

# Ordering of stacked bars is determined by levels of the factor

test_dat <- tibble(
  site = c(rep("Krabi", 4), rep("Nakorn", 4)),
  pool = rep(c(rep("Above", 2), rep("Below", 2)), 2),
  min_err = c(72.3, NA, -888.9 - 53.41, NA, 108 + 8.10, NA, -39.9 - 251 - 18 - 8.17, NA),
  max_err = c(72.3  + 18.88, NA, -888.9, NA, 108 + 8.10 + 23.2 + 2.85, NA, -39.9 - 251, NA)
)

p1 <- summary %>% 
  mutate(pool = factor(pool, levels = c("agc", "cwd", "soc", "bgc"))) %>%
  ggplot(aes(x = site, y = value, fill = pool)) +
  geom_bar(stat = "identity", width = 0.3) +
  geom_bar(stat = "identity", width = 0.3, color="black", show.legend=FALSE)+
  theme_bw() +
  xlab("Site") +
  ylab("Carbon storage (Mg C/ha)") +
  ylim(-1000, 200) +
  scale_fill_discrete(name="Ecosystem C Pools", 
                      breaks = c("agc", "cwd", "bgc", "soc"),
                      labels = c("Aboveground Biomass", 
                                 "Coarse Woody Debris",
                                 "Belowground Biomass",
                                 "Soil Organic Carbon")) +
  theme(text = element_text(size = 22))

p1 +  
  geom_linerange(aes(x = test_dat$site,
                     ymin = test_dat$min_err, 
                     ymax = test_dat$max_err)) +
  geom_segment(aes(x = 0.95, xend = 1.05, y = -888.9 - 53.41, yend = -888.9 - 53.41)) +
  geom_segment(aes(x = 0.95, xend = 1.05, y = 72.3  + 18.88, yend = 72.3  + 18.88)) +
  geom_segment(aes(x = 1.95, xend = 2.05, 
                   y = -39.9 - 251 - 18 - 8.17, 
                   yend = -39.9 - 251 - 18 - 8.17)) +
  geom_segment(aes(x = 1.95, xend = 2.05, 
                   y = 108 + 8.10 + 23.2 + 2.85, 
                   yend = 108 + 8.10 + 23.2 + 2.85))

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
  dplyr::select(site, plot, subplot, agb_tot, agb_ha, bgb_tot, bgb_ha, soc, tot_c, distance, n_sps) %>%
  distinct()

trees %>%
  filter(plot <= 7) %>%
  group_by(site, plot) %>%
  summarize(n = n(),
            agb_se = sd(agb)/sqrt(n) * 0.43,
            agb = 10000*sum(agb)/plot_size/1000 * 0.43,
            bgb_se = sd(bgb)/sqrt(n) * 0.43,
            bgb = 10000*sum(bgb)/plot_size/1000 * 0.43) %>%
  ggplot(aes(x = plot, y = agb + bgb, col = site)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = agb + bgb - agb_se - bgb_se, 
                    ymax = agb + bgb + agb_se + bgb_se), width = 0.15,
                size = 1.1) + 
  labs(y = "Biomass Carbon (Mg C/ha)",
       x = "Plot") +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  scale_colour_discrete(name="Site", 
                        labels = c("Krabi", "Nakorn"))