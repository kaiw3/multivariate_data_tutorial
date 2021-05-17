# Script by Kai Westwell
# 25th November 2020

# Set working directory to project directory

# Load packages
library(vegan)  # For performing the nmds
library(tidyverse)  # For data tidying functions and plotting


# Load data
barents <- read.csv("data/barents_data.csv")


# Separate environmental data and species data
barents_spp <- barents %>% select(Re_hi:Tr_spp)

barents_env_raw <- barents %>% select(ID_No:Temperature)


# Create ordinal groups for depth and temperature variables
barents_env <- barents_env_raw %>% 
  mutate(Depth_cat=cut(Depth, breaks=c(-Inf, 250, 350, Inf), labels=c("shallow","mid","deep"))) %>% 
  mutate(Temp_cat=cut(Temperature, breaks=c(-Inf, 1, 2, 3, Inf), labels=c("vc","c","m","w")))


# Perform nmds and fit environmental and species vectors
barents.mds <- metaMDS(barents_spp, distance = "bray", autotransform = FALSE)

barents.envfit <- envfit(barents.mds, barents_env, permutations = 999) # Fit environmental vectors

barents.sppfit <- envfit(barents.mds, barents_spp, permutations = 999) # Fit species vectors

barents.mds  # Stress value is less than 0.2, which is good. Shows how easy it was to condense multidimensional data into two dimensional space


# Save the results from the nmds and group the data by environmental variables
site.scrs <- as.data.frame(scores(barents.mds, display = "sites"))  # save NMDS results into dataframe

site.scrs <- cbind(site.scrs, Depth = barents_env$Depth_cat)  # make depth a grouping variable and save to dataframe

site.scrs <- cbind(site.scrs, Temperature = barents_env$Temp_cat) # make temperature a grouping variable and save to dataframe

head(site.scrs)  # View the dataframe


# Rename Environmental Factor levels
site.scrs$Temperature <- recode(site.scrs$Temperature, vc = "Very Cold", c = "Cold", m = "Medium", w = "Warm")

site.scrs$Depth <- recode(site.scrs$Depth, shallow = "Shallow", mid = "Mid", deep = "Deep")


# Save species data from nmds to dataframe
spp.scrs <- as.data.frame(scores(barents.sppfit, display = "vectors"))  # Save species values into dataframe

spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))  # Add species names to dataframe

spp.scrs <- cbind(spp.scrs, pval = barents.sppfit$vectors$pvals) # Add p values to dataframe so you can select species which are significant

sig.spp.scrs <- subset(spp.scrs, pval<=0.05) # Show only significant species (using 0.05 as cutoff)

head(spp.scrs)  # View dataframe


# Save environmental extrinsic variables
env.scores.barents <- as.data.frame(scores(barents.envfit, display = "vectors"))  # Extract nmds scores of all environmental variables from envifit dataframe

env.scores.barents <- cbind(env.scores.barents, env.variables = rownames(env.scores.barents))  # Name them 

env.scores.barents <- cbind(env.scores.barents, pval = barents.envfit$vectors$pvals) # Add p values to dataframe

sig.env.scrs <- subset(env.scores.barents, pval<=0.05) # Show only significant variables (using 0.05 as cutoff)

head(env.scores.barents)  # View dataframe


# Basic nmds
(nmds.plot.barents <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2)) + # Create the plot
  geom_point(aes(x = NMDS1, y = NMDS2, colour = factor(site.scrs$Depth), shape = factor(site.scrs$Temperature)), size = 2)+ # Add site points to plot with the shape determined by temperature and colour determined by depth
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid")) +
  labs(colour = "Depth", shape = "Temperature", title = "Does Fish Species Composition Change at\n Varying Water Depths and Temperatures") + # Add legend titles
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))) # Add legend

ggsave(filename = "outputs/nmds_1.png", nmds.plot.barents, height = 7, width = 14)


# Add species vector arrows
(nmds.plot.barents.2 <- nmds.plot.barents +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # Add vector arrows for significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25) # Add labels for species (with ggrepel::geom_text_repel so that labels do not overlap)
)

ggsave(filename = "outputs/nmds_2.png", nmds.plot.barents.2, height = 7, width = 14)


# Add environmental variable vector arrows
(nmds.plot.barents.3 <- nmds.plot.barents +
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # Add vector arrows of significant environmental variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25) # Add labels
)

ggsave(filename = "outputs/nmds_3.png", nmds.plot.barents.3, height = 7, width = 14)

# We can see that theres a distinction in samples based on depth and temperature. 
# Species info shows some species are more commonly found in different conditions.
# None of this is statistically verified yet, so a ANOSIM, SIMPER or indicator species analysis will need to be run to check.


# ANOSIM
mat_bar_spp <- as.matrix(barents_spp)


# ANOSIM with depth grouping
bar_depth <- anosim(mat_bar_spp, barents_env_raw$Depth, distance = "bray", permutations = 9999)

bar_depth  # ANOSIM significance is less than 0.05 so it is significant. R statistic is relatively low, suggesting the groups are similar to each other.


# ANOSIM with temperature grouping
bar_temp <- anosim(mat_bar_spp, barents_env_raw$Temperature, distance = "bray", permutations = 9999)

bar_temp  # Significance means this looks good too, and relatively low R statistic suggests similarity between groups.


# Challenge data
barents_env_chal <- barents_env_raw %>% 
  mutate(lat_cat=cut(Latitude, breaks=c(-Inf, 72, 74, Inf), labels=c("south","centre","north"))) %>% 
  mutate(long_cat=cut(Longitude, breaks=c(-Inf, 23, 30, Inf), labels=c("west","centre","east")))

