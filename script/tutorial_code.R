# Script by Kai Westwell
# 25th November 2020

# Set working directory to project directory

# Load packages
library(vegan)  # For performing the nmds
library(gplots)  # For heatmaps
library(tidyverse)  # For data tidying functions and plotting

# Load data
barents <- read.csv("data/barents_data.csv")

# Separate environmental data and species data
barents_spp <- barents %>% select(Re_hi:Tr_spp)
barents_env_raw <- barents %>% select(ID.No:Temperature)

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
(nmds.plot.barents <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ # Create the plot
  geom_point(aes(x = NMDS1, y = NMDS2, colour = factor(site.scrs$Depth), shape = factor(site.scrs$Temperature)), size = 2)+ # Add site points to plot with the shape determined by temperature and colour determined by depth
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Depth", shape = "Temperature", title = "Does Fish Species Composition Change at\n Varying Water Depths and Temperatures")+ # Add legend labels
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # Add legend
)
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








barents_pca$x
barents_pca <- prcomp(barents_env_raw)
summary(barents_pca)
screeplot(barents_pca, type = "lines")
(barents_pca$sdev)^2

y <- matrix(rnorm(89)) 
c <- cor(barents_spp, method="spearman")
d <- as.dist(1-c)
pear.dist <- as.dist((1-c)/2)
levelplot(c)
hp <- hclust(pear.dist, method="ward.D")
plot(hp, hang=-1, cex=0.2)

data.med.center <- sweep(barents_spp, MARGIN=1, STATS=apply(barents_spp,1,median), FUN = "-",) 
data.med.center[ data.med.center > 5 ] <- 5
data.med.center[ data.med.center < -5 ] <- -5

pear.cols <- cor(t(data.med.center), method="pearson")
pear.dist <- as.dist((1-pear.cols)/2)
hc.cols <- hclust(pear.dist, method="ward.D")
pear.rows <- cor(data.med.center, method="pearson")
pear.dist <- as.dist((1-pear.rows)/2)
hc.rows <- hclust(pear.dist, method="ward.D")

heatmap.2(t(data.med.center), 
          scale="none", 
          cluster.by.col=TRUE, 
          cluster.by.row=TRUE, 
          Rowv=TRUE,
          Colv=TRUE,
          hclust.col = hc.cols, 
          hclust.row = hc.rows, 
          labCol=NULL,
          cexRow = 0.5,
          )


# Failed NMDS
community_matrix <- matrix(
  sample(1:100,300,replace=T),nrow=10,
  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))

example_NMDS <- metaMDS(barents_spp, k = 2)
stressplot(example_NMDS)        # Produce shepherds diagram to test goodness of fit

plot(example_NMDS, "sites")   # Produces distance 
orditorp(example_NMDS, "sites")
colvec <- c("gray0", "gray0", "gray49", "gray49")   # Identifies colors for group assignments
pchvec <- c(21, 21, 22, 22)

ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)







