<center><img src="{{ site.baseurl }}/photos/tut_banner.png" alt="Img"></center>

Created by Kai Westwell

### Tutorial Aims

#### <a href="#section1"> 1. Understand the basics of multivariate analyses\n 
<a href="#section1a"> - Multivariate Stats Overview</a>
<a href="#section1b"> - NMDS</a>
<a href="#section1c"> - ANOSIM</a>

#### <a href="#section2"> 2. Make an NMDS plot</a>

#### <a href="#section3"> 3. Perform a SIMPER analysis</a>

This tutorial will cover some of the basic methods in analysing multivariate data. You will learn what multivariate data is, some of the ways you can analyse these kinds of data, and then you can practice it for yourself!
---------------------------

## Introduction

Ecological data can be complex. Often, large numbers of variables need to be studied in order to obtain an accurate picture of the system in question. This complexity can make the data hard to interpret, and even harder to analyse statistically. This tutorial will take you through the basics of understanding what multivariate data look like, and will introduce you to some of the ways we can use statistics to interpret these data. In order to transition into dealing with the wide array of ecological data you are likely to be presented with when looking at real systems, its important to be able to deal with complex data. And hopefully, as you take these skill forward, you will find more and more ways that to can apply multivariate stats into solving a problem!

You can get all of the resources for this tutorial from <a href="https://github.com/kaiw3/multivariate_data_tutorial" target="_blank">this GitHub repository</a>. Clone and download the repo as a zip file, then unzip it.

<a name="section1"></a>

## 1. Understand the basics of multivariate analyses

<a name="section1a"></a>

### Multivariate Stats Overview

So what do we mean when we talk about multivariate data? Well its pretty much what it sounds like. Multivariate data are data that contain more than 2 response variables (although theres usually quite a few more than just 3). Multivariate statistics are where you analyse these multiple variables simultaneously. In biology these data usually come about, either from counts of species abundances in assemblages, where each species acts as a variable; or from physical properties of environments (such as temperature, pH or habitat structure).

You will already be familiar with bivariate statistical tests (where ther are 2 response variables), such as t tests and one and two-way ANOVAs. So why look at multiple variables together when you can look at them separately? When you look at 2 variables in isolation, you are ignoring the many possible interactions between these variables, and the rest of the variables that you measured in the study. Multivariate statistics allow us to look at how multiple variables change together, and avoid the confounding effects that could occur from running the analyses separately. While looking at individual variables can sometimes be enough to answer your research question, if you are aiming to look at the overall effects of all variables, multivariate stats is the way to go. This can allow us to reveal patterns that would not be found by examining each variable separately. 

<a name="section1b"></a>

### NMDS




<a name="section2"></a>

## 2. Make an NMDS Plot

First, lets open R studio, make a new script and load the data and packages.

```
# Set working directory

# Load packages
library(vegan)  # For performing the nmds
library(tidyverse)  # For data tidying functions and plotting

# Load data
barents <- read.csv("data/barents_data.csv")
```
Now lets make the data easier to work with. Separate the species and environmental data and separate the environmental variables that we are interested in into ordinal groups.

```
# Separate environmental data and species data
barents_spp <- barents %>% select(Re_hi:Tr_spp)
barents_env_raw <- barents %>% select(ID.No:Temperature)

# Create ordinal groups for depth and temperature variables
barents_env <- barents_env_raw %>% 
  mutate(Depth_cat=cut(Depth, breaks=c(-Inf, 250, 350, Inf), labels=c("shallow","mid","deep"))) %>%
  mutate(Temp_cat=cut(Temperature, breaks=c(-Inf, 1, 2, 3, Inf), labels=c("vc","c","m","w")))
```

Now we have arranged the data, we can perform the nmds.

```
# Perform nmds and fit environmental and species vectors
barents.mds <- metaMDS(barents_spp, distance = "bray", autotransform = FALSE)
barents.envfit <- envfit(barents.mds, barents_env, permutations = 999) # Fit environmental vectors
barents.sppfit <- envfit(barents.mds, barents_spp, permutations = 999) # Fit species vectors
```
Have a look at the nmds output and check the stress. Sometimes he nmds cant represent all of the relationships between variables accurately. This is reflected by a high stress value. A general rule is if the stress value is below 0.2, the plot is generally ok.
```
barents.mds  # Stress value is less than 0.2, which is good. Shows how easy it was to condense multidimensional data into two dimensional space
```

After you have performed the nmds, you need to save the outputs so we can graph it later. Here you will also group the data by the environmental variables we are interested in : depth and temperature.
```
# Save the results from the nmds and group the data by environmental variables
site.scrs <- as.data.frame(scores(barents.mds, display = "sites"))  # save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Depth = barents_env$Depth_cat)  # make depth a grouping variable and save to dataframe
site.scrs <- cbind(site.scrs, Temperature = barents_env$Temp_cat) # make temperature a grouping variable and save to dataframe
head(site.scrs)  # View the dataframe
```

Take a look at the output you have saved to become familiar with the structure for future use of this dataframe.
```
head(spp.scrs)  # View dataframe
```

Now save the species specific data from the nmds analysis to a dataframe, so that we can plot this separately on the nmds plot later. Examine the dataframe so that you are comfortable with using it in later stages.
```
# Save species data from nmds to dataframe
spp.scrs <- as.data.frame(scores(barents.sppfit, display = "vectors"))  # Save species values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))  # Add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = barents.sppfit$vectors$pvals) # Add p values to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) # Show only significant species (using 0.05 as cutoff)
head(spp.scrs)  # View dataframe
```

Do the same thing with the environmental variables. This will let you plot the vector information of either the species or environmental groups, or both, on top of the nmds plot.
```
# Save environmental variables
env.scores.barents <- as.data.frame(scores(barents.envfit, display = "vectors"))  # Extract nmds scores of all environmental variables from envifit dataframe
env.scores.barents <- cbind(env.scores.barents, env.variables = rownames(env.scores.barents))  # Name them 
env.scores.barents <- cbind(env.scores.barents, pval = barents.envfit$vectors$pvals) # Add p values to dataframe
sig.env.scrs <- subset(env.scores.barents, pval<=0.05) # Show only significant variables (using 0.05 as cutoff)
head(env.scores.barents)  # View dataframe
```

Now we can get to the fun part, plotting our nmds data! We'll use ggplot2 which you should already be familiar with from <a href="https://ourcodingclub.github.io/tutorials/datavis/" target="_blank">previous coding club tutorials</a>. Here we will use all of the same ggplot2 functions you have already learned about, so this section should be relatively straightforward.
```
# Basic nmds
(nmds.plot.barents <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ # Create the plot
  geom_point(aes(x = NMDS1, y = NMDS2, colour = factor(site.scrs$Depth), shape = factor(site.scrs$Temperature)), size = 2)+ # Add site points to plot with the shape determined by temperature and colour determined by depth
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Depth", shape = "Temperature", title = "Does Fish Species Composition Change at\n Varying Water Depths and Temperatures")+ # Add legend labels
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # Add legend
)
```

<center> <img src="{{ site.baseurl }}/photos/nmds_1.png" alt="nmds" style="width: 800px;"/> </center>

And here we have an nmds plot! We can see that there are some different groupings going on here, with some samples being found in warmer temperatures or greater depths, for example. But we don't know whic species these groups relate to...it's luck we saved the species data from the nmds then!\n

Lets add an overlay with species vectors.
```
# Add species vector arrows
(nmds.plot.barents.2 <- nmds.plot.barents +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # Add vector arrows for significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25) # Add labels for species (with ggrepel::geom_text_repel so that labels do not overlap)
)
```

<center> <img src="{{ site.baseurl }}/photos/nmds_2.png" alt="nmds" style="width: 800px;"/> </center>

Great! Now we can see certain species group more in warmer water, or in colder water. We can also see how strong these relationships are based on the length of the arrows. While we can see that some environmental groupings exist, we may want to get a more clear idea of the directions these are acting in by overlaying the environmental nmds data that we also saved earlier. Lets include all of the measured variables just so that we can see what datasets with lots of variables would look like.
```
# Add environmental variable vector arrows
(nmds.plot.barents.3 <- nmds.plot.barents +
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + # Add vector arrows of significant environmental variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25) # Add labels
)
```

<center> <img src="{{ site.baseurl }}/photos/nmds_3.png" alt="nmds" style="width: 800px;"/> </center>

<a name="section3"></a>

## 3. Perform a SIMPER analysis

More text, code and images.

<a name="section4"></a>

## 4. Challenge

This is the end of the tutorial. Summarise what the student has learned, possibly even with a list of learning outcomes. In this tutorial we learned:

##### - how to generate fake bivariate data
##### - how to create a scatterplot in ggplot2
##### - some of the different plot methods in ggplot2


