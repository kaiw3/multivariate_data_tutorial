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




You can surround package names, functions, actions ("File/ New...") and small chunks of code with backticks, which defines them as inline code blocks and makes them stand out among the text, e.g. `ggplot2`.

When you have a larger chunk of code, you can paste the whole code in the `Markdown` document and add three backticks on the line before the code chunks starts and on the line after the code chunks ends. After the three backticks that go before your code chunk starts, you can specify in which language the code is written, in our case `R`.



At this point it would be a good idea to include an image of what the plot is meant to look like so students can check they've done it right. Replace `IMAGE_NAME.png` with your own image file:

<center> <img src="{{ site.baseurl }}/IMAGE_NAME.png" alt="Img" style="width: 800px;"/> </center>

<a name="section1"></a>

## 3. The third section

More text, code and images.

This is the end of the tutorial. Summarise what the student has learned, possibly even with a list of learning outcomes. In this tutorial we learned:

##### - how to generate fake bivariate data
##### - how to create a scatterplot in ggplot2
##### - some of the different plot methods in ggplot2

We can also provide some useful links, include a contact form and a way to send feedback.

For more on `ggplot2`, read the official <a href="https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf" target="_blank">ggplot2 cheatsheet</a>.

Everything below this is footer material - text and links that appears at the end of all of your tutorials.

<hr>
<hr>

#### Check out our <a href="https://ourcodingclub.github.io/links/" target="_blank">Useful links</a> page where you can find loads of guides and cheatsheets.

#### If you have any questions about completing this tutorial, please contact us on ourcodingclub@gmail.com

#### <a href="INSERT_SURVEY_LINK" target="_blank">We would love to hear your feedback on the tutorial, whether you did it in the classroom or online!</a>

<ul class="social-icons">
	<li>
		<h3>
			<a href="https://twitter.com/our_codingclub" target="_blank">&nbsp;Follow our coding adventures on Twitter! <i class="fa fa-twitter"></i></a>
		</h3>
	</li>
</ul>

### &nbsp;&nbsp;Subscribe to our mailing list:
<div class="container">
	<div class="block">
        <!-- subscribe form start -->
		<div class="form-group">
			<form action="https://getsimpleform.com/messages?form_api_token=de1ba2f2f947822946fb6e835437ec78" method="post">
			<div class="form-group">
				<input type='text' class="form-control" name='Email' placeholder="Email" required/>
			</div>
			<div>
                        	<button class="btn btn-default" type='submit'>Subscribe</button>
                    	</div>
                	</form>
		</div>
	</div>
</div>
