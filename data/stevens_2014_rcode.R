####################################################
## Data S2: R code
####################################################

###################################################
### stevens_rcode.R
### Created by Jeffrey R. Stevens on 3 Feb 2014 (jeffrey.r.stevens@gmail.com),
### Summary: This script calculates statistics and generates figures for the
###      analysis of evolutionary hypotheses explaining comparative intertemporal choice.
### Instructions: Place this file and the data files (stevens_data.csv)
### 	in the same directory.  Create a folder called "figures". Set the R
### 	working directory to this directory.  At the R command prompt, type 
### 	> source("stevens_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PDF versions of the figures in the figures directory.
### Uses: This script can be reproduced and modified for personal and scientific use.
###		Please cite the original source.
### Data files: Description of the data columns:
###  stevens_dataS1--intertemporal choice and comparative data
###   subject - name of subject if available
###   species - species (common name) of subject (e.g., Marmoset)
###   type - specific common name of subject (e.g., Common)
###   latin_name - Latin name of species
###   measure - type of data (indiff = indifference point, bodywt = body weight,
###		ecv = endocranial volume, groupsize = group size, homerange = home range,
###		lifespan = lifespan)
###   value - value for specific measure (units: indifference point = seconds,
###		body weight = grams, endocranial volume = cubic centimeters, group size =
###		individuals, home range = hectares, lifespan = years)
###   n - number of samples
###   sd - standard deviation of samples
###   reference - citation for source of data
###   notes - miscellaneous notes on the measure, value, or source
###################################################

 
######################
## Clear variables, load libraries, load data
######################
rm(list=ls())			# clears all variables
library(beeswarm) 		# needed for beeswarm plot
library(car)  			# needed for recode
library(caper)			# needed for pgls
library(foreach)		# needed for foreach
library(lattice)  		# needed for plots
library(latticeExtra)	# needed for axis scaling
library(psych)  		# needed for principal() for PCR
library(xtable)			# needed for xtable

## Load and prepare data
raw_data <- read.csv("stevens_data.csv")							# read in comparative data
raw_data$n <- recode(raw_data$n, "'?' = 1")							# if sample size (n) is unknown, recode it as 1
raw_data$n <- as.numeric(levels(raw_data$n)[as.integer(raw_data$n)])# convert samples sizes to numeric

species <- levels(raw_data$latin_name)		# extract species names
species_nospace <- gsub(" ", "_", species)	# remove space from species names
num_species <- length(species)				# calculate number of species
measures <- levels(raw_data$measure)		# extract measures
num_measures <- length(measures)			# calculate number of measures
line_counter <- 1							# initialize line counter
mean_data <- data.frame(matrix(rep(NA, num_species * num_measures * 3), ncol = 3))	# initialize mean_data data frame
names(mean_data) <- c("species", "measure", "value")	# name mean_data

######################
## Calculate weighted mean of different measures for each species
######################
foreach(j = measures) %do% {									# for each measure
  measures_subset <- subset(raw_data, measure == j)				# create a subset of the data for that measure
  foreach(k = species) %do% {									# for each species (within each measure)
    species_subset <- subset(measures_subset, latin_name == k)	# create a subset of the data for that measure
    weighted_mean <- sum(species_subset$value * species_subset$n, na.rm = T) / sum(species_subset$n, na.rm = T)	# calculate a weighted mean for that measure and species
    mean_data[line_counter, 1] <- k								# assign species to column 1
    mean_data[line_counter, 2] <- j								# assign measure to column 2
    mean_data[line_counter, 3] <- weighted_mean					# assign weighted mean to column 3
    line_counter <- line_counter + 1							# increment line counter
  }
}

######################
## Create primate data frame
######################
primate_data <- unstack(mean_data, form = "value ~ measure")					# unstack the data to create columns for each measure
primate_data <- data.frame(species = species_nospace, primate_data)				# create data frame with species names and measure data
primate_data$logindiff <- log(primate_data$indiff)								# create log indifference point data
primate_data$loglifespan <- log(primate_data$lifespan)							# create log lifespan data
primate_data$logbodywt <- log(primate_data$bodywt)								# create log body weight data
primate_data$logECV <- log(primate_data$ecv)									# create log brain size data
primate_data$loggroupsize <- log(primate_data$groupsize)						# create log group size data
primate_data$loghomerange <- log(primate_data$homerange)						# create log home range size data
primate_table <- data.frame(species, primate_data[, c(6, 2, 3, 7, 5, 4)])  	# create data frame for table
names(primate_table) <- c("Species", "Waiting time", "Body weight", "Brain size", "Lifespan", "Home range", "Group size")

######################
## Create phylogenetic tree
######################
primate_tree_nex <- "((((Macaca_fascicularis:5.045978,Macaca_mulatta:5.045978):24.954022,((Gorilla_gorilla:8.652233,(Pan_paniscus:2.333553,Pan_troglodytes:2.333553):6.318679):6.480222,Pongo_pygmaeus:15.132455):14.867546):16.811821,(Ateles_geoffroyi:21.321300,((Callithrix_jacchus:15.712251,Saguinus_oedipus:15.712253):4.380275,Sapajus_apella:20.092528):1.228773):25.490521):26.191197,(Eulemur_macaco:20.542803,(Varecia_rubra:1.270526,Varecia_variegata:1.270526):19.272276):52.460215);"	# read in phylogenetic data
primate_tree_raw <- read.tree(text = primate_tree_nex)		# convert to tree object
primate_tree <- rotate(primate_tree_raw, node = 16)			# rotate nodes to reorder tree
primate_tree <- rotate(primate_tree, node = 17)
primate_tree <- rotate(primate_tree, node = 19)
primate_tree <- rotate(primate_tree, node = 20)
primate_tree <- rotate(primate_tree, node = 22)
primate_tree <- rotate(primate_tree, node = 24)
primate_tree <- rotate(primate_tree, node = 25)

######################
## Integrate behavioral data with phylogenetic tree (preliminary to calculate relative brain size)
######################
prelim_primate_comp <- comparative.data(phy = primate_tree_raw, data = primate_data, names.col = species, 
	vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)		# combine tree with comparative data

######################
## Calculate relative brain size as residuals from regression of log(ecv) and log(bodywt)
######################
relative_brain <- data.frame(residuals(pgls(log(ecv) ~ log(bodywt), data = prelim_primate_comp)))	# calculate residuals from phylogenetic regression
relative_brain_species <- rownames(relative_brain)						# extract species names
rownames(relative_brain) <- NULL										# remove species names
relative_brain <- data.frame(relative_brain_species, relative_brain)	# create data frame of species names and residuals
names(relative_brain) <- c("species", "relative_brain")					# rename columns
primate_table$"Relative brain size" <- primate_data$relative_brain <- relative_brain$relative_brain			# append relative_brain to primate_data
primate_predictors <- primate_data[, c(2, 3, 7, 5, 14, 4, 6)]			# extract main predictors and indifference points
primate_logpredictors <- primate_data[, c(10, 11, 9, 13, 14, 12)]		# extract log predictors and indifference points

######################
## Analyze allometric variables
######################
# Calculate allometric correlations
allometric_vars <- primate_data[, c(9:11, 13)]		# extract body-size related columns
row.names(allometric_vars) <- primate_data$species	# assign species as row names
allometric_cor_matrix <- cor(allometric_vars)  		# create matrix of allometric correlation values
allometric_cor <- c(allometric_cor_matrix[2:4, 1], allometric_cor_matrix[3:4, 2], allometric_cor_matrix[4, 3])	# extract only single case of each correlation

# Conduct principal components analysis
allometric_pca <- principal(allometric_vars)		# calculate PCA on allometric variables
primate_data$allometric <- allometric_pca$scores	# extract PCA scores
allometric_propvar <- colSums(allometric_pca$loading*allometric_pca$loading)/dim(allometric_pca$loading)[1]	# calculate proportion of variance accounted for

predictors_pca <- principal(primate_logpredictors)

# Calculate correlations
allometric_itc_cor <- cor.test(log(primate_data$indiff), primate_data$allometric)			# calculate correlation of ITC and allometric
relative_brain_itc_cor <- cor.test(log(primate_data$indiff), primate_data$relative_brain)	# calculate correlation of ITC and relative brain size
group_size_itc_cor <- cor.test(log(primate_data$indiff), log(primate_data$groupsize))		# calculate correlation of ITC and group size

######################
## Integrate behavioral data with phylogenetic tree
######################
primate_comp <- comparative.data(phy = primate_tree_raw, data = primate_data, names.col = species, 
    vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)		# combine tree with comparative data

######################
## Analyze comparative data
######################
model_pgls <- pgls(log(indiff) ~ allometric + relative_brain + log(groupsize), 
  data = primate_comp, lambda ='ML')  	# calculate phylogenetic generalized least squares regression
summary(model_pgls)


# Create function that rounds p-values
round_pvalue <- function(pvalue) {
  if(pvalue < 0.01) {								# if p < 0.01
    new_pvalue <- "p < 0.01"						# print p < 0.01
  } else {				
    new_pvalue <- paste("p =", round(pvalue, 2))	# otherwise round p-value to 2 digits
  }
  new_pvalue										# return new p-value
}

# Calculate standardized coefficients (betas)
b <- summary(model_pgls)$coef[-1, 1]			# extract nonstandardized coefficients
sx1 <- sd(primate_data$allometric)				# calculate standard deviation of allometric data
sx2 <- sd(primate_data$relative_brain)			# calculate standard deviation of relative brain size data
sx3 <- sd(primate_data$loggroupsize)			# calculate standard deviation of log group size data
sy <- sd(primate_data$logindiff)				# calculate standard deviation of log indifference point data
allometric_beta <- round(b[1] * sx1/sy, 2)		# standardize allometric coefficient
relative_brain_beta <- round(b[2] * sx2/sy, 2)	# standardize relative brain size coefficient
groupsize_beta <- round(b[3] * sx3/sy, 2)		# standardize group size coefficient

# Extract t-values, F-statistics, p-values
model_fstat <- round(summary(model_pgls)$fstatistic[1],2)			# extract F statistic from overall model
model_pvalue <- round_pvalue(pf(summary(model_pgls)$fstatistic[1], summary(model_pgls)$fstatistic[2], summary(model_pgls)$fstatistic[3], lower = FALSE))		# calculate p-value for overall model
intercept_tvalue <- round(summary(model_pgls)$coefficients[1, 3], 2)		# extract t-values
allometric_tvalue <- round(summary(model_pgls)$coefficients[2, 3], 2)
relative_brain_tvalue <- round(summary(model_pgls)$coefficients[3, 3], 2)
groupsize_tvalue <- round(summary(model_pgls)$coefficients[4, 3], 2)
intercept_pvalue <- round_pvalue(summary(model_pgls)$coefficients[1, 4])	# extract p-values
allometric_pvalue <- round_pvalue(summary(model_pgls)$coefficients[2, 4])
relative_brain_pvalue <- round_pvalue(summary(model_pgls)$coefficients[3, 4])
groupsize_pvalue <- round_pvalue(summary(model_pgls)$coefficients[4, 4])

######################
## Calculate phylogenetic signal estimate (lambda)
######################
lambda <- round(summary(model_pgls)$param[2], 2)							# extract lambda
lambda_pvalues <- round(summary(model_pgls)$param.CI$lambda$bounds.p, 2)	# extract lambda estimate p-values for significant difference from 0 and 1

######################
## Analysis without orangutan data
######################
primate_data_no_orang <- subset(primate_data, species != "Pongo_pygmaeus")
primate_comp_no_orang <- comparative.data(phy = primate_tree_raw, data = primate_data_no_orang, names.col = species, vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE)  	# combine tree with comparative data
group_size_itc_cor_no_orang <- cor.test(log(primate_data_no_orang$indiff), log(primate_data_no_orang$groupsize))		# correlate ITC and group size
model_pgls_no_orang <- pgls(log(indiff) ~ allometric + relative_brain + log(groupsize), 
  data = primate_comp_no_orang, lambda ='ML')   								# calculate phylogenetic generalized least squares regression
summary(model_pgls_no_orang)
b_no_orang <- summary(model_pgls_no_orang)$coef[-1, 1]							# extract nonstandardized coefficients
sx_no_orang <- sd(primate_data_no_orang$loggroupsize)							# calculate standard deviation of log group size data
sy_no_orang <- sd(primate_data_no_orang$logindiff)								# calculate standard deviation of log indifference point data
groupsize_beta_no_orang <- round(b_no_orang[3] * sx_no_orang/sy_no_orang, 2)	# standardize group size coefficient
groupsize_pvalue_no_orang <- round_pvalue(summary(model_pgls_no_orang)$coefficients[4, 4])	# extract p-value for group size

######################
## Plot data
######################
## Plot hypothesis correlations
# Allometric
xy_allometric <- xyplot(indiff ~ allometric, data = primate_data, scales = list(y = list(log = TRUE, equispaced.log = FALSE, tck = c(1,0))), 
  xlab = "Allometric score", ylab = "Waiting time (s)", aspect = 1,  
  par.settings = list(par.ylab.text = list(cex = 2), par.xlab.text = list(cex = 2), 
    axis.text = list(cex = 1.5), layout.widths = list(ylab.axis.padding = 1, ylab = 1, 
    axis.left = 0.8, panel = 1, between = 1, axis.right = 0, ylab.right = 0, 
    right.padding = 0), axis.components = list(left = list(pad2 = 0))),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, type = c("p", "r"), cex = 2)
    panel.text(x = -1.7, y = log10(120), labels = "(a)", cex = 2)
  }
)
pdf(file = "figures/itc_allometric.pdf", width = 7, height = 5.5)
print(xy_allometric)
dev.off()

# Relative brain size
xy_relative_brain <- xyplot(indiff ~ relative_brain, data = primate_data,  scales = list(y = list(log = TRUE, equispaced.log = FALSE, labels = "", tck = c(0, 1))), 
  xlab = "Relative brain size", ylab = NULL, aspect = 1,  
  par.settings = list(par.ylab.text = list(cex = 2), par.xlab.text = list(cex = 2), 
    axis.text = list(cex = 1.5), layout.widths = list(ylab.axis.padding = 0, ylab = 0, axis.left = 0, panel = 1, 
    between = 0, axis.right = 0, ylab.right = 0, right.padding = 0), axis.components = list(left = list(pad2 = 0))),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, cex = 2)
    panel.text(x = -0.3, y = log10(120), labels = "(b)", cex = 2)
  }
)
pdf(file = "figures/itc_relative_brain.pdf", width = 7, height = 5.5)
print(xy_relative_brain)
dev.off()

# Group size
xy_groupsize <- xyplot(indiff ~ groupsize, data = primate_data, 
  xlab = "Group size", ylab = NULL, cex = 2, aspect = 1,
  scales = list(log = TRUE, equispaced.log = FALSE, y = list(labels = "")),
  par.settings = list(par.ylab.text = list(cex = 2), par.xlab.text = list(cex = 2), 
    axis.text = list(cex = 1.5), layout.widths = list(ylab.axis.padding = 0, ylab = 0, axis.left = 0, panel = 1, 
    between = 0, axis.right = 0, ylab.right = 0, right.padding = 0), axis.components = list(left = list(pad2 = 0))),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, cex = 2)
    panel.text(x = log10(2), y = log10(120), labels = "(c)", cex = 2)
  }                
)
pdf(file = "figures/itc_groupsize.pdf", width = 7, height = 5.5)
print(xy_groupsize)
dev.off()

# Combine separate plots into figure with three subfigures
pdf(file = "figures/itc_correlates.pdf", width = 12, height = 5)
print(xy_allometric, split = c(1,1,3,1), position = c(0.03, 0, 1, 1), panel.width = list(x = 3.55, unit = "in"), more=T)
print(xy_relative_brain, split = c(2,1,3,1), position = c(0.075, 0, 1, 1), panel.width = list(x = 3.5, unit = "in"), more=T)
print(xy_groupsize, split = c(3,1,3,1), position = c(0.05, 0, 1, 1), panel.width = list(x = 3.5, unit = "in"), more=F)
dev.off()

## Plot indifference points
indiff_data <- subset(raw_data, measure == "indiff")
indiff_data$latin_name <- factor(as.character(indiff_data$latin_name),  	# reorder species
  levels = c("Eulemur macaco", "Varecia rubra", "Varecia variegata", "Saguinus oedipus", 
    "Callithrix jacchus", "Sapajus apella", "Ateles geoffroyi", "Macaca fascicularis", 
    "Macaca mulatta", "Pongo pygmaeus", "Gorilla gorilla", "Pan paniscus", "Pan troglodytes"))

pdf("figures/indifference_points.pdf", width = 9)
par(mai = c(2, 0.9, 0, 0))
indifference_plot <- boxplot(value ~ latin_name, data = indiff_data, range = 0, col="#0000ff22",
  ylab = "Adjusted delay at indifference (s)", xlab = NULL, ylim = c(0, 160), las = 2, 
  cex = 1.5, cex.axis = 1.2, cex.lab = 1.5, lwd = 0.5, boxwex = 0.5, outlier = F, color = T
)          							# generate boxplot of species indifference points for 2v6 condition
mean_indiff <- tapply(indiff_data$value, indiff_data$latin_name, mean)			# calculate means
points(seq(indifference_plot$n), mean_indiff, pch = 17, cex = 1.5, lwd = 2)		# overlay means on boxplot
beeswarm(value ~ latin_name, data = indiff_data, pch = 1, cex = 1.15, col = "grey30", add = T)	# overlay beeswarm plot of individual data points
dev.off()

## Plot phylogenetic tree
pdf(file = "figures/phylogeny.pdf", width = 7, height = 5.5)# print tree to file
plot(primate_tree)
dev.off()

## Plot all correlations
# Create function for upper panels of correlation plot
panel.cor <- function(x, y, ...) {
  # Calculate center of x dimension and assign x to log or linear scale
  if(min(x) > 0) {    		# if the minimum value is positive, this variable is plotted on a log scale (hack)
    rx <- sqrt(min(x) * max(x)) # calculate midpoint of log scale
    scaledx <- log(x) 			# log transform x
  } else { 						# if the minimum value is negative, this variable is plotted on a linear scale (hack)
    rx <- mean(x) 				# calculate midpoint of linear scale
    scaledx <- x   				# keep x in linear scale
  }
  # Calculate center of y dimension and assign y to log or linear scale
  if(min(y) > 0) {  			# if the minimum value is positive, this variable is plotted on a log scale (hack)
    ry <- sqrt(min(y) * max(y)) # calculate midpoint of log scale
    scaledy <- log(y) 			# log transform y
  } else {  					# if the minimum value is negative, this variable is plotted on a linear scale (hack)
    ry <- mean(y) 				# calculate midpoint of linear scale
    scaledy <- y   				# keep y in linear scale
  }
  # Calculate correlation coefficients for each panel
  r <- round(abs(cor(scaledx, scaledy)), 2) # calculate correlation coefficient
  if(r == 0) r = "0.00"						# if correlation coefficient = 0, print "0.00"
  # Color background of allometric variables
  if(r > 0.9) {								# if correlation coefficient > 0.9 (therefore, is one of the allometric variables)
    minx <- min(x) - abs(min(x)) * 0.50		# calculate minimum of x-axis
    maxx <- max(x) + abs(max(x)) * 0.50		# calculate maximum of x-axis
    miny <- min(y) - abs(min(y)) * 0.50		# calculate minimum of y-axis	
    maxy <- max(y) + abs(max(y)) * 0.50		# calculate maximum of y-axis
    rect(minx, miny, maxx, maxy, col = "gray90", border = NA)	# print gray rectangle as background for allometric variables
  }
  text(rx, ry, r, cex = 2)					# print correlation coefficient
}

# Create function for lower panels of correlation plot
panel.mypoints = function(x, y, ...) {
  # Determine whether scale should be log or linear
  if(min(x) > 0) {  			# if the minimum value is positive, this variable is plotted on a log scale (hack)
    rx <- sqrt(min(x) * max(x)) # calculate midpoint of log scale
    scaledx <- log(x) 			# log transform x
  } else {  					# if the minimum value is negative, this variable is plotted on a linear scale (hack)
    rx <- mean(x) 				# calculate midpoint of linear scale
    scaledx <- x   				# keep x in linear scale
  }
  # Calculate center of y dimension and assign y to log or linear scale
  if(min(y) > 0) {  			# if the minimum value is positive, this variable is plotted on a log scale (hack)
    ry <- sqrt(min(y) * max(y)) # calculate midpoint of log scale
    scaledy <- log(y) 			# log transform y
  } else {  					# if the minimum value is negative, this variable is plotted on a linear scale (hack)
    ry <- mean(y) 				# calculate midpoint of linear scale
    scaledy <- y   				# keep y in linear scale
  }
  # Calculate correlation coefficients for each panel
  r <- round(abs(cor(scaledx, scaledy)), 2) # calculate correlation coefficient
  if(r == 0) r = "0.00"						# if correlation coefficient = 0, print "0.00"
  # Color background of allometric variables
  if(r > 0.9) {								# if correlation coefficient > 0.9 (therefore, is one of the allometric variables)
    minx <- min(x) - abs(min(x)) * 0.50		# calculate minimum of x-axis
    maxx <- max(x) + abs(max(x)) * 0.50		# calculate maximum of x-axis
    miny <- min(y) - abs(min(y)) * 0.50		# calculate minimum of y-axis	
    maxy <- max(y) + abs(max(y)) * 0.50		# calculate maximum of y-axis
    rect(minx, miny, maxx, maxy, col = "gray90", border = NA)	# print gray rectangle as background for allometric variables
  }
  points(x, y)	# plot points
  if(r > 0.9) {  							# if correlation coefficient > 0.9 (therefore, is one of the allometric variables)
  abline(coef = coef(lm(log10(y) ~ log10(x))))
  }
}

# Create function for diagonal panels of correlation plot
diag_labels <- c("Body mass", "Absolute\nbrain volume", "Lifespan", "Home range")	# list labels for allometric variables
diag_colors <-  c("gray80", "gray80", "gray80", "gray80", "white", "white", "white")# assign background colors
panel.diag <- function(x, y, labels, ...) {
  minx <- min(x)		# calculate minimum of x-axis
  maxx <- max(x)		# calculate maximum of x-axis
  miny <- min(y)		# calculate minimum of y-axis
  maxy <- max(y)		# calculate maximum of y-axis
  if(labels %in% diag_labels) {	# if label is allometric label
    rect(0.1, 0.1, 10, 10, col = "gray90", border = NA)		# print gray rectangle as background for allometric variables
  } 
  text(minx, miny, labels, cex = 1.85)	# print label
}

# Print correlations plot
options(scipen=3)    						# set options to remove scientific notation from axes
pdf(file = "figures/all_correlates.pdf", width = 11, height = 8)	# save as PDF
pairs(primate_predictors, log = c(1:4, 6:7), upper.panel = panel.cor, lower.panel = panel.mypoints,  text.panel = panel.diag,
      labels = c("Body mass", "Absolute\nbrain volume", "Lifespan", "Home range", "Relative\nbrain size", "Group size", "Indifference\npoints"), 
      cex.labels = 1.9, cex.axis = 1.35, cex = 1.25, gap = 1.5
)												# create correlation graphs
dev.off()
