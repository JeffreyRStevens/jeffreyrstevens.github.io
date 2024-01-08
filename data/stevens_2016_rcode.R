###################################################
### stevens_dataS4.R
### Created by Jeffrey R. Stevens on 13 May 2011 (jeffrey.r.stevens@gmail.com),
###	finalized on 21 Jan 2015
### Summary: This script calculates descriptive statistics, runs model comparisons,
###     and generates figures for the analysis of human intertemporal choice data.
### Instructions: Place this file and the data files (stevens_expt[1/2]_data.csv)
### 	in the same directory.  Create a folder called "figures". Set the R
### 	working directory to this directory.  At the R command prompt, type 
### 	> source("stevens_dataS4.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PDF versions of the figures in the figures directory.
### Uses: This script can be reproduced and modified for personal and scientific use.
### Data files: 
###  stevens_expt[1, 2, or 3]_data--choice data for experiment 1/2/3
### Description of the data columns:
###  subject - participant number
###  question - question number
###  type - type of question (binary, staircase, similarity)
###  small_amount - small amount for question (in euros)
###  large_amount - large amount for question (in euros)
###  short_delay - short delay for question (in days)
###  long_delay - long delay for question (in days)
###  choice - choice for smaller, sooner (1) or larger, later (2)
###  rt - decision time from start of trial until choice is made
###################################################

##############################
### Load libraries, R version, and define functions
##############################
rm(list=ls())					# clear all variables	
library(bbmle)  			# needed for mle2
library(car)					# needed for leveneTest
library(epicalc)			# needed for aggregate with multiple functions
library(Hmisc)				# needed for xYplot
library(foreach)			# needed for foreach
library(lattice)  		# needed for lattice plots
library(latticeExtra) # needed for layer
library(plyr)					# needed for ddply
library(xtable) 			# needed for xtable
library(zoo)  				# needed for rollmean
ver <- getRversion()	# get R version

###############
## Define functions
###############

#####
## Calculate within-subjects confidence intervals
##  from Cookbook for R by Winston Chang (http://wiki.stdout.org/rcookbook/)
#####

# Summarize data.
# Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#   data: a data frame.
#   measurevar: the name of a column that contains the variable to be summariezed
#   groupvars: a vector containing names of columns that contain grouping variables
#   na.rm: a boolean that indicates whether to ignore NA's
#   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval=.95, .drop = TRUE) {
  length2 <- function (x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,	.fun= function(xx, col, na.rm) {	# This is does the summary; it's not easy to understand...
		c(N = length2(xx[,col], na.rm=na.rm), mean = mean(xx[,col], na.rm=na.rm), sd = sd(xx[,col], na.rm=na.rm))
	}, measurevar,	na.rm)
  datac <- rename(datac, c("mean"=measurevar))	# Rename the "mean" column 
  datac$se <- datac$sd / sqrt(datac$N)  				# Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1) # Calculate t-statistic for confidence interval
  datac$ci <- datac$se * ciMult
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data = NULL, idvar, measurevar, betweenvars = NULL, na.rm = FALSE, .drop = TRUE) {
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop, .fun = function(xx, col, na.rm) {
		c(subjMean = mean(xx[,col], na.rm=na.rm))
	}, measurevar, na.rm)
  data <- merge(data, data.subjMean)  										# Put the subject means with original data
  measureNormedVar <- paste(measurevar, "Normed", sep="")	# Get the normalized data in a new column
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  data$subjMean <- NULL 																	# Remove this subject mean column
  return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data = NULL, measurevar, betweenvars = NULL, withinvars = NULL, idvar = NULL, na.rm = FALSE, conf.interval = 0.95, .drop = TRUE) {
  factorvars <- sapply(data[, c(betweenvars, withinvars), drop = FALSE], FUN = is.factor)  # Ensure that the betweenvars and withinvars are factors
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ", paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  data <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)					# Norm each subject's data 
  measureNormedVar <- paste(measurevar, "Normed", sep="")  																	# This is the name of the new column
  data[,measurevar] <- data[,measureNormedVar]  																						# Replace the original data column with the normed one
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars), na.rm=na.rm, 	# Collapse the normed data - now we can treat between and within vars the same,        
                     conf.interval=conf.interval, .drop=.drop)
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(sapply(datac[,withinvars, drop = FALSE], FUN=nlevels))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  datac$sd <- datac$sd * correctionFactor  # Apply the correction factor
  datac$se <- datac$se * correctionFactor
  datac$ci <- datac$ci * correctionFactor
  return(datac)
}

#####
## Creates a function that is the inverse of %in%, that is it finds items that do not match the vector
#####
"%notin%" <- function(x, table) {  # create function analogous to %in% that searches for items not in a vector
  match(x, table, nomatch = 0) == 0
}

##############################
## Experiment 1--Magnitude effect
##############################

################
## Load and prepare data
################

# Load all data
all_data1 <- read.csv("stevens_dataS1.csv")

subject_nums1 <- unique(all_data1$subject)	# find unique subject numbers
num_subjects1 <- length(subject_nums1)  		# find number of subjects

## Staircase (training data)
staircase1 <- subset(all_data1, type == "staircase" & question < 18 & small_amount != 10)	# subset staircase data with adjusting amount (question < 18) and with different small and large amounts (small_amount != 10)

# Plot staircase data
staircase_means1 <- aggregate(choice ~ long_delay * short_delay * subject, staircase1, mean)		# aggregate staircase choice data by long delay, short delay, and subject
staircase_means1$short_delay <- as.factor(staircase_means1$short_delay)													# convert short_delay to factor
col.blind <- c("#0072B2", "#D55E00", "#009E73")																									# assign color-blind safe colors for plots
staircase_plot1 <- xyplot(choice * 100 ~ long_delay | as.factor(subject), group = short_delay, data = staircase_means1, 
	xlab = "Long delay (days)", ylab = "Mean percent choice for LL", type = "b", col = col.blind, as.table = TRUE,
  key = list(space = "top", padding.text = 3, cex = 1, title = "Short delay", cex.title = 1.25, # create key
    text = list(c("0 days", "6 days", "12 days"), adj = 1), 
    lines = list(col = col.blind), pch = 1, lwd = 1.5, type = "b", divide = 1),
  par.settings = list(par.xlab.text = list(cex = 1.65), par.ylab.text = list(cex = 1.65))
)
cairo_pdf(file = "figures/staircase_plot1.pdf", width = 8, height = 11)													# create PDF
plot(staircase_plot1)																																			# print plot
dev.off()																																												# turn off PDF device

## Similarity data
sim_data1 <- subset(all_data1, type == "similarity" & question != 8 & question != 11 & question != 20 & question != 34 & question != 48) # subset all similarity questions that are not duplicates (questions 8, 11, 34, 48,) or errors (question 20)
sim_data1$attribute <- ifelse(sim_data1$question < 24, "amount", "delay")	# create attribute factor for whether amount or delay similarity judgment
names(sim_data1)[8] <- "similarity"																				# rename similarity column
sim_amt1 <- subset(sim_data1, attribute == "amount")											# subset similarity judgments for amount
sim_amt1 <- sim_amt1[, c(1, 4, 5, 8)]																			# extract relevant columns
sim_delay1 <- subset(sim_data1, attribute == "delay")											# subset similarity judgments for delay
sim_delay1 <- sim_delay1[, c(1, 6, 7, 8)]																	# extract relevant columns

## Binary (validation) data
binary1 <- subset(all_data1, type == "binary")																															# subset binary choice data
validation1 <- merge(binary1, sim_amt1, by = c("subject", "small_amount", "large_amount"), all.x = TRUE)		# merge binary choice data with similarity amount judgment data
validation1 <- validation1[!duplicated(validation1), ]																											# remove duplicated entries from merge
names(validation1)[10] <- "sim_amt"																																					# rename sim_amt column
validation1 <- merge(validation1, sim_delay1, by = c("subject", "short_delay", "long_delay"), all.x = TRUE)	# merge binary choice data with similarity delay judgment data
validation1 <- validation1[!duplicated(validation1), ]																											# remove duplicated entries from merge
names(validation1)[11] <- "sim_delay"																																				# rename sim_amt column
validation1$sim_diff <- ifelse(validation1$sim_amt == 1 & validation1$sim_delay == 1, 1, ifelse(validation1$sim_amt == 1 & validation1$sim_delay == 0, 2, ifelse(validation1$sim_amt == 0 & validation1$sim_delay == 1, 3, ifelse(validation1$sim_amt == 0 & validation1$sim_delay == 0, 4, NA))))  # categorize difference in similarity judgments (1 = both similar, 2 = amounts similar, delays dissimilar, 3 = amounts dissimilar, delays similar, 4 = both dissimilar)

###############
## Participant info
###############
## Age data experiment 1
age1 <- c(26, 23, 28, 26, 23, 24, 25, 25, 27, 31, 28, 22, 27, 26, 22, 23, 25, 29, 31, 21, 21, 28, 24, 30, 27, 25, 30, 23, 26, 26, 24, 27, 22, 21, 21, 28, 33, 19, 25, 25, 26, 25, 27, 21, 24, 24, 28, 26, 26, 32, 26, 29, 29, 25, 27, 25, 28, 25, 30, 24, 28, 22, 28, 26) 
age1_m <- mean(age1)	# calculate mean age
age1_sd <- sd(age1)		# calculate standard deviation of age
age1_min <- min(age1)	# calculate minimum age
age1_max <- max(age1)	# calculate maximum age

## Payment data experiment 1
pay1 <- c(8, 8, 7, 5, 9, 9, 12, 7, 7, 8, 8, 12, 10, 7, 5.5, 8, 8.5, 3, 9, NA, 6, 10, 8, 6, 8, 9, 9, 8, 8, 7, 6.7, 7, 6, 4, 6, 5, 8, 8, 8, 5, 10, 8, 4.9, 6, 8, 9, 6, 8, 1, 9, 10, 6, 7, 5, 8, 1, 4, 6, 10, 3, 10, 4.1, 7, 15) 
pay1_m <- mean(pay1, na.rm = TRUE)	# calculate mean age
pay1_sd <- sd(pay1, na.rm = TRUE)		# calculate standard deviation of age
pay1_min <- min(pay1, na.rm = TRUE)	# calculate minimum age
pay1_max <- max(pay1, na.rm = TRUE)	# calculate maximum age

## Sex data experiment 1
males1 <- c(9, 10, 11, 12, 13, 25, 27, 29, 31, 34, 38, 44, 47, 49, 53, 55, 56, 58, 60, 62, 63, 64, 65, 66, 69, 73, 75, 77, 80)
females1 <- c(22, 15, 16, 26, 28, 30, 32, 33, 35, 36, 37, 39, 40, 41, 42, 43, 45, 46, 48, 50, 51, 52, 54, 57, 59, 61, 67, 68, 70, 71, 72, 74, 76, 78, 79)

###############
## Model selection
###############

## Models
# Inverse logit for error function
invlogit <- function(z)
{
  return(1 / (1 + exp(-z)))
}
# Exponential discounting
expon <- function(a, t, delta) {
  return(a * exp(-delta * t))
}
# Hyperbolic (Mazur) discounting
hyper_m <- function(a, t, k) {
  return(a / (1 + k * t))
}
# Hyperbolic (Rachlin) discounting
hyper_r <- function(a, t, k, s) {
  return(a / (1 + k * t ^ s))
}
# Hyperbolic (Kirby) discounting
hyper_k <- function(a, t, k, m) {
  return(a / (1 + k * abs(a) ^ m * t))
}
# Hyperbolic (Loewenstein & Prelec) discounting
hyper_lp <- function(a, t, alpha, beta) {
  return(a / ((1 + alpha * t) ^ (beta / alpha)))
}
# Additive discounting (Killeen 2009)
additive <- function(a, t, lambda) {
  return(a - lambda * t)
}

####
## Fit staircase (training) data with maximum likelihood estimation (MLE)
####

## Conduct subject-wise maximum likelihood estimation
# Prepare vectors and data frames
good_subjects1 <- NA # initiate vector of 'good' subjects (subject for whom the MLE converges)

fitsNAs <- rep(NA, num_subjects1)	# initiate vector of NAs to fill empty data frame
mle_fits1 <- data.frame(subject = subject_nums1, expon_AICc = fitsNAs, delta = fitsNAs, hyper_m_AICc = fitsNAs, k_m = fitsNAs, hyper_r_AICc = fitsNAs, k_r = fitsNAs, sigma = fitsNAs, hyper_k_AICc = fitsNAs, k_k = fitsNAs, mu = fitsNAs, hyper_lp_AICc = fitsNAs, alpha = fitsNAs, beta = fitsNAs, add_AICc = fitsNAs, lambda = fitsNAs, expon_diffs = fitsNAs, hyper_m_diffs = fitsNAs, hyper_r_diffs = fitsNAs, hyper_k_diffs = fitsNAs, hyper_lp_diffs = fitsNAs, add_diffs = fitsNAs)	# initiate data frame for MLE analysis

# Conduct MLE for each subject
foreach(subj = subject_nums1) %do% {								# for each subject
  curr_subj <- subset(staircase1, subject == subj)	# assign curr_subj to the current subject's data

	# Calculate MLE for each subject
  if(sum(curr_subj$choice) < 147 & subj != 36) {		# remove subjects whose MLE do not converge
    good_subjects1 <- c(good_subjects1, subj)  			# append current subject to good_subjects1
    expon_fit <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = curr_subj, start = list(delta = 0.05))			# fit exponential discounting model
    hyper_m_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = curr_subj, start = list(k = 0.05))					# fit hyperbolic (Mazur) discounting model
    hyper_r_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, s) - hyper_r(small_amount, short_delay, k, s)), size = 1), data = curr_subj, start = list(k = 0.05, s = 0.5))  # fit hyperbolic (Rachlin) discounting model
    hyper_k_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, m) - hyper_k(small_amount, short_delay, k, m)), size = 1), data = curr_subj, start = list(k = 0.05, m = -0.25))  # fit hyperbolic (Kirby) discounting model
    if(subj == 25) {							# if subject = 25, use different starting parameters
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.5))				# fit hyperbolic (L&P) discounting model
    } else if(subj == 77) {				# if subject = 77, use different starting parameters
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.5, beta = 0.5))    		# fit hyperbolic (L&P) discounting model
    } else {
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.05))			# fit hyperbolic (L&P) discounting model
    }
    add_fit <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = curr_subj, start = list(lambda = 0.05))	# fit additive discounting model

    # Calculate AICc for each model and extract parameters    
    mle_fits1$expon_AICc[which(mle_fits1$subject == subj)[1]] <- expon_aic <- AICc(expon_fit, nobs = 1, k = 2)					# calculate AICc for exponential discounting
    mle_fits1$delta[which(mle_fits1$subject == subj)[1]] <- coef(expon_fit)																							# extract delta
    
    mle_fits1$hyper_m_AICc[which(mle_fits1$subject == subj)[1]] <- hyper_m_aic <- AICc(hyper_m_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur) discounting
    mle_fits1$k_m[which(mle_fits1$subject == subj)[1]] <- coef(hyper_m_fit)																							# extract k
    
    mle_fits1$hyper_r_AICc[which(mle_fits1$subject == subj)[1]] <- hyper_r_aic <- AICc(hyper_r_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Rachlin) discounting
    mle_fits1$k_r[which(mle_fits1$subject == subj)[1]] <- coef(hyper_r_fit)[1]																					# extract k
    mle_fits1$sigma[which(mle_fits1$subject == subj)[1]] <- coef(hyper_r_fit)[2]																				# extract sigma
    
    mle_fits1$hyper_k_AICc[which(mle_fits1$subject == subj)[1]] <- hyper_k_aic <- AICc(hyper_k_fit, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby) discounting
    mle_fits1$k_k[which(mle_fits1$subject == subj)[1]] <- coef(hyper_k_fit)[1]																					# extract k
    mle_fits1$mu[which(mle_fits1$subject == subj)[1]] <- coef(hyper_k_fit)[2]																						# extract mu

		mle_fits1$hyper_lp_AICc[which(mle_fits1$subject == subj)[1]] <- hyper_lp_aic <- AICc(hyper_lp_fit, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P) discounting
    mle_fits1$alpha[which(mle_fits1$subject == subj)[1]] <- coef(hyper_lp_fit)[1]																				# extract alpha
    mle_fits1$beta[which(mle_fits1$subject == subj)[1]] <- coef(hyper_lp_fit)[2]																				# extract beta
    
    mle_fits1$add_AICc[which(mle_fits1$subject == subj)[1]] <- add_aic <- AICc(add_fit, nobs = 1, k = 2)								# calculate AICc for additive discounting
    mle_fits1$lambda[which(mle_fits1$subject == subj)[1]] <- coef(add_fit)																							# extract lambda

    subj_aics <- c(expon_aic, hyper_m_aic, hyper_r_aic, hyper_k_aic, hyper_lp_aic, add_aic)															# concatenate AICc values for all models
    aic_diffs <- subj_aics - min(subj_aics)																																							# calculate AICc differences
    aic_weights <- (exp(-0.5 * aic_diffs)) / (sum(exp(-0.5 * aic_diffs)))																								# calculate AICc weights
    mle_fits1[which(mle_fits1$subject == subj)[1], 17:22] <- aic_diffs																									# assign AICc differences to mle_first1
  }
}

good_subjects1 <- good_subjects1[-1]					# remove initial NA
num_good_subjects1 <- length(good_subjects1)	# find number of good subjects

aic_individual_median1 <- apply(mle_fits1[, c(2, 4, 6, 9, 12, 15)], 2, median, na.rm = TRUE)  # calculate median AICc values across subjects
aic_diff_individual_median1 <- apply(mle_fits1[, c(17:22)], 2, median, na.rm = TRUE)  				# calculate median AICc differences across subjects

parameters1 <- mle_fits1[, c(3, 5, 7, 8, 10, 11, 13, 14, 16)]	# extract parameters for each subject and model
parameters1 <- subset(parameters1, !is.na(delta))							# remove subjects with no parameter estimates

## Conduct MLE for aggregate data
train1 <- subset(staircase1, subject %in% good_subjects1)	# subset staircase data from 'good' subjects

expon_fit_group1 <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = train1, start = list(delta = 0.05))	# calculate MLE for exponential discounting model

hyper_m_fit_group1 <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = train1, start = list(k = 0.05))	# calculate MLE for hyperbolic (Mazur) discounting model

hyper_r_fit_group1 <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, sigma) - hyper_r(small_amount, short_delay, k, sigma)), size = 1), data = train1, start = list(k = 0.05, sigma = 0.5))  # calculate MLE for hyperbolic (Rachlin) discounting model

hyper_k_fit_group1 <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, mu) - hyper_k(small_amount, short_delay, k, mu)), size = 1), data = train1, start = list(k = 0.05, mu = -0.5))  # calculate MLE for hyperbolic (Kirby) discounting model

hyper_lp_fit_group1 <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = train1, start = list(alpha = 0.05, beta = 0.05))	# calculate MLE for hyperbolic (L&P) discounting model

add_fit_group1 <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = train1, start = list(lambda = 0.05))	# calculate MLE for exponential discounting model

expon_aic_group1 <- AICc(expon_fit_group1, nobs = 1, k = 2)				# calculate AICc for exponential discounting
hyper_m_aic_group1 <- AICc(hyper_m_fit_group1, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur)  discounting
hyper_r_aic_group1 <- AICc(hyper_r_fit_group1, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Rachlin)  discounting
hyper_k_aic_group1 <- AICc(hyper_k_fit_group1, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby)  discounting
hyper_lp_aic_group1 <- AICc(hyper_lp_fit_group1, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P)  discounting
add_aic_group1 <- AICc(add_fit_group1, nobs = 1, k = 2)						# calculate AICc for additive discounting

aics_group1 <- c(expon_aic_group1, hyper_m_aic_group1, hyper_r_aic_group1, hyper_k_aic_group1, hyper_lp_aic_group1, add_aic_group1)	# concatenate AICc values for all models
aic_diffs_group1 <- aics_group1 - min(aics_group1)																																									# calculate AICc differences
aic_weights_group1 <- (exp(-0.5 * aic_diffs_group1)) / (sum(exp(-0.5 * aic_diffs_group1)))																					# calculate AICc weights
evidence_ratios_group1 <- c(aic_weights_group1[3] / aic_weights_group1[1], aic_weights_group1[3] / aic_weights_group1[2], aic_weights_group1[3] / aic_weights_group1[4], aic_weights_group1[3] / aic_weights_group1[5], aic_weights_group1[3] / aic_weights_group1[6])																																# calculate evidence ratios

####
## Predict binary choice (validation) data using parameters fit from staircase (training) data
####

## Prepare data
validation1$sim_domain <- ifelse(is.na(validation1$sim_delay), NA, ifelse(is.na(validation1$sim_amt), NA, ifelse(validation1$sim_diff == 1 | validation1$sim_diff == 4, 0, 1))) 	# create column with domain, non-domain, and unrated questions
valid1 <- subset(validation1, subject %in% good_subjects1 & !is.na(sim_diff))		# create subset of questions from good subjects (with variance in choices) that have similarity judgments
valid1$subject <- factor(valid1$subject)

# Establish number of models and parameters
num_models1 <- 6  												# assign total number of models fitted (not including similarity models)
total_num_models1 <- 2 * num_models1 + 1	# include similarity models in total number of models
num_params1 <- dim(parameters1)[2]				# assign total number of parameters estimated

## Create xtable of questions and mean responses
valid1_aggregation <- aggregate(valid1$choice, by = list(valid1$small_amount, valid1$large_amount, valid1$short_delay, valid1$long_delay, valid1$question), FUN = c("mean", "sd")) 												# calculate mean and sd for each question
valid1_aggregation <- valid1_aggregation[, -5]    							# remove question column
valid1_aggregation <- valid1_aggregation[, c(3, 4, 1, 2, 5, 6)] # reorder columns with delays first
names(valid1_aggregation) <- c("Short delay", "Long delay", "Small amount", "Large amount", "Mean choice for LL", "Standard deviation")
valid1_aggregation <- valid1_aggregation[order(valid1_aggregation$"Short delay", valid1_aggregation$"Long delay", valid1_aggregation$"Small amount", valid1_aggregation$"Large amount", valid1_aggregation$"Mean choice for LL"), ]																				# reorder rows
valid1_xtable <- xtable(valid1_aggregation, caption = 'Questions and Mean Responses for Experiment 1', label = "tab:EXP1QUEST")	# create xtable for import into LaTeX

## Predict binary choice data (generalization/validation)
# Prepare data frames
preds_trials1 <- data.frame(valid1$subject, matrix(rep(NA, length(valid1$subject) * (total_num_models1 - 1)), nrow = length(valid1$subject)))	# create data frame of NAs for trial predictions
names(preds_trials1) <- c("subject", "expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith")
preds_subj1 <- data.frame(good_subjects1, matrix(rep(NA, num_good_subjects1 * (total_num_models1 - 1)), nrow = num_good_subjects1))						# create data frame of NAs for subject predictions
names(preds_subj1) <- names(preds_trials1) 
num_questions1 <- length(unique(valid1$question))	# find the number of unique questions

# Exponential discounting model
valid1$delta <- rep(parameters1$delta, each = num_questions1)   	# create column of consecutives instances of each subject's fitted parameter for each question
exp_pred1 <- ifelse(valid1$small_amount * exp(-valid1$delta * valid1$short_delay) < valid1$large_amount * exp(-valid1$delta * valid1$long_delay), 1, ifelse(valid1$small_amount * exp(-valid1$delta * valid1$short_delay) > valid1$large_amount * exp(-valid1$delta * valid1$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials1$expo <- ifelse(exp_pred1 == valid1$choice, 1, 0)
preds_subj1$expo <- aggregate(preds_trials1$expo, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Hyperbolic discounting model Mazur
valid1$k_m <- rep(parameters1$k_m, each = num_questions1)
hyp_m_pred1 <- ifelse(valid1$small_amount / (1 + valid1$k_m * valid1$short_delay) < valid1$large_amount / (1 + valid1$k_m * valid1$long_delay), 1, ifelse(valid1$small_amount / (1 + valid1$k_m * valid1$short_delay) > valid1$large_amount / (1 + valid1$k_m * valid1$long_delay), 0, NA))
preds_trials1$hyper_m <- ifelse(hyp_m_pred1 == valid1$choice, 1, 0)
preds_subj1$hyper_m <- aggregate(preds_trials1$hyper_m, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Hyperbolic discounting model Rachlin
valid1$k_r <- rep(parameters1$k_r, each = num_questions1)
valid1$sigma <- rep(parameters1$sigma, each = num_questions1)
hyp_r_pred1 <- ifelse(valid1$small_amount / (1 + valid1$k_r * valid1$short_delay ^ valid1$sigma)  < valid1$large_amount / (1 + valid1$k_r * valid1$long_delay ^ valid1$sigma), 1, ifelse(valid1$small_amount / (1 + valid1$k_r * valid1$short_delay ^ valid1$sigma) > valid1$large_amount / (1 + valid1$k_r * valid1$long_delay ^ valid1$sigma), 0, NA))
preds_trials1$hyper_r <- ifelse(hyp_r_pred1 == valid1$choice, 1, 0)
preds_subj1$hyper_r <- aggregate(preds_trials1$hyper_r, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Hyperbolic discounting model Kirby
valid1$k_k <- rep(parameters1$k_k, each = num_questions1)
valid1$mu <- rep(parameters1$mu, each = num_questions1)
hyp_k_pred1 <- ifelse(valid1$small_amount / (1 + valid1$k_k * valid1$small_amount ^ valid1$mu * valid1$short_delay) < valid1$large_amount / (1 + valid1$k_k * valid1$large_amount ^ valid1$mu * valid1$long_delay), 1, ifelse(valid1$small_amount / (1 + valid1$k_k * valid1$small_amount ^ valid1$mu * valid1$short_delay) > valid1$large_amount / (1 + valid1$k_k * valid1$large_amount ^ valid1$mu * valid1$long_delay), 0, NA))
preds_trials1$hyper_k <- ifelse(hyp_k_pred1 == valid1$choice, 1, 0)
preds_subj1$hyper_k <- aggregate(preds_trials1$hyper_k, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Hyperbolic discounting model L&P
valid1$alpha <- rep(parameters1$alpha, each = num_questions1)
valid1$beta <- rep(parameters1$beta, each = num_questions1)
hyp_lp_pred1 <- ifelse(valid1$small_amount / (1 + valid1$alpha * valid1$short_delay) ^ (valid1$beta / valid1$alpha)  < valid1$large_amount / (1 + valid1$alpha * valid1$long_delay) ^ (valid1$beta / valid1$alpha), 1, ifelse(valid1$small_amount / (1 + valid1$alpha * valid1$short_delay) ^ (valid1$beta / valid1$alpha) > valid1$large_amount / (1 + valid1$alpha * valid1$long_delay) ^ (valid1$beta / valid1$alpha), 0, NA))
preds_trials1$hyper_lp <- ifelse(hyp_lp_pred1 == valid1$choice, 1, 0)
preds_subj1$hyper_lp <- aggregate(preds_trials1$hyper_lp, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Arithmetic model
valid1$lambda <- rep(parameters1$lambda, each = num_questions1)
arith_pred1 <- ifelse(valid1$small_amount - valid1$lambda * valid1$short_delay < valid1$large_amount - valid1$lambda * valid1$long_delay, 1,  ifelse(valid1$small_amount - valid1$lambda - valid1$short_delay > valid1$large_amount - valid1$lambda - valid1$long_delay, 0, NA))
preds_trials1$arith <- ifelse(arith_pred1 == valid1$choice, 1, 0)
preds_subj1$arith <- aggregate(preds_trials1$arith, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (exponential)--first look for similarity differences, then use exponential discounting if similarity doesn't distinguish
sim_e_pred1 <- ifelse(valid1$sim_diff == 2, 0, ifelse(valid1$sim_diff == 3, 1, exp_pred1))
preds_trials1$sim_e <- ifelse(sim_e_pred1 == valid1$choice, 1, 0)
preds_subj1$sim_e <- aggregate(preds_trials1$sim_e, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--Mazur)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hm_pred1 <- ifelse(valid1$sim_diff == 2, 0,	ifelse(valid1$sim_diff == 3, 1, hyp_m_pred1))
preds_trials1$sim_h_m <- ifelse(sim_hm_pred1 == valid1$choice, 1, 0)
preds_subj1$sim_h_m <- aggregate(preds_trials1$sim_h_m, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--Rachlin)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hr_pred1 <- ifelse(valid1$sim_diff == 2, 0,  ifelse(valid1$sim_diff == 3, 1, hyp_r_pred1))
preds_trials1$sim_h_r <- ifelse(sim_hr_pred1 == valid1$choice, 1, 0)
preds_subj1$sim_h_r <- aggregate(preds_trials1$sim_h_r, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--Kirby)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hk_pred1 <- ifelse(valid1$sim_diff == 2, 0,  ifelse(valid1$sim_diff == 3, 1, hyp_k_pred1))
preds_trials1$sim_h_k <- ifelse(sim_hk_pred1 == valid1$choice, 1, 0)
preds_subj1$sim_h_k <- aggregate(preds_trials1$sim_h_k, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--L&P)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hlp_pred1 <- ifelse(valid1$sim_diff == 2, 0,  ifelse(valid1$sim_diff == 3, 1, hyp_lp_pred1))
preds_trials1$sim_h_lp <- ifelse(sim_hlp_pred1 == valid1$choice, 1, 0)
preds_subj1$sim_h_lp <- aggregate(preds_trials1$sim_h_lp, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (Arithmetic)--first look for similarity differences, then use Arithmetic discounting if similarity doesn't distinguish
sim_a_pred1 <- ifelse(valid1$sim_diff == 2, 0, ifelse(valid1$sim_diff == 3, 1, arith_pred1))
preds_trials1$sim_arith <- ifelse(sim_a_pred1 == valid1$choice, 1, 0)
preds_subj1$sim_arith <- aggregate(preds_trials1$sim_arith, by = list(preds_trials1$subject), FUN = mean)[, 2]

# Similarity model (Leland)
sim_l_pred1 <- ifelse(valid1$sim_diff == 2, 0, ifelse(valid1$sim_diff == 3, 1, NA))	# create array differentiating between higher sim for amount, higher sim for delay, and equal sim

# Add Leland similarity model and other columns to preds_trials
preds_trials1$sim_l <- ifelse(sim_l_pred1 == valid1$choice, 1, 0)	# create column of correct predictions for similarity
preds_trials1$sim_diff <- valid1$sim_diff													# append similarity difference values
preds_trials1$sim_domain <- ifelse(is.na(valid1$sim_delay), NA, ifelse(is.na(valid1$sim_amt), NA, ifelse(valid1$sim_diff == 1 | valid1$sim_diff == 4, 0, 1)))	# create column with domain, non-domain, and unrated questions
preds_trials1$sim_preds <- sim_l_pred1														# append similarity predictions
preds_trials1$choice <- valid1$choice															# append observed choices

# Select questions within similarity domain
subjects_df1 <- data.frame(subject = good_subjects1)	# create data frame of good subjects
domain1 <- subset(preds_trials1, sim_domain == 1)			# select domain questions
domain_sim1 <- aggregate(domain1$sim_l, by = list(domain1$subject), FUN = c("mean", "sum", "length"), na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(domain_sim1) <- c("subject", "acc", "numcorr", "numquestions")
domain_sim1 <- merge(subjects_df1, domain_sim1, all.x = TRUE) # merge with subject list to include subjects with no choices in similarity domain
domain_models_subj1 <- aggregate(domain1[, c(2:14)], by = list(domain1$subject), FUN = c("mean"), na.rm = TRUE)		# aggregate similarity predictive accuracy by subject
names(domain_models_subj1)[1] <- "subject"			
domain_models1 <- colMeans(domain_models_subj1[, -1])	# calculate mean predictive accuracy for models in domain

# Select questions outside of similarity domain
bothsim1 <- subset(preds_trials1, sim_diff == 1)  									# select both similar questions
bothsim1_sim <- aggregate(bothsim1$choice, by = list(bothsim1$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothsim1_sim) <- c("subject", "acc", "numLL", "numquestions")
bothsim1_subj <- factor(bothsim1_sim$subject)
bothsim1_sim <- merge(subjects_df1, bothsim1_sim, all.x = TRUE)			# merge subjects data frame with both_sim data
bothsim1_sim$obs <- bothsim1_sim$numLL / bothsim1_sim$numquestions 	# calculate observed proportion choice for LL
bothsim1_sim$deviation <- abs(bothsim1_sim$obs - 0.5)								# calculate deviation from chance (0.5)
bothsim1_sim$acc <- bothsim1_sim$deviation / 0.5										# calculate predictive accuracy for non-domain questions
bothdis1 <- subset(preds_trials1, sim_diff == 4) 										# select both dissimilar questions
bothdis1_sim <- aggregate(bothdis1$choice, by = list(bothdis1$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothdis1_sim) <- c("subject", "acc", "numLL", "numquestions")
bothdis1_sim$obs <- bothdis1_sim$numLL / bothdis1_sim$numquestions 	# calculate observed proportion choice for LL
bothdis1_sim$deviation <- abs(bothdis1_sim$obs - 0.5)								# calculate deviation from chance (0.5)
bothdis1_sim$acc <- bothdis1_sim$deviation / 0.5										# calculate predictive accuracy for non-domain questions

# Create overall predictive accuracy for similarity (Leland)
preds_all1 <- preds_subj1					# copy preds_subj1 to preds_all1
for(i in 1:num_good_subjects1) {	# calculate predictive accuracy as weighted average of accuracy within similarity domain and deviation from chance outside of domain
  if(preds_all1$subject[i] %in% bothsim1_subj) {	# if subjects are in the list of subjects with both similarity questions
    preds_all1$sim_l[i] <- (domain_sim1$numquestions[i] * (domain_sim1$numcorr[i] / domain_sim1$numquestions[i]) + bothsim1_sim$numquestions[i] * (1 - (abs((bothsim1_sim$numLL[i] / bothsim1_sim$numquestions[i]) - 0.5) / 0.5)) + bothdis1_sim$numquestions[i] * (1 - (abs((bothdis1_sim$numLL[i] / bothdis1_sim$numquestions[i]) - 0.5) / 0.5))) / (domain_sim1$numquestions[i] + bothsim1_sim$numquestions[i] + bothdis1_sim$numquestions[i])		# calculate Leland's similarity model predictive accuracy
  }
  else {	# if subjects are NOT in the list of subjects with both similarity questions
    preds_all1$sim_l[i] <- (domain_sim1$numquestions[i] * (domain_sim1$numcorr[i] / domain_sim1$numquestions[i]) + bothdis1_sim$numquestions[i] * (1 - (abs((bothdis1_sim$numLL[i] / bothdis1_sim$numquestions[i]) - 0.5) / 0.5))) / (domain_sim1$numquestions[i] + bothdis1_sim$numquestions[i])		# calculate Leland's similarity model predictive accuracy without "both similar" questions
  }
}
preds_all1$sim_l[34] <- (bothdis1_sim$numquestions[34] * (1 - (abs((bothdis1_sim$numLL[34] / bothdis1_sim$numquestions[34]) - 0.5) / 0.5))) / bothdis1_sim$numquestions[34] # calculate Leland's similarity model predictive accuracy only for outside domain for subject with NA within similarity domain

parameters_medians1 <- apply(parameters1, 2, FUN = median)	# calculate median parameter values

## Aggregate analysis
accuracy_all1 <- colMeans(preds_all1[, 2:(total_num_models1 + 1)], na.rm = TRUE) * 100   # calculate mean predictive accuracies for all models in all questions
accuracy_all1 <- accuracy_all1[c(1:6, 13, 7:12)]

# Find number of trials with no delay similarity judgments
num_trials1 <- length(valid1[, 1])																		# find total number of trials
num_rated_trials1 <-length(subset(valid1, !is.na(sim_domain))[, 1])		# find number of trials in which there was no rating of delay similarity
valid1a <- subset(validation1, subject %in% good_subjects1)  	# create subset of questions 
num_questions_all1 <- length(unique(valid1a$question))				# find total of questions
num_questions_dom1 <- length(unique(valid1$question))					# find number of questions with similarity judgments
simpercent1 <- num_questions_dom1 / num_questions_all1 * 100	# calculate percent of trials with similarity judgments

## Plot boxplot of predictive accuracy for each model
# Prepare data for boxplot
pred_acc_subj1 <- stack(preds_all1[c(2:14)])  				# reshape data to long/stacked format
names(pred_acc_subj1) <- c("pred_acc", "model")
pred_acc_subj1$model <- factor(pred_acc_subj1$model, levels = c("expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_l", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith"))	# reorder model levels
model_names_pred <- c("Exponential", "Hyperbolic (Mazur)", "Hyperbolic (Rachlin)", "Hyperbolic (Kirby)", "Hyperbolic (L&P)", "Arithmetic", "Similarity (Leland)", "Similarity+exponential", "Similarity+Mazur", "Similarity+Rachlin", "Similarity+Kirby", "Similarity+L&P", "Similarity+arithmetic")
pred_acc_subj1$model <- factor(pred_acc_subj1$model, labels = model_names_pred)	# rename model labels
pred_acc_subj1$subject <- rep(unique(preds_all1$subject), total_num_models1)		# create column of subject numbers

pred_acc_means1 <- aggregate(pred_acc_subj1$pred_acc, by=list(pred_acc_subj1$model), FUN = "mean", na.rm = TRUE)	# aggregate predictive accuracy by subject
names(pred_acc_means1) <- c("model", "pred_acc")
pred_acc_mean1 <- normDataWithin(data = pred_acc_subj1, idvar = "subject", measurevar="pred_acc", na.rm = TRUE)		# calculate normalized data for within-subjects CIs
pred_acc_norm1 <- summarySEwithin(data = pred_acc_mean1, measurevar="pred_accNormed", withinvars = "model", idvar = "subject", na.rm = TRUE)	# calculate within-subjects CIs
pred_acc_norm1$pred_acc <- pred_acc_means1$pred_acc									# copy over actual means
pred_acc_norm1$lci <- pred_acc_norm1$pred_acc - pred_acc_norm1$ci  	# add CIs
pred_acc_norm1$uci <- pred_acc_norm1$pred_acc + pred_acc_norm1$ci		# add CIs

white2blue <- colorRampPalette(c("white", "#145078"))(7)			# create color ramp from white to blue
box_colors <- c(white2blue[2:7], white2blue)						# duplicate color ramp for boxplot boxes

# Plot predictive accuracy boxplots for each model
pred_acc_bw <- bwplot(pred_acc * 100 ~ model, data = pred_acc_subj1,
  aspect = 0.8, coef = 0, ylab = "Predictive accuracy (%)", ylim = c(-5, 105),
  scales = list(x = list(rot = 90)), subscripts = T,
  par.settings = list(axis.text = list(cex = 1.8), par.ylab.text = list(cex = 2.2), 
    layout.heights = list(strip = 1.8), box.umbrella = list(lty = 1, col = "black", lwd = 1), 
    box.rectangle = list(lwd = 1, col = "black")),
  panel = function(x, y, ...) {
    x2 <- tapply(as.numeric(x), x, mean)
    mean.values <<- tapply(y, x, mean, na.rm = TRUE)						# calculates means
    panel.bwplot(x, y, pch = "|", horizontal = F, coef = 0,			# generate boxplot
      fill = box_colors)
    panel.points(mean.values, pch = 18, cex = 1.35, col = "black")	# plot means as diamonds
    panel.abline(h = max(mean.values), lty = 2)
    panel.abline(v = 6.5, lty = 1)  								# plot separation line between similarity and non-similarity models
    panel.text(x = 3.5, y = 0, "Discounting models", cex = 1.8)
    panel.text(x = 10, y = 0, "Similarity models", cex = 1.8)
  }
)
addWithinCI <- layer(panel.segments(x0 = model, y0 = lci[subscripts] * 100, x1 = model, y1 = uci[subscripts] * 100, col = "black", length = 0, unit = "native", angle = 0, code = 3, lwd = 2), data = pred_acc_norm1)	# create layer with within-subject CIs
cairo_pdf("figures/pred_acc1.pdf", width = 9.5, height = 11)
plot(pred_acc_bw + addWithinCI)
dev.off()

## Test difference between Mazur and similarity+Mazur
hyper_sim_diff_subj1 <- preds_all1$sim_h_m - preds_all1$hyper_m												# calculate difference between predictive accuracy of Mazur and similarity+Mazur
hyper_sim_diff_ci1 <- ci(hyper_sim_diff_subj1)																				# calculate mean, sd, CIs of difference
hyper_sim_diff_ci1$ci <- hyper_sim_diff_ci1$mean - hyper_sim_diff_ci1$lower95ci				# find CI difference
hyper_sim_diff_ci1$effect_size_d <-  hyper_sim_diff_ci1$mean / hyper_sim_diff_ci1$sd	# calculate Cohen's d effect size

## Find the number of questions in the similarity domain for each subject
sim_domain <- subset(valid1, sim_diff == 2 | sim_diff == 3)										# select data in the similarity domain
sim_no_domain <- subset(valid1, sim_diff == 1 | sim_diff == 4)								# select data outside the similarity domain
sim_domain_subj <- aggregate(choice ~ subject, sim_domain, length)						# find number of questions in the similarity domain for each subject
sim_domain_subj$pchoice <- sim_domain_subj$choice / num_questions_dom1 * 100	# calculate the proportion of choices in the similarity domain for each subject
sim_domain_subj <- merge(sim_domain_subj, data.frame(subject = as.factor(good_subjects1)), by = "subject", all.y = TRUE)	# merge the similarity domain data with the good subjects
sim_domain_mean <- mean(sim_domain_subj$pchoice, na.rm = TRUE)								# calculate mean proportion of choices in the similarity domain
sim_domain_median <- median(sim_domain_subj$pchoice, na.rm = TRUE)						# calculate median proportion of choices in the similarity domain
sim_domain_min <- min(sim_domain_subj$pchoice, na.rm = TRUE)									# calculate minimum proportion of choices in the similarity domain
sim_domain_max <- max(sim_domain_subj$pchoice, na.rm = TRUE)									# calculate maximum proportion of choices in the similarity domain

###############
## Magnitude effect
###############
## Prepare magnitude effect data
magnitude <- validation1[validation1$question >= 21 & validation1$question < 39, ]													# extract subset of magnitude effect questions (21-38)
magnitude$amt_ratio <- as.factor(ifelse(magnitude$question > 32, 1, ifelse(magnitude$question > 26, 2, 3)))	# assign amount ratio blocks to 1, 2, or 3

amt_ratiolab <- c("Amount ratio=0.80", "Amount ratio=0.67", "Amount ratio=0.50")			# create amount ratio labels
magnitude$immed <- as.factor(ifelse(magnitude$short_delay == 0, 1, 0))								# code short delays of 0 as immediate and > 0 as delayed
immlab <- c("Immediate", "Delayed")																										# create short delay labels
magnitude <- magnitude[order(magnitude$subject, magnitude$question), ]								# reorder data by subject and question
magnitude$adiff_num <- as.factor(rep(rep(1:3, 6), length(unique(magnitude$subject))))	# code amount differences (1, 2, 3)
magnitude$choice100 <- magnitude$choice * 100																					# convert choice to percentage

magnitude_means <- aggregate(magnitude$choice * 100, by = list(magnitude$large_amount, magnitude$immed, magnitude$amt_ratio), FUN = "mean")	# aggregate choice by large_amount, immed, and amt_ratio
names(magnitude_means) <- c("large_amount", "immed", "amt_ratio", "choice")
magnitude_mean <- normDataWithin(data = magnitude, idvar = "subject", measurevar="choice100")	# calculate normalized data for within-subjects CIs
magnitude_mean$large_amount <- as.factor(magnitude_mean$large_amount) 	# convert large_amount to factor
magnitude_norm <- summarySEwithin(data = magnitude_mean, measurevar="choice100Normed", withinvars = c("amt_ratio", "immed", "large_amount"), idvar = "subject")	# calculate within-subjects CIs
magnitude_norm$large_amount1 <- as.numeric(levels(magnitude_norm$large_amount)[as.integer(magnitude_norm$large_amount)])	# convert large_amount from factor to numeric
magnitude_norm$choice <- magnitude_means$choice													# copy actual means
magnitude_norm$uci <- magnitude_norm$choice + magnitude_norm$ci					# add CIs
magnitude_norm$lci <- magnitude_norm$choice - magnitude_norm$ci					# add CIs

## Plot magnitude effect choice data
my_pch <- c(16, 17)												# assign values for plot characters
my_col <- c("#0072B2", "#D55E00")					# assign values for plot colors
magnitude_xy <- xYplot(Cbind(choice, lci, uci) ~ large_amount1 | factor(amt_ratio, levels = levels(amt_ratio), labels = amt_ratiolab),
  groups = immed, data = magnitude_norm, type = "b", layout = c(3, 1), aspect = 1.8,
  xlab = "Large amount magnitude", ylim = c(-5, 105), ylab = "Percent choosing larger, later option",
  label.curves = F, between = list(x = 0.5), cex = 1.5, pch = my_pch,
  lty = c(1, 2), lwd = 3, col = my_col, par.strip.text = list(cex = 1.5),
  scales = list(x = list(relation ="free", at = list(c(5, 10, 15), c(3, 6, 12), c(2, 10, 18)))),
  par.settings = list(axis.text = list(cex = 1.4), par.ylab.text = list(cex = 1.8), par.xlab.text = list(cex = 1.8)), 
  key = list(text = list(immlab), cex = 1.5, corner = c(0, 0.9), points = list(pch = my_pch, col = my_col))
)
cairo_pdf(file = "figures/magnitude_tests1.pdf", width = 8, height = 5.5)
plot(magnitude_xy)
dev.off()

## Calculate similarity judgments for magnitude effect data
magnitude_sim_m <- aggregate(magnitude$sim_amt, by = list(magnitude$adiff_num), FUN = c("mean", "sd", "length"))	# aggregate amount similarity judgments by amount difference
names(magnitude_sim_m) <- c("large_amount", "sim_diff", "sd", "N")
magnitude_sim <- aggregate(magnitude$sim_amt * 100, by = list(magnitude$large_amount, magnitude$immed, magnitude$amt_ratio), FUN = "mean")	# aggregate amount similarity judgments by amount difference, immediacy, and k values
names(magnitude_sim) <- c("large_amount", "immed", "amt_ratio", "sim_amt")
magnitude$sim_amt100 <- magnitude$sim_amt * 100		# convert choice to percentage

magnitude_sim_mean <- normDataWithin(data = magnitude, idvar = "subject", measurevar="sim_amt100")	# calculate normalized data for within-subjects CIs
magnitude_sim_mean$large_amount <- as.factor(magnitude_sim_mean$large_amount) 											# convert large_amount to factor
magnitude_sim_norm <- summarySEwithin(data = magnitude_sim_mean, measurevar="sim_amt100Normed", withinvars = c("amt_ratio", "immed", "large_amount"), idvar = "subject")	# calculate within-subjects CIs
magnitude_sim_norm$large_amount1 <- as.numeric(levels(magnitude_sim_norm$large_amount)[as.integer(magnitude_sim_norm$large_amount)])	# convert large_amount from factor to numeric
magnitude_sim_norm$sim_amt <- magnitude_sim$sim_amt														# copy actual means
magnitude_sim_norm$uci <- magnitude_sim_norm$sim_amt + magnitude_sim_norm$ci	# add CIs
magnitude_sim_norm$lci <- magnitude_sim_norm$sim_amt - magnitude_sim_norm$ci	# add CIs

## Plot magnitude effect similarity judgment data
magnitude_sim <- xYplot(Cbind(sim_amt, lci, uci) ~ large_amount1 | factor(amt_ratio, levels = levels(amt_ratio), labels = amt_ratiolab),
  data = subset(magnitude_sim_norm, immed == 0), type = "b", layout = c(3, 1), aspect = 1.8,
  xlab = "Large amount magnitude", ylim = c(-5, 105), ylab = "Percent judging amounts as similar",
  label.curves = F, between = list(x = 0.5), cex = 1.5, pch = my_pch,
  lty = c(1, 2), lwd = 3, col = my_col, par.strip.text = list(cex = 1.5),
  scales = list(x = list(relation ="free", at = list(c(5, 10, 15), c(3, 6, 12), c(2, 10, 18)))),
  par.settings = list(axis.text = list(cex = 1.4), par.ylab.text = list(cex = 1.8), par.xlab.text = list(cex = 1.8)), 
)
cairo_pdf(file = "figures/magnitude_tests_sim1.pdf", width = 8, height = 5.5)
plot(magnitude_sim)
dev.off()

##############################
## Experiment 2--Delay magnitude effect
##############################

################
## Load and prepare data
################

# Load all data
all_data3 <- read.csv("stevens_dataS2.csv")  # load data
all_data3 <- subset(all_data3, question != 26)    		# remove question accidentally included

subject_nums3 <- unique(all_data3$subject)	# find unique subject numbers
num_subjects3 <- length(subject_nums3)  		# find number of subjects

## Staircase (training data)
staircase3 <- subset(all_data3, type == "staircase" & small_amount != 10)	# subset staircase data with adjusting amount and with different small and large amounts (small_amount != 10)

# Plot staircase data
staircase_means3 <- aggregate(choice ~ long_delay * short_delay * subject, staircase3, mean)		# aggregate staircase choice data by long delay, short delay, and subject
staircase_means3$short_delay <- as.factor(staircase_means3$short_delay)													# convert short_delay to factor
col.blind <- c("#0072B2", "#D55E00", "#009E73")																									# assign color-blind safe colors for plots
staircase_plot3 <- xyplot(choice * 100 ~ long_delay | as.factor(subject), data = staircase_means3, 
  xlab = "Long delay (days)", ylab = "Mean percent choice for LL", type = "b", as.table = TRUE,
  par.settings = list(par.xlab.text = list(cex = 1.65), par.ylab.text = list(cex = 1.65))
)
cairo_pdf(file = "figures/staircase_plot2.pdf", width = 8, height = 11)										# create PDF
plot(staircase_plot3)																																			# print plot
dev.off()																																												# turn off PDF device

## Similarity data
sim_data3 <- subset(all_data3, type == "amount" | type == "delay") 	# subset all similarity questions
names(sim_data3)[8] <- "similarity"																	# rename similarity column
sim_amt3 <- subset(sim_data3, type == "amount")											# subset similarity judgments for amount
sim_amt3 <- sim_amt3[, c(1, 4, 5, 8)]																# extract relevant columns
sim_delay3 <- subset(sim_data3, type == "delay")										# subset similarity judgments for delay
sim_delay3 <- sim_delay3[, c(1, 6, 7, 8)]														# extract relevant columns

## Binary (validation) data
binary3 <- subset(all_data3, type == "binary")																															# subset binary choice data
validation3 <- merge(binary3, sim_amt3, by = c("subject", "small_amount", "large_amount"), all.x = TRUE)		# merge binary choice data with similarity amount judgment data
names(validation3)[9] <- "sim_amt"																																					# rename sim_amt column
validation3 <- merge(validation3, sim_delay3, by = c("subject", "short_delay", "long_delay"), all.x = TRUE)	# merge binary choice data with similarity delay judgment data
names(validation3)[10] <- "sim_delay"																																				# rename sim_amt column
validation3$sim_diff <- ifelse(validation3$sim_amt == 1 & validation3$sim_delay == 1, 1, ifelse(validation3$sim_amt == 1 & validation3$sim_delay == 0, 2, ifelse(validation3$sim_amt == 0 & validation3$sim_delay == 1, 3, ifelse(validation3$sim_amt == 0 & validation3$sim_delay == 0, 4, NA))))  # categorize difference in similarity judgments (1 = both similar, 2 = amounts similar, delays dissimilar, 3 = amounts dissimilar, delays similar, 4 = both dissimilar)
validation3$k <- (validation3$small_amount - validation3$large_amount)/(-validation3$small_amount * validation3$long_delay + validation3$large_amount * validation3$short_delay)		# calculate k value

###############
## Participant info
###############
## Age data experiment 3
age3 <- c(18, 21, 19, 19, 19, 18, 20, 20, 21, 19, 21, 20, 19, 19, 20, 18, 18, 20, 22, 18, 19, 21, 22, 21, 18, 19, 21, 20, 18, 19, 22, 18, 19, 20, 21, 20, 22, 20, 19, 21, 18, 18, 19, 20, 18, 19, 45, 21, 20, 20, 19, 20, 22, 23, 20, 20, 21, 20, 20, 19, 18, 18) 
age3_m <- mean(age3)	# calculate mean age
age3_sd <- sd(age3)		# calculate standard deviation of age
age3_min <- min(age3)	# calculate minimum age
age3_max <- max(age3)	# calculate maximum age

## Sex data experiment 3
sex3 <- c(1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2, 1)
males3 <- table(sex3)[1]
females3 <- table(sex3)[2]

###############
## Model selection
###############

####
## Fit staircase (training) data with maximum likelihood estimation (MLE)
####

## Conduct subject-wise maximum likelihood estimation
# Prepare vectors and data frames
good_subjects3 <- NA # initiate vector of 'good' subjects (subject for whom the MLE converges)

fitsNAs <- rep(NA, num_subjects3)	# initiate vector of NAs to fill empty data frame
mle_fits3 <- data.frame(subject = subject_nums3, expon_AICc = fitsNAs, delta = fitsNAs, hyper_m_AICc = fitsNAs, k_m = fitsNAs, hyper_r_AICc = fitsNAs, k_r = fitsNAs, sigma = fitsNAs, hyper_k_AICc = fitsNAs, k_k = fitsNAs, mu = fitsNAs, hyper_lp_AICc = fitsNAs, alpha = fitsNAs, beta = fitsNAs, add_AICc = fitsNAs, lambda = fitsNAs, expon_diffs = fitsNAs, hyper_m_diffs = fitsNAs, hyper_r_diffs = fitsNAs, hyper_k_diffs = fitsNAs, hyper_lp_diffs = fitsNAs, add_diffs = fitsNAs)	# initiate data frame for MLE analysis

# Conduct MLE for each subject
foreach(subj = subject_nums3) %do% {								# for each subject
  curr_subj <- subset(staircase3, subject == subj)	# assign curr_subj to the current subject's data
  
  # Calculate MLE for each subject
  if(sum(curr_subj$choice) < 70 & subj != 6 & subj != 17 & subj != 56) {  	# remove subjects whose MLE do not converge
    good_subjects3 <- c(good_subjects3, subj)  			# append current subject to good_subjects3
    expon_fit <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = curr_subj, start = list(delta = 0.05))			# fit exponential discounting model
    if(subj == 10 | subj == 20 | subj == 21 | subj == 39 | subj == 41 | subj == 44 | subj == 61) {  # for some subjects, use different starting parameters
      hyper_m_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = curr_subj, start = list(k = 0.5))  				# fit hyperbolic (Mazur) discounting model
    } else {
      hyper_m_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = curr_subj, start = list(k = 0.05))					# fit hyperbolic (Mazur) discounting model
    }
    if(subj == 6) {    						# if subject = 6, use different starting parameters
      hyper_r_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, s) - hyper_r(small_amount, short_delay, k, s)), size = 1), data = curr_subj, start = list(k = 0.05, s = 0.05))  # fit hyperbolic (Rachlin) discounting model
    } else {
      hyper_r_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, s) - hyper_r(small_amount, short_delay, k, s)), size = 1), data = curr_subj, start = list(k = 0.05, s = 0.5))  # fit hyperbolic (Rachlin) discounting model
    }
    if(subj == 9) {    						# if subject = 9, use different starting parameters
      hyper_k_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, m) - hyper_k(small_amount, short_delay, k, m)), size = 1), data = curr_subj, start = list(k = 0.5, m = -0.25))  # fit hyperbolic (Kirby) discounting model      
    } else {
      hyper_k_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, m) - hyper_k(small_amount, short_delay, k, m)), size = 1), data = curr_subj, start = list(k = 0.05, m = -0.25))  # fit hyperbolic (Kirby) discounting model
    }
    if(subj == 32 | subj == 44) {	# if subject = 32 or 44, use different starting parameters
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.5, beta = 0.05))				# fit hyperbolic (L&P) discounting model
    } else if(subj == 41 | subj == 47 | subj == 51 | subj == 61) {  			# if subject = 41, 47, 51, or 61, use different starting parameters
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 5, beta = 0.05))    		# fit hyperbolic (L&P) discounting model
    } else if(subj == 57) {  			# if subject = 57, use different starting parameters
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.05), method = "Nelder-Mead")    		# fit hyperbolic (L&P) discounting model
    } else {
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.05))			# fit hyperbolic (L&P) discounting model
    }
    add_fit <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = curr_subj, start = list(lambda = 0.05))	# fit additive discounting model
    
    # Calculate AICc for each model and extract parameters    
    mle_fits3$expon_AICc[which(mle_fits3$subject == subj)[1]] <- expon_aic <- AICc(expon_fit, nobs = 1, k = 2)					# calculate AICc for exponential discounting
    mle_fits3$delta[which(mle_fits3$subject == subj)[1]] <- coef(expon_fit)																							# extract delta
    
    mle_fits3$hyper_m_AICc[which(mle_fits3$subject == subj)[1]] <- hyper_m_aic <- AICc(hyper_m_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur) discounting
    mle_fits3$k_m[which(mle_fits3$subject == subj)[1]] <- coef(hyper_m_fit)																							# extract k
    
    mle_fits3$hyper_r_AICc[which(mle_fits3$subject == subj)[1]] <- hyper_r_aic <- AICc(hyper_r_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Rachlin) discounting
    mle_fits3$k_r[which(mle_fits3$subject == subj)[1]] <- coef(hyper_r_fit)[1]																					# extract k
    mle_fits3$sigma[which(mle_fits3$subject == subj)[1]] <- coef(hyper_r_fit)[2]																				# extract sigma
    
    mle_fits3$hyper_k_AICc[which(mle_fits3$subject == subj)[1]] <- hyper_k_aic <- AICc(hyper_k_fit, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby) discounting
    mle_fits3$k_k[which(mle_fits3$subject == subj)[1]] <- coef(hyper_k_fit)[1]																					# extract k
    mle_fits3$mu[which(mle_fits3$subject == subj)[1]] <- coef(hyper_k_fit)[2]																						# extract mu
    
    mle_fits3$hyper_lp_AICc[which(mle_fits3$subject == subj)[1]] <- hyper_lp_aic <- AICc(hyper_lp_fit, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P) discounting
    mle_fits3$alpha[which(mle_fits3$subject == subj)[1]] <- coef(hyper_lp_fit)[1]																				# extract alpha
    mle_fits3$beta[which(mle_fits3$subject == subj)[1]] <- coef(hyper_lp_fit)[2]																				# extract beta
    
    mle_fits3$add_AICc[which(mle_fits3$subject == subj)[1]] <- add_aic <- AICc(add_fit, nobs = 1, k = 2)								# calculate AICc for additive discounting
    mle_fits3$lambda[which(mle_fits3$subject == subj)[1]] <- coef(add_fit)																							# extract lambda
    
    subj_aics <- c(expon_aic, hyper_m_aic, hyper_r_aic, hyper_k_aic, hyper_lp_aic, add_aic)															# concatenate AICc values for all models
    aic_diffs <- subj_aics - min(subj_aics)		# calculate AICc differences
    aic_weights <- (exp(-0.5 * aic_diffs)) / (sum(exp(-0.5 * aic_diffs)))	# calculate AICc weights
    mle_fits3[which(mle_fits3$subject == subj)[1], 17:22] <- aic_diffs		# assign AICc differences to mle_first3
  }
}

good_subjects3 <- good_subjects3[-1]					# remove initial NA
num_good_subjects3 <- length(good_subjects3)	# find number of good subjects

aic_individual_median3 <- apply(mle_fits3[, c(2, 4, 6, 9, 12, 15)], 2, median, na.rm = TRUE)  # calculate median AICc values across subjects
aic_diff_individual_median3 <- apply(mle_fits3[, c(17:22)], 2, median, na.rm = TRUE)  				# calculate median AICc differences across subjects

parameters3 <- mle_fits3[, c(3, 5, 7, 8, 10, 11, 13, 14, 16)]	# extract parameters for each subject and model
parameters3 <- subset(parameters3, !is.na(delta))							# remove subjects with no parameter estimates

## Conduct MLE for aggregate data
train3 <- subset(staircase3, subject %in% good_subjects3)	# subset staircase data from 'good' subjects

expon_fit_group3 <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = train3, start = list(delta = 0.05))	# calculate MLE for exponential discounting model

hyper_m_fit_group3 <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = train3, start = list(k = 0.05))	# calculate MLE for hyperbolic (Mazur) discounting model

hyper_r_fit_group3 <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, sigma) - hyper_r(small_amount, short_delay, k, sigma)), size = 1), data = train3, start = list(k = 0.05, sigma = 0.5))  # calculate MLE for hyperbolic (Rachlin) discounting model

hyper_k_fit_group3 <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, mu) - hyper_k(small_amount, short_delay, k, mu)), size = 1), data = train3, start = list(k = 0.05, mu = -0.5))  # calculate MLE for hyperbolic (Kirby) discounting model

hyper_lp_fit_group3 <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = train3, start = list(alpha = 0.05, beta = 0.05))	# calculate MLE for hyperbolic (L&P) discounting model

add_fit_group3 <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = train3, start = list(lambda = 0.05))	# calculate MLE for exponential discounting model

expon_aic_group3 <- AICc(expon_fit_group3, nobs = 1, k = 2)				# calculate AICc for exponential discounting
hyper_m_aic_group3 <- AICc(hyper_m_fit_group3, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur)  discounting
hyper_r_aic_group3 <- AICc(hyper_r_fit_group3, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Rachlin)  discounting
hyper_k_aic_group3 <- AICc(hyper_k_fit_group3, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby)  discounting
hyper_lp_aic_group3 <- AICc(hyper_lp_fit_group3, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P)  discounting
add_aic_group3 <- AICc(add_fit_group3, nobs = 1, k = 2)						# calculate AICc for additive discounting

aics_group3 <- c(expon_aic_group3, hyper_m_aic_group3, hyper_r_aic_group3, hyper_k_aic_group3, hyper_lp_aic_group3, add_aic_group3)	# concatenate AICc values for all models
aic_diffs_group3 <- aics_group3 - min(aics_group3)																																									# calculate AICc differences
aic_weights_group3 <- (exp(-0.5 * aic_diffs_group3)) / (sum(exp(-0.5 * aic_diffs_group3)))																					# calculate AICc weights
evidence_ratios_group3 <- c(aic_weights_group3[3] / aic_weights_group3[1], aic_weights_group3[3] / aic_weights_group3[2], aic_weights_group3[3] / aic_weights_group3[4], aic_weights_group3[3] / aic_weights_group3[5], aic_weights_group3[3] / aic_weights_group3[6])																																# calculate evidence ratios

####
## Predict binary choice (validation) data using parameters fit from staircase (training) data
####

## Prepare data
validation3$sim_domain <- ifelse(is.na(validation3$sim_delay), NA, ifelse(is.na(validation3$sim_amt), NA, ifelse(validation3$sim_diff == 1 | validation3$sim_diff == 4, 0, 1))) 	# create column with domain, non-domain, and unrated questions
valid3 <- subset(validation3, subject %in% good_subjects3 & !is.na(sim_diff) & small_amount == 7)		# create subset of questions from good subjects (with variance in choices) that have similarity judgments
valid3$subject <- factor(valid3$subject)

# Establish number of models and parameters
num_models3 <- 6  												# assign total number of models fitted (not including similarity models)
total_num_models3 <- 2 * num_models3 + 1	# include similarity models in total number of models
num_params3 <- dim(parameters3)[2]				# assign total number of parameters estimated

## Create xtable of questions and mean responses
valid3_aggregation <- aggregate(valid3$choice, by = list(valid3$small_amount, valid3$large_amount, valid3$short_delay, valid3$long_delay, valid3$question, valid3$k), FUN = c("mean", "sd")) 												# calculate mean and sd for each question
valid3_aggregation <- valid3_aggregation[, -5]    							# remove question column
valid3_aggregation <- valid3_aggregation[, c(3, 4, 1, 2, 5:7)] # reorder columns with delays first
names(valid3_aggregation) <- c("Short delay", "Long delay", "Small amount", "Large amount", "k", "Mean choice for LL", "Standard deviation")
valid3_aggregation <- valid3_aggregation[order(valid3_aggregation$"k", valid3_aggregation$"Short delay", valid3_aggregation$"Long delay", valid3_aggregation$"Small amount", valid3_aggregation$"Large amount", valid3_aggregation$"Mean choice for LL"), ]																				# reorder rows
valid3_xtable <- xtable(valid3_aggregation, caption = 'Questions and Mean Responses for Experiment 2', label = "tab:EXP2QUEST")	# create xtable for import into LaTeX

## Predict binary choice data (generalization/validation)
# Prepare data frames
preds_trials3 <- data.frame(valid3$subject, matrix(rep(NA, length(valid3$subject) * (total_num_models3 - 1)), nrow = length(valid3$subject)))	# create data frame of NAs for trial predictions
names(preds_trials3) <- c("subject", "expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith")
preds_subj3 <- data.frame(good_subjects3, matrix(rep(NA, num_good_subjects3 * (total_num_models3 - 1)), nrow = num_good_subjects3))						# create data frame of NAs for subject predictions
names(preds_subj3) <- names(preds_trials3) 
num_questions3 <- length(unique(valid3$question))	# find the number of unique questions

# Exponential discounting model
valid3$delta <- rep(parameters3$delta, each = num_questions3)   	# create column of consecutives instances of each subject's fitted parameter for each question
exp_pred3 <- ifelse(valid3$small_amount * exp(-valid3$delta * valid3$short_delay) < valid3$large_amount * exp(-valid3$delta * valid3$long_delay), 1, ifelse(valid3$small_amount * exp(-valid3$delta * valid3$short_delay) > valid3$large_amount * exp(-valid3$delta * valid3$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials3$expo <- ifelse(exp_pred3 == valid3$choice, 1, 0)
preds_subj3$expo <- aggregate(preds_trials3$expo, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Hyperbolic discounting model Mazur
valid3$k_m <- rep(parameters3$k_m, each = num_questions3)
hyp_m_pred3 <- ifelse(valid3$small_amount / (1 + valid3$k_m * valid3$short_delay) < valid3$large_amount / (1 + valid3$k_m * valid3$long_delay), 1, ifelse(valid3$small_amount / (1 + valid3$k_m * valid3$short_delay) > valid3$large_amount / (1 + valid3$k_m * valid3$long_delay), 0, NA))
preds_trials3$hyper_m <- ifelse(hyp_m_pred3 == valid3$choice, 1, 0)
preds_subj3$hyper_m <- aggregate(preds_trials3$hyper_m, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Hyperbolic discounting model Rachlin
valid3$k_r <- rep(parameters3$k_r, each = num_questions3)
valid3$sigma <- rep(parameters3$sigma, each = num_questions3)
hyp_r_pred3 <- ifelse(valid3$small_amount / (1 + valid3$k_r * valid3$short_delay ^ valid3$sigma)  < valid3$large_amount / (1 + valid3$k_r * valid3$long_delay ^ valid3$sigma), 1, ifelse(valid3$small_amount / (1 + valid3$k_r * valid3$short_delay ^ valid3$sigma) > valid3$large_amount / (1 + valid3$k_r * valid3$long_delay ^ valid3$sigma), 0, NA))
preds_trials3$hyper_r <- ifelse(hyp_r_pred3 == valid3$choice, 1, 0)
preds_subj3$hyper_r <- aggregate(preds_trials3$hyper_r, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Hyperbolic discounting model Kirby
valid3$k_k <- rep(parameters3$k_k, each = num_questions3)
valid3$mu <- rep(parameters3$mu, each = num_questions3)
hyp_k_pred3 <- ifelse(valid3$small_amount / (1 + valid3$k_k * valid3$small_amount ^ valid3$mu * valid3$short_delay) < valid3$large_amount / (1 + valid3$k_k * valid3$large_amount ^ valid3$mu * valid3$long_delay), 1, ifelse(valid3$small_amount / (1 + valid3$k_k * valid3$small_amount ^ valid3$mu * valid3$short_delay) > valid3$large_amount / (1 + valid3$k_k * valid3$large_amount ^ valid3$mu * valid3$long_delay), 0, NA))
preds_trials3$hyper_k <- ifelse(hyp_k_pred3 == valid3$choice, 1, 0)
preds_subj3$hyper_k <- aggregate(preds_trials3$hyper_k, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Hyperbolic discounting model L&P
valid3$alpha <- rep(parameters3$alpha, each = num_questions3)
valid3$beta <- rep(parameters3$beta, each = num_questions3)
hyp_lp_pred3 <- ifelse(valid3$small_amount / (1 + valid3$alpha * valid3$short_delay) ^ (valid3$beta / valid3$alpha)  < valid3$large_amount / (1 + valid3$alpha * valid3$long_delay) ^ (valid3$beta / valid3$alpha), 1, ifelse(valid3$small_amount / (1 + valid3$alpha * valid3$short_delay) ^ (valid3$beta / valid3$alpha) > valid3$large_amount / (1 + valid3$alpha * valid3$long_delay) ^ (valid3$beta / valid3$alpha), 0, NA))
preds_trials3$hyper_lp <- ifelse(hyp_lp_pred3 == valid3$choice, 1, 0)
preds_subj3$hyper_lp <- aggregate(preds_trials3$hyper_lp, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Arithmetic model
valid3$lambda <- rep(parameters3$lambda, each = num_questions3)
arith_pred3 <- ifelse(valid3$small_amount - valid3$lambda * valid3$short_delay < valid3$large_amount - valid3$lambda * valid3$long_delay, 1,  ifelse(valid3$small_amount - valid3$lambda - valid3$short_delay > valid3$large_amount - valid3$lambda - valid3$long_delay, 0, NA))
preds_trials3$arith <- ifelse(arith_pred3 == valid3$choice, 1, 0)
preds_subj3$arith <- aggregate(preds_trials3$arith, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (exponential)--first look for similarity differences, then use exponential discounting if similarity doesn't distinguish
sim_e_pred3 <- ifelse(valid3$sim_diff == 2, 0, ifelse(valid3$sim_diff == 3, 1, exp_pred3))
preds_trials3$sim_e <- ifelse(sim_e_pred3 == valid3$choice, 1, 0)
preds_subj3$sim_e <- aggregate(preds_trials3$sim_e, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--Mazur)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hm_pred3 <- ifelse(valid3$sim_diff == 2, 0,	ifelse(valid3$sim_diff == 3, 1, hyp_m_pred3))
preds_trials3$sim_h_m <- ifelse(sim_hm_pred3 == valid3$choice, 1, 0)
preds_subj3$sim_h_m <- aggregate(preds_trials3$sim_h_m, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--Rachlin)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hr_pred3 <- ifelse(valid3$sim_diff == 2, 0,  ifelse(valid3$sim_diff == 3, 1, hyp_r_pred3))
preds_trials3$sim_h_r <- ifelse(sim_hr_pred3 == valid3$choice, 1, 0)
preds_subj3$sim_h_r <- aggregate(preds_trials3$sim_h_r, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--Kirby)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hk_pred3 <- ifelse(valid3$sim_diff == 2, 0,  ifelse(valid3$sim_diff == 3, 1, hyp_k_pred3))
preds_trials3$sim_h_k <- ifelse(sim_hk_pred3 == valid3$choice, 1, 0)
preds_subj3$sim_h_k <- aggregate(preds_trials3$sim_h_k, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (hyperbolic--L&P)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hlp_pred3 <- ifelse(valid3$sim_diff == 2, 0,  ifelse(valid3$sim_diff == 3, 1, hyp_lp_pred3))
preds_trials3$sim_h_lp <- ifelse(sim_hlp_pred3 == valid3$choice, 1, 0)
preds_subj3$sim_h_lp <- aggregate(preds_trials3$sim_h_lp, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (Arithmetic)--first look for similarity differences, then use Arithmetic discounting if similarity doesn't distinguish
sim_a_pred3 <- ifelse(valid3$sim_diff == 2, 0, ifelse(valid3$sim_diff == 3, 1, arith_pred3))
preds_trials3$sim_arith <- ifelse(sim_a_pred3 == valid3$choice, 1, 0)
preds_subj3$sim_arith <- aggregate(preds_trials3$sim_arith, by = list(preds_trials3$subject), FUN = mean)[, 2]

# Similarity model (Leland)
sim_l_pred3 <- ifelse(valid3$sim_diff == 2, 0, ifelse(valid3$sim_diff == 3, 1, NA))	# create array differentiating between higher sim for amount, higher sim for delay, and equal sim

# Add Leland similarity model and other columns to preds_trials
preds_trials3$sim_l <- ifelse(sim_l_pred3 == valid3$choice, 1, 0)	# create column of correct predictions for similarity
preds_trials3$sim_diff <- valid3$sim_diff													# append similarity difference values
preds_trials3$sim_domain <- ifelse(is.na(valid3$sim_delay), NA, ifelse(is.na(valid3$sim_amt), NA, ifelse(valid3$sim_diff == 1 | valid3$sim_diff == 4, 0, 1)))	# create column with domain, non-domain, and unrated questions
preds_trials3$sim_preds <- sim_l_pred3														# append similarity predictions
preds_trials3$choice <- valid3$choice															# append observed choices

# Select questions within similarity domain
subjects_df3 <- data.frame(subject = good_subjects3)	# create data frame of good subjects
domain3 <- subset(preds_trials3, sim_domain == 1)			# select domain questions
domain_sim3 <- aggregate(domain3$sim_l, by = list(domain3$subject), FUN = c("mean", "sum", "length"), na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(domain_sim3) <- c("subject", "acc", "numcorr", "numquestions")
domain_sim3 <- merge(subjects_df3, domain_sim3, all.x = TRUE) # merge with subject list to include subjects with no choices in similarity domain
domain_models_subj3 <- aggregate(domain3[, c(2:14)], by = list(domain3$subject), FUN = c("mean"), na.rm = TRUE)		# aggregate similarity predictive accuracy by subject
names(domain_models_subj3)[1] <- "subject"			
domain_models3 <- colMeans(domain_models_subj3[, -1])	# calculate mean predictive accuracy for models in domain

# Select questions outside of similarity domain
bothsim3 <- subset(preds_trials3, sim_diff == 1)  									# select both similar questions
bothsim3_sim <- aggregate(bothsim3$choice, by = list(bothsim3$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothsim3_sim) <- c("subject", "acc", "numLL", "numquestions")
bothsim3_subj <- factor(bothsim3_sim$subject)
bothsim3_sim <- merge(subjects_df3, bothsim3_sim, all.x = TRUE)			# merge subjects data frame with both_sim data
bothsim3_sim$obs <- bothsim3_sim$numLL / bothsim3_sim$numquestions 	# calculate observed proportion choice for LL
bothsim3_sim$deviation <- abs(bothsim3_sim$obs - 0.5)								# calculate deviation from chance (0.5)
bothsim3_sim$acc <- bothsim3_sim$deviation / 0.5										# calculate predictive accuracy for non-domain questions
bothdis3 <- subset(preds_trials3, sim_diff == 4) 										# select both dissimilar questions
bothdis3_sim <- aggregate(bothdis3$choice, by = list(bothdis3$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothdis3_sim) <- c("subject", "acc", "numLL", "numquestions")
bothdis3_sim$obs <- bothdis3_sim$numLL / bothdis3_sim$numquestions 	# calculate observed proportion choice for LL
bothdis3_sim$deviation <- abs(bothdis3_sim$obs - 0.5)								# calculate deviation from chance (0.5)
bothdis3_sim$acc <- bothdis3_sim$deviation / 0.5										# calculate predictive accuracy for non-domain questions

# Create overall predictive accuracy for similarity (Leland)
preds_all3 <- preds_subj3					# copy preds_subj3 to preds_all3
for(i in 1:num_good_subjects3) {	# calculate predictive accuracy as weighted average of accuracy within similarity domain and deviation from chance outside of domain
  if(preds_all3$subject[i] %in% bothsim3_subj) {	# if subjects are in the list of subjects with both similarity questions
    preds_all3$sim_l[i] <- (domain_sim3$numquestions[i] * (domain_sim3$numcorr[i] / domain_sim3$numquestions[i]) + bothsim3_sim$numquestions[i] * (1 - (abs((bothsim3_sim$numLL[i] / bothsim3_sim$numquestions[i]) - 0.5) / 0.5)) + bothdis3_sim$numquestions[i] * (1 - (abs((bothdis3_sim$numLL[i] / bothdis3_sim$numquestions[i]) - 0.5) / 0.5))) / (domain_sim3$numquestions[i] + bothsim3_sim$numquestions[i] + bothdis3_sim$numquestions[i])		# calculate Leland's similarity model predictive accuracy
  }
  else {	# if subjects are NOT in the list of subjects with both similarity questions
    preds_all3$sim_l[i] <- (domain_sim3$numquestions[i] * (domain_sim3$numcorr[i] / domain_sim3$numquestions[i]) + bothdis3_sim$numquestions[i] * (1 - (abs((bothdis3_sim$numLL[i] / bothdis3_sim$numquestions[i]) - 0.5) / 0.5))) / (domain_sim3$numquestions[i] + bothdis3_sim$numquestions[i])		# calculate Leland's similarity model predictive accuracy without "both similar" questions
  }
}

parameters_medians3 <- apply(parameters3, 2, FUN = median)	# calculate median parameter values

## Aggregate analysis
accuracy_all3 <- colMeans(preds_all3[, 2:(total_num_models3 + 1)], na.rm = TRUE) * 100   # calculate mean predictive accuracies for all models in all questions
accuracy_all3 <- accuracy_all3[c(1:6, 13, 7:12)]

# Find number of trials with no delay similarity judgments
num_trials3 <- length(valid3[, 1])																		# find total number of trials
num_rated_trials3 <-length(subset(valid3, !is.na(sim_domain))[, 1])		# find number of trials in which there was no rating of delay similarity
valid3a <- subset(validation3, subject %in% good_subjects3)  	# create subset of questions 
num_questions_all3 <- length(unique(valid3a$question))				# find total of questions
num_questions_dom3 <- length(unique(valid3$question))					# find number of questions with similarity judgments
simpercent3 <- num_questions_dom3 / num_questions_all3 * 100	# calculate percent of trials with similarity judgments

## Plot boxplot of predictive accuracy for each model
# Prepare data for boxplot
pred_acc_subj3 <- stack(preds_all3[c(2:14)])  				# reshape data to long/stacked format
names(pred_acc_subj3) <- c("pred_acc", "model")
pred_acc_subj3$model <- factor(pred_acc_subj3$model, levels = c("expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_l", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith"))	# reorder model levels
model_names_pred <- c("Exponential", "Hyperbolic (Mazur)", "Hyperbolic (Rachlin)", "Hyperbolic (Kirby)", "Hyperbolic (L&P)", "Arithmetic", "Similarity (Leland)", "Similarity+exponential", "Similarity+Mazur", "Similarity+Rachlin", "Similarity+Kirby", "Similarity+L&P", "Similarity+arithmetic")
pred_acc_subj3$model <- factor(pred_acc_subj3$model, labels = model_names_pred)	# rename model labels
pred_acc_subj3$subject <- rep(unique(preds_all3$subject), total_num_models3)		# create column of subject numbers

pred_acc_means3 <- aggregate(pred_acc_subj3$pred_acc, by=list(pred_acc_subj3$model), FUN = "mean", na.rm = TRUE)	# aggregate predictive accuracy by subject
names(pred_acc_means3) <- c("model", "pred_acc")
pred_acc_mean3 <- normDataWithin(data = pred_acc_subj3, idvar = "subject", measurevar="pred_acc", na.rm = TRUE)		# calculate normalized data for within-subjects CIs
pred_acc_norm3 <- summarySEwithin(data = pred_acc_mean3, measurevar="pred_accNormed", withinvars = "model", idvar = "subject", na.rm = TRUE)	# calculate within-subjects CIs
pred_acc_norm3$pred_acc <- pred_acc_means3$pred_acc									# copy over actual means
pred_acc_norm3$lci <- pred_acc_norm3$pred_acc - pred_acc_norm3$ci  	# add CIs
pred_acc_norm3$uci <- pred_acc_norm3$pred_acc + pred_acc_norm3$ci		# add CIs

white2blue <- colorRampPalette(c("white", "#145078"))(7)			# create color ramp from white to blue
box_colors <- c(white2blue[2:7], white2blue)						# duplicate color ramp for boxplot boxes

# Plot predictive accuracy boxplots for each model
pred_acc_bw <- bwplot(pred_acc * 100 ~ model, data = pred_acc_subj3,
  aspect = 0.8, coef = 0, ylab = "Predictive accuracy (%)", ylim = c(-10, 105),
  scales = list(x = list(rot = 90)), subscripts = T,
  par.settings = list(axis.text = list(cex = 1.8), par.ylab.text = list(cex = 2.2), 
    layout.heights = list(strip = 1.8), box.umbrella = list(lty = 1, col = "black", lwd = 1), 
    box.rectangle = list(lwd = 1, col = "black")),
  panel = function(x, y, ...) {
    x2 <- tapply(as.numeric(x), x, mean)
    mean.values <<- tapply(y, x, mean, na.rm = TRUE)						# calculates means
    panel.bwplot(x, y, pch = "|", horizontal = F, coef = 0,	fill = box_colors)		# generate boxplot
    panel.points(mean.values, pch = 18, cex = 1.35, col = "black")	# plot means as diamonds
    panel.abline(h = max(mean.values), lty = 2)
    panel.abline(v = 6.5, lty = 1)  								# plot separation line between similarity and non-similarity models
    panel.text(x = 3.5, y = -5, "Discounting models", cex = 1.8)
    panel.text(x = 10, y = -5, "Similarity models", cex = 1.8)
  }
)
addWithinCI <- layer(panel.segments(x0 = model, y0 = lci[subscripts] * 100, x1 = model, y1 = uci[subscripts] * 100, col = "black", length = 0, unit = "native", angle = 0, code = 3, lwd = 2), data = pred_acc_norm3)	# create layer with within-subject CIs
cairo_pdf("figures/pred_acc2.pdf", width = 9.5, height = 11)
plot(pred_acc_bw + addWithinCI)
dev.off()

## Test difference between Mazur and similarity+Mazur
hyper_sim_diff_subj3 <- preds_all3$sim_h_m - preds_all3$hyper_m												# calculate difference between predictive accuracy of Mazur and similarity+Mazur
hyper_sim_diff_ci3 <- ci(hyper_sim_diff_subj3)																				# calculate mean, sd, CIs of difference
hyper_sim_diff_ci3$ci <- hyper_sim_diff_ci3$mean - hyper_sim_diff_ci3$lower95ci				# find CI difference
hyper_sim_diff_ci3$effect_size_d <-  hyper_sim_diff_ci3$mean / hyper_sim_diff_ci3$sd	# calculate Cohen's d effect size

## Find the number of questions in the similarity domain for each subject
sim_domain <- subset(valid3, sim_diff == 2 | sim_diff == 3)										# select data in the similarity domain
sim_no_domain <- subset(valid3, sim_diff == 1 | sim_diff == 4)								# select data outside the similarity domain
sim_domain_subj <- aggregate(choice ~ subject, sim_domain, length)						# find number of questions in the similarity domain for each subject
sim_domain_subj$pchoice <- sim_domain_subj$choice / num_questions_dom3 * 100	# calculate the proportion of choices in the similarity domain for each subject
sim_domain_subj <- merge(sim_domain_subj, data.frame(subject = as.factor(good_subjects3)), by = "subject", all.y = TRUE)	# merge the similarity domain data with the good subjects
sim_domain_mean <- mean(sim_domain_subj$pchoice, na.rm = TRUE)								# calculate mean proportion of choices in the similarity domain
sim_domain_median <- median(sim_domain_subj$pchoice, na.rm = TRUE)						# calculate median proportion of choices in the similarity domain
sim_domain_min <- min(sim_domain_subj$pchoice, na.rm = TRUE)									# calculate minimum proportion of choices in the similarity domain
sim_domain_max <- max(sim_domain_subj$pchoice, na.rm = TRUE)									# calculate maximum proportion of choices in the similarity domain

validation3$sim_pred <- ifelse(validation3$sim_diff == 2, 0, ifelse(validation3$sim_diff == 3, 1, NA))
validation3$sim_acc <- ifelse(is.na(validation3$sim_pred), NA, ifelse(validation3$choice == validation3$sim_pred, 1, 0))

###############
## Delay magnitude effect
###############
## Prepare delay magnitude effect data
magnitude3 <- subset(validation3, small_amount == 7)												# extract subset of magnitude effect questions (21-38)

klab <- c("k=0.333", "k=0.5", "k=0.6", "k=0.75", "k=1")			# create k labels
magnitude3$choice100 <- magnitude3$choice * 100																					# convert choice to percentage

magnitude_means3 <- aggregate(magnitude3$choice * 100, by = list(magnitude3$short_delay, magnitude3$k), FUN = "mean")	# aggregate choice by large_amount, immed, and amt_ratio
names(magnitude_means3) <- c("short_delay", "k", "choice")
magnitude_mean3 <- normDataWithin(data = magnitude3, idvar = "subject", measurevar="choice100")	# calculate normalized data for within-subjects CIs
magnitude_mean3$large_amount <- as.factor(magnitude_mean3$large_amount) 	# convert large_amount to factor
magnitude_norm3 <- summarySEwithin(data = magnitude_mean3, measurevar="choice100Normed", withinvars = c("k", "short_delay"), idvar = "subject")	# calculate within-subjects CIs
magnitude_norm3$short_delay3 <- as.numeric(levels(magnitude_norm3$short_delay)[as.integer(magnitude_norm3$short_delay)])	# convert short_delay from factor to numeric
magnitude_norm3$choice <- magnitude_means3$choice													# copy actual means
magnitude_norm3$uci <- magnitude_norm3$choice + magnitude_norm3$ci					# add CIs
magnitude_norm3$lci <- magnitude_norm3$choice - magnitude_norm3$ci					# add CIs

## Plot magnitude effect choice data
my_pch <- c(36, 17)												# assign values for plot characters
my_col <- c("#0072B2", "#D55E00")					# assign values for plot colors
magnitude_xy <- xYplot(Cbind(choice, lci, uci) ~ short_delay3 | factor(k, labels = klab),
  data = magnitude_norm3, type = "b", layout = c(5, 1), aspect = 1.8,
  xlab = "Short delay", ylab = "Percent choosing larger, later option", ylim = c(18, 88),
  label.curves = F, between = list(x = 0.2, y = 0.2), cex = 1.5, pch = 16,
  lty = c(1, 2), lwd = 3, col = my_col, par.strip.text = list(cex = 1.5),
  par.settings = list(axis.text = list(cex = 1.4), par.ylab.text = list(cex = 1.8), par.xlab.text = list(cex = 1.8)), 
)
cairo_pdf(file = "figures/magnitude_tests2.pdf", width = 12, height = 5.5)
plot(magnitude_xy)
dev.off()

## Calculate similarity judgments for magnitude effect data
magnitude_sim3 <- aggregate(magnitude3$sim_delay * 100, by = list(magnitude3$short_delay, magnitude3$k), FUN = "mean")	# aggregate amount similarity judgments by amount difference, immediacy, and k values
names(magnitude_sim3) <- c("short_delay", "k", "sim_delay")
magnitude3$sim_delay100 <- magnitude3$sim_delay * 100		# convert choice to percentage

magnitude_sim_mean3 <- normDataWithin(data = magnitude3, idvar = "subject", measurevar="sim_delay100")	# calculate normalized data for within-subjects CIs
magnitude_sim_mean3$short_delay <- as.factor(magnitude_sim_mean3$short_delay) 											# convert short_delay to factor
magnitude_sim_norm3 <- summarySEwithin(data = magnitude_sim_mean3, measurevar="sim_delay100Normed", withinvars = c("k", "short_delay"), idvar = "subject")	# calculate within-subjects CIs
magnitude_sim_norm3$short_delay3 <- as.numeric(levels(magnitude_sim_norm3$short_delay)[as.integer(magnitude_sim_norm3$short_delay)])	# convert large_amount from factor to numeric
magnitude_sim_norm3$sim_delay <- magnitude_sim3$sim_delay														# copy actual means
magnitude_sim_norm3$uci <- magnitude_sim_norm3$sim_delay + magnitude_sim_norm3$ci	# add CIs
magnitude_sim_norm3$lci <- magnitude_sim_norm3$sim_delay - magnitude_sim_norm3$ci	# add CIs

## Plot magnitude effect similarity judgment data
magnitude_sim3_plot <- xYplot(Cbind(sim_delay, lci, uci) ~ short_delay3 | factor(k, levels = levels(k), labels = klab),
  data = magnitude_sim_norm3, type = "b", layout = c(5, 1), aspect = 1.8,
  xlab = "Short delay", ylab = "Percent judging delays as similar", ylim = c(5, 80),
  label.curves = F, between = list(x = 0.2, y = 0.2), cex = 1.5, pch = 16,
  lty = c(1, 2), lwd = 3, col = my_col, par.strip.text = list(cex = 1.5),
  par.settings = list(axis.text = list(cex = 1.4), par.ylab.text = list(cex = 1.8), par.xlab.text = list(cex = 1.8)), 
)
cairo_pdf(file = "figures/magnitude_tests_sim2.pdf", width = 12, height = 5.5)
plot(magnitude_sim3_plot)
dev.off()

detach("package:plyr", unload = TRUE)

##############################
## Experiment 3--Gains and losses
##############################

################
## Load and prepare data
################

all_data2 <- read.csv("stevens_dataS3.csv")

subject_nums2 <- unique(all_data2$subject)	# find unique subject numbers
num_subjects2 <- length(subject_nums2)  		# find number of subjects

## Staircase (training data)
staircase2 <- subset(all_data2, type == "staircase" & small_amount != 10)		# subset staircase data EXCEPT when both small amount and large amount are $10
staircase_gain2 <- subset(staircase2, gainloss == "gain")										# subset staircase gains data
staircase_loss2 <- subset(staircase2, gainloss == "loss")										# subset staircase loss data

# Plot staircase data
staircase_means_gain2 <- aggregate(choice ~ long_delay * short_delay * subject, staircase_gain2, mean)	# calculate mean choice by delays and subject for gains
staircase_means_gain2$short_delay <- as.factor(staircase_means_gain2$short_delay)
staircase_plot_gain2 <- xyplot(choice * 100 ~ long_delay | as.factor(subject), group = short_delay, data = staircase_means_gain2, xlab = "Long delay (days)", ylab = "Mean percent choice for LL", type = "b", as.table = TRUE, par.settings = list(par.xlab.text = list(cex = 1.65), par.ylab.text = list(cex = 1.65))
)
cairo_pdf(file = "figures/staircase_plot_gain3.pdf", width = 8, height = 11)
plot(staircase_plot_gain2)
dev.off()

staircase_means_loss2 <- aggregate(choice ~ long_delay * short_delay * subject, staircase_loss2, mean)	# calculate mean choice by delays and subject for losses
staircase_means_loss2$short_delay <- as.factor(staircase_means_loss2$short_delay)
staircase_plot_loss2 <- xyplot(choice * 100 ~ long_delay | as.factor(subject), group = short_delay, data = staircase_means_loss2, xlab = "Long delay (days)", ylab = "Mean percent choice for LL", type = "b", as.table = TRUE, par.settings = list(par.xlab.text = list(cex = 1.65), par.ylab.text = list(cex = 1.65))
)
cairo_pdf(file = "figures/staircase_plot_loss3.pdf", width = 8, height = 11)
plot(staircase_plot_loss2)
dev.off()

# Determine direction of preference (advancing or delaying gains and losses)
staircase_equal_amt2 <- subset(all_data2, type == "staircase" & small_amount == 10)  # subset staircase data ONLY when both small amount and large amount are $10
preference_direction2 <- aggregate(choice ~ subject * gainloss, data = staircase_equal_amt2, FUN = "mean")
preference_direction2$direction <- ifelse(preference_direction2$choice > 0.5, "L", ifelse(preference_direction2$choice < 0.5, "S", "N"))


## Similarity data
sim_data2 <- subset(all_data2, type == "similarity") 												# select all similarity questions
sim_data2 <- subset(sim_data2, question != 35 & question != 36 & question != 41 & question != 83 & question != 84 & question != 89)  # remove duplicate question
sim_data2$attribute <- ifelse(sim_data2$question < 97, "amount", "delay")		# assign whether question is for amount or delay
names(sim_data2)[9] <- "similarity"
sim_amt2 <- subset(sim_data2, attribute == "amount")												# select amount judgment questions
sim_amt2 <- sim_amt2[, c(1, 3, 5, 6, 9)]																		# select amount judgment columns
sim_delay2 <- subset(sim_data2, attribute == "delay")												# select delay judgment questions
sim_delay2 <- sim_delay2[, c(1, 7, 8, 9)]																		# select delay judgment columns

## Binary (validation) data
binary2 <- subset(all_data2, type == "binary")					# select binary choice questions
validation2 <- merge(binary2, sim_amt2, by = c("subject", "small_amount", "large_amount", "gainloss"), all.x = TRUE)	# merge similarity amount judgments with binary choice data
validation2 <- validation2[!duplicated(validation2), ]	# remove duplicates
names(validation2)[10] <- "sim_amt"
validation2 <- merge(validation2, sim_delay2, by = c("subject", "short_delay", "long_delay"), all.x = TRUE)	# merge similarity delay judgments with binary choice data
validation2 <- validation2[!duplicated(validation2), ]	# remove duplicates
names(validation2)[11] <- "sim_delay"
validation2$sim_diff <- ifelse(validation2$sim_amt == 1 & validation2$sim_delay == 1, 1, ifelse(validation2$sim_amt == 1 & validation2$sim_delay == 0, 2, ifelse(validation2$sim_amt == 0 & validation2$sim_delay == 1, 3, ifelse(validation2$sim_amt == 0 & validation2$sim_delay == 0, 4, NA))))  # categorize difference in similarity judgments (1 = both similar, 2 = amounts similar, delays dissimilar, 3 = amounts dissimilar, delays similar, 4 = both dissimilar)
validation2 <- validation2[order(validation2$subject, validation2$question),]	# reorder validation data

###############
## Participant info
###############

age2 <- c(19, 19, 18, 19, 18, 20, 19, 20, 21, 21, 18, 25, 20, 39, 18, 22, 18, 19, 18, 21, 22, 18, 18, 18, 19, 18, 22, 18, 18, 20, 20, 19, 19, 19, 18, 18, 19, 17, 22, 18, 19, 23, 20, 21, 19, 18, 19, 21, 18, 18, 22, 19, 20, 18, 20, 20, 19, 19, 18, 20, 18, 21, 19, 20, 21, 21, 21, 22)
age2_m <- mean(age2)		# calculate mean age
age2_sd <- sd(age2)			# calculate sd of age
age2_min <- min(age2)		# calculate minimum age
age2_max <- max(age2)		# calculate maximum age

################
## Calculate model fit
################

#######
## Gains data
#######
# Prepare vectors and data frames for subject-wise maximum likelihood estimation (MLE)
good_subjects_gain2 <- NA																# initiate vector of 'good' subjects (subject for whom the MLE converges)
num_good_subjects_gain2 <- length(good_subjects_gain2)	# find number of good subjects
fitsNAs <- rep(NA, num_subjects2)												# initiate vector of NAs to fill empty data frame
mle_fits_gain2 <- data.frame(subject = subject_nums2, expon_AICc = fitsNAs, delta = fitsNAs, hyper_m_AICc = fitsNAs, k_m = fitsNAs, hyper_r_AICc = fitsNAs, k_r = fitsNAs, sigma = fitsNAs, hyper_k_AICc = fitsNAs, k_k = fitsNAs, mu = fitsNAs, hyper_lp_AICc = fitsNAs, alpha = fitsNAs, beta = fitsNAs, add_AICc = fitsNAs, lambda = fitsNAs, expon_diffs = fitsNAs, hyper_m_diffs = fitsNAs, hyper_r_diffs = fitsNAs, hyper_k_diffs = fitsNAs, hyper_lp_diffs = fitsNAs, add_diffs = fitsNAs)	# initiate data frame for MLE analysis

# Conduct MLE for each subject
foreach(subj = subject_nums2) %do% {														# for each subject
  curr_subj <- subset(staircase_gain2, subject == subj)					# assign curr_subj to the current subject's data
 # Calculate MLE for each model
  if(sum(curr_subj$choice) > 1 & sum(curr_subj$choice) < 51) {  # remove subjects whose MLE do not converge
    good_subjects_gain2 <- c(good_subjects_gain2, subj)  					# append current subject as good subject
    expon_fit <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = curr_subj, start = list(delta = 0.05))			# fit exponential discounting model
    hyper_m_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = curr_subj, start = list(k = 0.05))					# fit hyperbolic (Mazur) discounting model
    hyper_r_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, sigma) - hyper_r(small_amount, short_delay, k, sigma)), size = 1), data = curr_subj, start = list(k = 0.05, sigma = 0.5))	# fit hyperbolic (Rachlin) discounting model
    hyper_k_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, mu) - hyper_k(small_amount, short_delay, k, mu)), size = 1), data = curr_subj, start = list(k = 0.05, mu = -0.5))  # fit hyperbolic (Kirby) discounting model
    if(subj %in% c(3)) {
    hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.005, beta = 0.5))  			# fit hyperbolic (L&P) discounting model
  } else if(subj %in% c(18, 52)){
    hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.05))  			# fit hyperbolic (L&P) discounting model
  } else if(subj %in% c(21, 31, 36, 41, 54)){
    hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.5, beta = 0.005))    		# fit hyperbolic (L&P) discounting model
  } else {
    hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.5))				# fit hyperbolic (L&P) discounting model
}
    add_fit <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = curr_subj, start = list(lambda = 0.05))	# fit additive discounting model
    
# Calculate AICc for each model and extract parameters    
    mle_fits_gain2$expon_AICc[which(mle_fits_gain2$subject == subj)[1]] <- expon_aic <- AICc(expon_fit, nobs = 1, k = 2)					# calculate AICc for exponential discounting
    mle_fits_gain2$delta[which(mle_fits_gain2$subject == subj)[1]] <- coef(expon_fit)																							# extract delta
    
    mle_fits_gain2$hyper_m_AICc[which(mle_fits_gain2$subject == subj)[1]] <- hyper_m_aic <- AICc(hyper_m_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur) discounting
    mle_fits_gain2$k_m[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_m_fit)																							# extract k
    
    mle_fits_gain2$hyper_r_AICc[which(mle_fits_gain2$subject == subj)[1]] <- hyper_r_aic <- AICc(hyper_r_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Rachlin) discounting
    mle_fits_gain2$k_r[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_r_fit)[1]																					# extract k
    mle_fits_gain2$sigma[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_r_fit)[2]																				# extract sigma
    
    mle_fits_gain2$hyper_k_AICc[which(mle_fits_gain2$subject == subj)[1]] <- hyper_k_aic <- AICc(hyper_k_fit, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby) discounting
    mle_fits_gain2$k_k[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_k_fit)[1]																					# extract k
    mle_fits_gain2$mu[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_k_fit)[2]																						# extract mu

    mle_fits_gain2$hyper_lp_AICc[which(mle_fits_gain2$subject == subj)[1]] <- hyper_lp_aic <- AICc(hyper_lp_fit, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P) discounting
    mle_fits_gain2$alpha[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_lp_fit)[1]																				# extract alpha
    mle_fits_gain2$beta[which(mle_fits_gain2$subject == subj)[1]] <- coef(hyper_lp_fit)[2]																				# extract beta
    
    mle_fits_gain2$add_AICc[which(mle_fits_gain2$subject == subj)[1]] <- add_aic <- AICc(add_fit, nobs = 1, k = 2)								# calculate AICc for additive discounting
    mle_fits_gain2$lambda[which(mle_fits_gain2$subject == subj)[1]] <- coef(add_fit)																							# extract lambda
    
    subj_aics <- c(expon_aic, hyper_m_aic, hyper_r_aic, hyper_k_aic, hyper_lp_aic, add_aic)		# concatenate AICc values for all models
    aic_diffs <- subj_aics - min(subj_aics)																										# calculate AICc differences
    mle_fits_gain2[which(mle_fits_gain2$subject == subj)[1], 17:22] <- aic_diffs							# assign AICc differences to mle_first_gain2
  }
}

good_subjects_gain2 <- good_subjects_gain2[-1]					# remove initial NA
num_good_subjects_gain2 <- length(good_subjects_gain2)	# find number of good subjects

aic_individual_median_gain2 <- apply(mle_fits_gain2[, c(2, 4, 6, 9, 12, 15)], 2, median, na.rm = TRUE)  # calculate median AICc values across subjects
aic_diff_individual_median_gain2 <- apply(mle_fits_gain2[, c(17:22)], 2, median, na.rm = TRUE)  		# calculate median AICc differences across subjects

parameters_gain2 <- mle_fits_gain2[, c(1, 3, 5, 7, 8, 10, 11, 13, 14, 16)]	# extract parameters for each subject and model
parameters_gain2 <- subset(parameters_gain2, !is.na(delta))									# remove subjects with no parameter estimates

# Conduct MLE for aggregate data
train_gain2 <- subset(staircase_gain2, subject %in% good_subjects_gain2)	# subset staircase data from 'good' subjects

expon_fit_group_gain2 <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = train_gain2, start = list(delta = 0.05))	# calculate MLE for exponential discounting model
expon_par_group_gain2 <- coef(expon_fit_group_gain2)			# extract delta
expon_ci_group_gain2 <- confint(expon_fit_group_gain2)  	# calculate delta CIs
expon_ci_diff_group_gain2 <- expon_par_group_gain2 - expon_ci_group_gain2[1]  	# calculate delta CIs

hyper_m_fit_group_gain2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = train_gain2, start = list(k = 0.05))	# calculate MLE for hyperbolic (Mazur) discounting model

hyper_r_fit_group_gain2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, s) - hyper_r(small_amount, short_delay, k, s)), size = 1), data = train_gain2, start = list(k = 0.05, s = 0.5))	# calculate MLE for hyperbolic (Rachlin) discounting model

hyper_k_fit_group_gain2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, mu) - hyper_k(small_amount, short_delay, k, mu)), size = 1), data = train_gain2, start = list(k = 0.05, mu = -0.5))  # calculate MLE for hyperbolic (Kirby) discounting model

hyper_lp_fit_group_gain2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = train_gain2, start = list(alpha = 0.05, beta = 0.05))	# calculate MLE for hyperbolic (L&P) discounting model

add_fit_group_gain2 <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = train_gain2, start = list(lambda = 0.05))	# calculate MLE for exponential discounting model

expon_aic_group_gain2 <- AICc(expon_fit_group_gain2, nobs = 1, k = 2)				# calculate AICc for exponential discounting
hyper_m_aic_group_gain2 <- AICc(hyper_m_fit_group_gain2, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur) discounting
hyper_r_aic_group_gain2 <- AICc(hyper_r_fit_group_gain2, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Rachlin) discounting
hyper_k_aic_group_gain2 <- AICc(hyper_k_fit_group_gain2, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Rachlin) discounting
hyper_lp_aic_group_gain2 <- AICc(hyper_lp_fit_group_gain2, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P) discounting
add_aic_group_gain2 <- AICc(add_fit_group_gain2, nobs = 1, k = 2)						# calculate AICc for additive discounting

aics_group_gain2 <- c(expon_aic_group_gain2, hyper_m_aic_group_gain2, hyper_r_aic_group_gain2, hyper_k_aic_group_gain2, hyper_lp_aic_group_gain2, add_aic_group_gain2)	# concatenate AICc values for all models
aic_diffs_group_gain2 <- aics_group_gain2 - min(aics_group_gain2)						# calculate AICc differences
aic_weights_group_gain2 <- (exp(-0.5 * aic_diffs_group_gain2)) / (sum(exp(-0.5 * aic_diffs_group_gain2)))	# calculate AICc weights
evidence_ratios_group_gain2 <- c(aic_weights_group_gain2[3] / aic_weights_group_gain2[1], aic_weights_group_gain2[3] / aic_weights_group_gain2[2], aic_weights_group_gain2[3] / aic_weights_group_gain2[4], aic_weights_group_gain2[3] / aic_weights_group_gain2[5], aic_weights_group_gain2[3] / aic_weights_group_gain2[6])	# calculate evidence ratios

#######
## Losses data
#######
# Prepare vectors and data frames for subject-wise maximum likelihood estimation (MLE)
staircase_loss2$small_amount <- -staircase_loss2$small_amount	# make small amounts negative
staircase_loss2$large_amount <- -staircase_loss2$large_amount	# make large amounts negative
good_subjects_loss2 <- NA 																		# create vector of 'good' subjects (subject for whom the MLE converges)
mle_fits_loss2 <- data.frame(subject = subject_nums2, expon_AICc = fitsNAs, delta = fitsNAs, hyper_m_AICc = fitsNAs, k_m = fitsNAs, hyper_r_AICc = fitsNAs, k_r = fitsNAs, sigma = fitsNAs, hyper_k_AICc = fitsNAs, k_k = fitsNAs, mu = fitsNAs, hyper_lp_AICc = fitsNAs, alpha = fitsNAs, beta = fitsNAs, add_AICc = fitsNAs, lambda = fitsNAs, expon_diffs = fitsNAs, hyper_m_diffs = fitsNAs, hyper_r_diffs = fitsNAs, hyper_k_diffs = fitsNAs, hyper_lp_diffs = fitsNAs, add_diffs = fitsNAs)	# initiate data frame for MLE analysis

# Conduct MLE for each subject
foreach(subj = subject_nums2) %do% {				# for each subject
  curr_subj <- subset(staircase_loss2, subject == subj)	# assign curr_subj to the current subject's data
# Calculate MLE for each model
  if(sum(curr_subj$choice) > 1 & sum(curr_subj$choice) < 53) {	# remove subjects whose MLE do not converge
    good_subjects_loss2 <- c(good_subjects_loss2, subj) 				# append current subject to good subjects
    expon_fit <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = curr_subj, start = list(delta = 0.05))			# fit exponential discounting model
    hyper_m_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = curr_subj, start = list(k = 0.05))					# fit hyperbolic (Mazur) discounting model
    hyper_r_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, sigma) - hyper_r(small_amount, short_delay, k, sigma)), size = 1), data = curr_subj, start = list(k = 0.05, sigma = 0.5))  # fit hyperbolic (Rachlin) discounting model
    hyper_k_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, mu) - hyper_k(small_amount, short_delay, k, mu)), size = 1), data = curr_subj, start = list(k = 0.05, mu = -0.5))  # fit hyperbolic (Kirby) discounting model
    if(subj == 49) {	# if subject = 25 or 46, use different starting parameters
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.005, beta = 0.05))				# fit hyperbolic (L&P) discounting model
    } else {
      hyper_lp_fit <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = curr_subj, start = list(alpha = 0.05, beta = 0.05))				# fit hyperbolic (L&P) discounting model
    }
    add_fit <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = curr_subj, start = list(lambda = 0.05))	# fit additive discounting model
    
# Calculate AICc for each model and extract parameters    
    mle_fits_loss2$expon_AICc[which(mle_fits_loss2$subject == subj)[1]] <- expon_aic <- AICc(expon_fit, nobs = 1, k = 2)					# calculate AICc for exponential discounting
    mle_fits_loss2$delta[which(mle_fits_loss2$subject == subj)[1]] <- coef(expon_fit)																							# extract delta
    
    mle_fits_loss2$hyper_m_AICc[which(mle_fits_loss2$subject == subj)[1]] <- hyper_m_aic <- AICc(hyper_m_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur) discounting
    mle_fits_loss2$k_m[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_m_fit)																							# extract k
    
    mle_fits_loss2$hyper_r_AICc[which(mle_fits_loss2$subject == subj)[1]] <- hyper_r_aic <- AICc(hyper_r_fit, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Rachlin) discounting
    mle_fits_loss2$k_r[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_r_fit)[1]																					# extract k
    mle_fits_loss2$sigma[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_r_fit)[2]																				# extract sigma
    
    mle_fits_loss2$hyper_k_AICc[which(mle_fits_loss2$subject == subj)[1]] <- hyper_k_aic <- AICc(hyper_k_fit, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby) discounting
    mle_fits_loss2$k_k[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_k_fit)[1]																					# extract k
    mle_fits_loss2$mu[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_k_fit)[2]																						# extract mu

    mle_fits_loss2$hyper_lp_AICc[which(mle_fits_loss2$subject == subj)[1]] <- hyper_lp_aic <- AICc(hyper_lp_fit, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P) discounting
    mle_fits_loss2$alpha[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_lp_fit)[1]																				# extract alpha
    mle_fits_loss2$beta[which(mle_fits_loss2$subject == subj)[1]] <- coef(hyper_lp_fit)[2]																				# extract beta
    
    mle_fits_loss2$add_AICc[which(mle_fits_loss2$subject == subj)[1]] <- add_aic <- AICc(add_fit, nobs = 1, k = 2)								# calculate AICc for additive discounting
    mle_fits_loss2$lambda[which(mle_fits_loss2$subject == subj)[1]] <- coef(add_fit)																							# extract lambda
    
    subj_aics <- c(expon_aic, hyper_m_aic, hyper_r_aic, hyper_k_aic, hyper_lp_aic, add_aic)		# concatenate AICc values for all models
    aic_diffs <- subj_aics - min(subj_aics)																										# calculate AICc differences
    mle_fits_loss2[which(mle_fits_loss2$subject == subj)[1], 17:22] <- aic_diffs							# assign AICc differences to mle_first_loss2
  }
}
good_subjects_loss2 <- good_subjects_loss2[-1]					# remove initial NA
num_good_subjects_loss2 <- length(good_subjects_loss2)	# find number of good subjects

aic_individual_median_loss2 <- apply(mle_fits_loss2[, c(2, 4, 6, 9, 12, 15)], 2, median, na.rm = TRUE)  # calculate median AICc values across subjects
aic_diff_individual_median_loss2 <- apply(mle_fits_loss2[, c(17:22)], 2, median, na.rm = TRUE) 					# calculate median AICc differences across subjects

parameters_loss2 <- mle_fits_loss2[, c(1, 3, 5, 7, 8, 10, 11, 13, 14, 16)]	# extract parameters for each subject and model
parameters_loss2 <- subset(parameters_loss2, !is.na(delta))									# remove subjects with no parameter estimates

# Conduct MLE for aggregate data
train_loss2 <- subset(staircase_loss2, subject %in% good_subjects_loss2)	# subset staircase data from 'good' subjects

expon_fit_group_loss2 <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = train_loss2, start = list(delta = 0.05))	# calculate MLE for exponential discounting model
expon_par_group_loss2 <- coef(expon_fit_group_loss2)				# extract delta
expon_ci_group_loss2 <- confint(expon_fit_group_loss2)			# calculate delta CIs
expon_ci_diff_group_loss2 <- expon_par_group_loss2 - expon_ci_group_loss2[1]    # calculate delta CIs

hyper_m_fit_group_loss2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_m(large_amount, long_delay, k) - hyper_m(small_amount, short_delay, k)), size = 1), data = train_loss2, start = list(k = 0.05))	# calculate MLE for hyperbolic (Mazur) discounting model

hyper_r_fit_group_loss2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_r(large_amount, long_delay, k, s) - hyper_r(small_amount, short_delay, k, s)), size = 1), data = train_loss2, start = list(k = 0.05, s = 0.5))	# calculate MLE for hyperbolic (Rachlin) discounting model

hyper_k_fit_group_loss2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_k(large_amount, long_delay, k, mu) - hyper_k(small_amount, short_delay, k, mu)), size = 1), data = train_loss2, start = list(k = 0.05, mu = -0.5))  # calculate MLE for hyperbolic (Kirby) discounting model

hyper_lp_fit_group_loss2 <- mle2(choice ~ dbinom(prob = invlogit(hyper_lp(large_amount, long_delay, alpha, beta) - hyper_lp(small_amount, short_delay, alpha, beta)), size = 1), data = train_loss2, start = list(alpha = 0.05, beta = 0.05))	# calculate MLE for hyperbolic (L&P) discounting model

add_fit_group_loss2 <- mle2(choice ~ dbinom(prob = invlogit(additive(large_amount, long_delay, lambda) - additive(small_amount, short_delay, lambda)), size = 1), data = train_loss2, start = list(lambda = 0.05))	# calculate MLE for exponential discounting model

expon_aic_group_loss2 <- AICc(expon_fit_group_loss2, nobs = 1, k = 2)				# calculate AICc for exponential discounting
hyper_m_aic_group_loss2 <- AICc(hyper_m_fit_group_loss2, nobs = 1, k = 2)		# calculate AICc for hyperbolic (Mazur)  discounting
hyper_r_aic_group_loss2 <- AICc(hyper_r_fit_group_loss2, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Rachlin)  discounting
hyper_k_aic_group_loss2 <- AICc(hyper_k_fit_group_loss2, nobs = 1, k = 2)  	# calculate AICc for hyperbolic (Kirby)  discounting
hyper_lp_aic_group_loss2 <- AICc(hyper_lp_fit_group_loss2, nobs = 1, k = 2)	# calculate AICc for hyperbolic (L&P)  discounting
add_aic_group_loss2 <- AICc(add_fit_group_loss2, nobs = 1, k = 2)						# calculate AICc for additive discounting

aics_group_loss2 <- c(expon_aic_group_loss2, hyper_m_aic_group_loss2, hyper_r_aic_group_loss2, hyper_k_aic_group_loss2, hyper_lp_aic_group_loss2, add_aic_group_loss2)	# concatenate AICc values for all models
aic_diffs_group_loss2 <- aics_group_loss2 - min(aics_group_loss2)	# calculate AICc differences
aic_weights_group_loss2 <- (exp(-0.5 * aic_diffs_group_loss2)) / (sum(exp(-0.5 * aic_diffs_group_loss2)))	# calculate AICc weights
evidence_ratios_group_loss2 <- c(aic_weights_group_loss2[3] / aic_weights_group_loss2[1], aic_weights_group_loss2[3] / aic_weights_group_loss2[2], aic_weights_group_loss2[3] / aic_weights_group_loss2[4], aic_weights_group_loss2[3] / aic_weights_group_loss2[5], aic_weights_group_loss2[3] / aic_weights_group_loss2[6])			# calculate evidence ratios

###################
## Calculate predictive accuracy
###################

#####
## Calculate predictive accuracy for gain data
#####
## Prepare data
validation_gain2 <- subset(validation2, gainloss == "gain" & subject %in% good_subjects_gain2) # extract gain data for subjects with complete training data
num_models2 <- 6    											# assign total number of models fitted (not including similarity models)
total_num_models2 <- 2 * num_models2 + 1	# include similarity models in total number of models
num_params2 <- dim(parameters_gain2)[2]		# assign total number of parameters estimated

## Create data frame of questions and mean responses
validation_gain2_aggregation <- aggregate(validation_gain2$choice, by = list(validation_gain2$small_amount, validation_gain2$large_amount, validation_gain2$short_delay, validation_gain2$long_delay, validation_gain2$question), FUN = c("mean", "sd")) # calculate mean and sd for each question
validation_gain2_aggregation <- validation_gain2_aggregation[, -5]  								# remove question column
validation_gain2_aggregation <- validation_gain2_aggregation[, c(3, 4, 1, 2, 5, 6)] # reorder columns with delays first
names(validation_gain2_aggregation) <- c("Short delay", "Long delay", "Small amount", "Large amount", "Mean choice for LL", "Standard deviation")
validation_gain2_aggregation <- validation_gain2_aggregation[order(validation_gain2_aggregation$"Short delay", validation_gain2_aggregation$"Long delay", validation_gain2_aggregation$"Small amount", validation_gain2_aggregation$"Large amount"), ]																			# sort columns by delays and amounts

## Similarity data in similarity domain
validation_gain2$sim_pred <- ifelse(validation_gain2$sim_diff == 2, 0, ifelse(validation_gain2$sim_diff == 3, 1, NA)) 	# generate predictions when similarity can make a prediction (only one of two attributes is similar)
validation_gain2$sim_acc <- ifelse(!is.na(validation_gain2$sim_pred), ifelse(validation_gain2$choice == validation_gain2$sim_pred, 1, 0), NA) # assess whether prediction matches observed choice
validation_gain2 <- validation_gain2[order(validation_gain2$subject, validation_gain2$question),]	# sort by subject and question

library(plyr)

## Predict binary choice data (generalization/validation)
preds_trials_gain2 <- data.frame(validation_gain2$subject, matrix(rep(NA, length(validation_gain2$subject) * (total_num_models2 - 1)), nrow = length(validation_gain2$subject)))  # create data frame of NAs for trial predictions
names(preds_trials_gain2) <- c("subject", "expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith")
preds_subj_gain2 <- data.frame(good_subjects_gain2, matrix(rep(NA, num_good_subjects_gain2 * (total_num_models2 - 1)), nrow = num_good_subjects_gain2))	# create data frame of NAs for subject predictions
names(preds_subj_gain2) <- names(preds_trials_gain2) 
column_num1 <- column_num2 <- 2															# assign initial column number
num_questions2 <- length(unique(validation_gain2$question))	# find number of questions

# Exponential discounting model
validation_gain2$delta <- rep(parameters_gain2$delta, each = num_questions2)   	# create column of each subject's fitted parameter
exp_pred_gain2 <- ifelse(validation_gain2$small_amount * exp(-validation_gain2$delta * validation_gain2$short_delay) < validation_gain2$large_amount * exp(-validation_gain2$delta * validation_gain2$long_delay), 1, 
	ifelse(validation_gain2$small_amount * exp(-validation_gain2$delta * validation_gain2$short_delay) > validation_gain2$large_amount * exp(-validation_gain2$delta * validation_gain2$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$expo <- ifelse(exp_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$expo <- aggregate(preds_trials_gain2$expo, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model Mazur
validation_gain2$k_m <- rep(parameters_gain2$k_m, each = num_questions2)   	# create column of each subject's fitted parameter
hyp_m_pred_gain2 <- ifelse(validation_gain2$small_amount / (1 + validation_gain2$k_m * validation_gain2$short_delay) < validation_gain2$large_amount / (1 + validation_gain2$k_m * validation_gain2$long_delay), 1, 
	ifelse(validation_gain2$small_amount / (1 + validation_gain2$k_m * validation_gain2$short_delay) > validation_gain2$large_amount / (1 + validation_gain2$k_m * validation_gain2$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$hyper_m <- ifelse(hyp_m_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$hyper_m <- aggregate(preds_trials_gain2$hyper_m, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model Rachlin
validation_gain2$k_r <- rep(parameters_gain2$k_r, each = num_questions2)   			# create column of each subject's fitted parameter
validation_gain2$sigma <- rep(parameters_gain2$sigma, each = num_questions2)   	# create column of each subject's fitted parameter
hyp_r_pred_gain2 <- ifelse(validation_gain2$small_amount / (1 + validation_gain2$k_r * validation_gain2$short_delay ^ validation_gain2$sigma)  < validation_gain2$large_amount / (1 + validation_gain2$k_r * validation_gain2$long_delay ^ validation_gain2$sigma), 1, ifelse(validation_gain2$small_amount / (1 + validation_gain2$k_r * validation_gain2$short_delay ^ validation_gain2$sigma) > validation_gain2$large_amount / (1 + validation_gain2$k_r * validation_gain2$long_delay ^ validation_gain2$sigma), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$hyper_r <- ifelse(hyp_r_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$hyper_r <- aggregate(preds_trials_gain2$hyper_r, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model Kirby
validation_gain2$k_k <- rep(parameters_gain2$k_k, each = num_questions2)  			# create column of each subject's fitted parameter
validation_gain2$mu <- rep(parameters_gain2$mu, each = num_questions2)  				# create column of each subject's fitted parameter
hyp_k_pred_gain2 <- ifelse(validation_gain2$small_amount / (1 + validation_gain2$k_k * validation_gain2$small_amount ^ validation_gain2$mu * validation_gain2$short_delay) < validation_gain2$large_amount / (1 + validation_gain2$k_k * validation_gain2$large_amount ^ validation_gain2$mu * validation_gain2$long_delay), 1, ifelse(validation_gain2$small_amount / (1 + validation_gain2$k_k * validation_gain2$small_amount ^ validation_gain2$mu * validation_gain2$short_delay) > validation_gain2$large_amount / (1 + validation_gain2$k_k * validation_gain2$large_amount ^ validation_gain2$mu * validation_gain2$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$hyper_k <- ifelse(hyp_k_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$hyper_k <- aggregate(preds_trials_gain2$hyper_k, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model L&P
validation_gain2$alpha <- rep(parameters_gain2$alpha, each = num_questions2)   	# create column of each subject's fitted parameter
validation_gain2$beta <- rep(parameters_gain2$beta, each = num_questions2)   		# create column of each subject's fitted parameter
hyp_lp_pred_gain2 <- ifelse(validation_gain2$small_amount / (1 + validation_gain2$alpha * validation_gain2$short_delay) ^ (validation_gain2$beta / validation_gain2$alpha)  < validation_gain2$large_amount / (1 + validation_gain2$alpha * validation_gain2$long_delay) ^ (validation_gain2$beta / validation_gain2$alpha), 1, ifelse(validation_gain2$small_amount / (1 + validation_gain2$alpha * validation_gain2$short_delay) ^ (validation_gain2$beta / validation_gain2$alpha) > validation_gain2$large_amount / (1 + validation_gain2$alpha * validation_gain2$long_delay) ^ (validation_gain2$beta / validation_gain2$alpha), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$hyper_lp <- ifelse(hyp_lp_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$hyper_lp <- aggregate(preds_trials_gain2$hyper_lp, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Arithmetic model
validation_gain2$lambda <- rep(parameters_gain2$lambda, each = num_questions2)   	# create column of each subject's fitted parameter
arith_pred_gain2 <- ifelse(validation_gain2$small_amount - validation_gain2$lambda * validation_gain2$short_delay < validation_gain2$large_amount - validation_gain2$lambda * validation_gain2$long_delay, 1, ifelse(validation_gain2$small_amount - validation_gain2$lambda - validation_gain2$short_delay > validation_gain2$large_amount - validation_gain2$lambda - validation_gain2$long_delay, 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$arith <- ifelse(arith_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$arith <- aggregate(preds_trials_gain2$arith, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (exponential)--first look for similarity differences, then use exponential discounting if similarity doesn't distinguish
sim_e_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, exp_pred_gain2)	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$sim_e <- ifelse(sim_e_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$sim_e <- aggregate(preds_trials_gain2$sim_e, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--Mazur)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hm_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, hyp_m_pred_gain2)	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$sim_h_m <- ifelse(sim_hm_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$sim_h_m <- aggregate(preds_trials_gain2$sim_h_m, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--Rachlin)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hr_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, hyp_r_pred_gain2)	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$sim_h_r <- ifelse(sim_hr_pred_gain2 == validation_gain2$choice, 1, 0)
preds_subj_gain2$sim_h_r <- aggregate(preds_trials_gain2$sim_h_r, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--Kirby)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hk_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, hyp_k_pred_gain2)  # determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$sim_h_k <- ifelse(sim_hk_pred_gain2 == validation_gain2$choice, 1, 0)
preds_subj_gain2$sim_h_k <- aggregate(preds_trials_gain2$sim_h_k, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--L&P)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hlp_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, hyp_lp_pred_gain2)	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$sim_h_lp <- ifelse(sim_hlp_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$sim_h_lp <- aggregate(preds_trials_gain2$sim_h_lp, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (Arithmetic)--first look for similarity differences, then use Arithmetic discounting if similarity doesn't distinguish
sim_a_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, arith_pred_gain2)	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_gain2$sim_arith <- ifelse(sim_a_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_gain2$sim_arith <- aggregate(preds_trials_gain2$sim_arith, by = list(preds_trials_gain2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity (Leland) model
sim_l_pred_gain2 <- ifelse(validation_gain2$sim_diff == 2 | validation_gain2$sim_diff == 3, validation_gain2$sim_pred, NA)	# create array differentiating between higher sim for amount, higher sim for delay, and equal sim

# Add similarity (Leland) model and other columns to preds_trials_gain2
preds_trials_gain2$sim_l <- ifelse(sim_l_pred_gain2 == validation_gain2$choice, 1, 0)	# create column of correct predictions for similarity
preds_trials_gain2$sim_diff <- validation_gain2$sim_diff															# append similarity difference values
preds_trials_gain2$sim_domain <- ifelse(is.na(validation_gain2$sim_delay), NA, ifelse(is.na(validation_gain2$sim_amt), NA, ifelse(validation_gain2$sim_diff == 1 | validation_gain2$sim_diff == 4, 0, 1)))	# create column with domain, non-domain, and unrated questions
preds_trials_gain2$sim_preds <- sim_l_pred_gain2																			# append similarity predictions
preds_trials_gain2$choice <- validation_gain2$choice																	# append observed choices

# Select questions within similarity domain
subjects_df_gain2 <- data.frame(subject = good_subjects_gain2)	# create data frame of good subjects
domain_gain2 <- subset(preds_trials_gain2, sim_domain == 1)			# select domain questions
domain_sim_gain2 <- aggregate(domain_gain2$sim_l, by = list(domain_gain2$subject), FUN = c("mean", "sum", "length"), na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(domain_sim_gain2) <- c("subject", "acc", "numcorr", "numquestions")
domain_sim_gain2 <-merge(subjects_df_gain2, domain_sim_gain2, all.x = TRUE) 			# merge with subject list to include subjects with no choices in similarity domain
domain_models_subj_gain2 <- aggregate(domain_gain2[, c(2:7, 14)], by = list(domain_gain2$subject), FUN = c("mean"), na.rm = TRUE)		# aggregate similarity predictive accuracy by subject
names(domain_models_subj_gain2)[1] <- "subject"			
domain_models_gain2 <- colMeans(domain_models_subj_gain2[, -1], na.rm = TRUE)			# calculate mean predictive accuracy for models in domain

# Select questions outside of similarity domain
bothsim_gain2 <- subset(preds_trials_gain2, sim_diff == 1)  											# select both similar questions
bothsim_sim_gain2 <- aggregate(bothsim_gain2$choice, by = list(bothsim_gain2$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothsim_sim_gain2) <- c("subject", "acc", "numLL", "numquestions")
bothsim_sim_gain2 <- merge(subjects_df_gain2, bothsim_sim_gain2, all.x = TRUE)		# merge subjects list with subjects with both similar data
bothsim_sim_gain2$obs <- bothsim_sim_gain2$numLL / bothsim_sim_gain2$numquestions # calculate observed proportion choice for LL
bothsim_sim_gain2$deviation <- abs(bothsim_sim_gain2$obs - 0.5)										# calculate deviation from chance (0.5)
bothsim_sim_gain2$acc <- bothsim_sim_gain2$deviation / 0.5												# calculate predictive accuracy for non-domain questions
bothdis_gain2 <- subset(preds_trials_gain2, sim_diff == 4) 												# select both dissimilar questions
bothdis_sim_gain2 <- aggregate(bothdis_gain2$choice, by = list(bothdis_gain2$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothdis_sim_gain2) <- c("subject", "acc", "numLL", "numquestions")
bothdis_sim_gain2$obs <- bothdis_sim_gain2$numLL / bothdis_sim_gain2$numquestions # calculate observed proportion choice for LL
bothdis_sim_gain2$deviation <- abs(bothdis_sim_gain2$obs - 0.5)										# calculate deviation from chance (0.5)
bothdis_sim_gain2$acc <- bothdis_sim_gain2$deviation / 0.5												# calculate predictive accuracy for non-domain questions

# Create overall predictive accuracy for similarity
preds_subj_gain2$sim_l <- rep(NA, num_good_subjects_gain2)	# create vector of NAs for similarity (Leland) predictive accuracy
for(i in 1:num_good_subjects_gain2) {		# calculate predictive accuracy as weighted average of accuracy within similarity domain and deviation from chance outside of domain
    preds_subj_gain2$sim_l[i] <- (domain_sim_gain2$numquestions[i] * (domain_sim_gain2$numcorr[i] / domain_sim_gain2$numquestions[i]) + bothsim_sim_gain2$numquestions[i] * (1 - (abs((bothsim_sim_gain2$numLL[i] / bothsim_sim_gain2$numquestions[i]) - 0.5) / 0.5)) + bothdis_sim_gain2$numquestions[i] * (1 - (abs((bothdis_sim_gain2$numLL[i] / bothdis_sim_gain2$numquestions[i]) - 0.5) / 0.5))) / (domain_sim_gain2$numquestions[i] + bothsim_sim_gain2$numquestions[i] + bothdis_sim_gain2$numquestions[i])
}

parameters_medians_gain2 <- apply(parameters_gain2[, c(2:10)], 2, FUN = median)	# calculate mean parameter values

## Aggregate analysis
accuracy_all_gain2 <- colMeans(preds_subj_gain2[, 2:(total_num_models2 + 1)], na.rm = TRUE) * 100   # calculate mean predictive accuracies for all models in all questions
accuracy_all_gain2 <- accuracy_all_gain2[c(1:6, 13, 7:12)]		# reorder columns to include similarity (Leland) data before other similarity models

## Calculate CIs and effect sizes
pred_acc_subj_gain2 <- stack(preds_subj_gain2[c(2:14)])    		# reshape data to long/stacked format
names(pred_acc_subj_gain2) <- c("pred_acc", "model")
pred_acc_subj_gain2$model <- factor(pred_acc_subj_gain2$model, levels = c("expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_l", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith"))	# reorder model levels
model_names_pred <- c("Exponential", "Hyperbolic (Mazur)", "Hyperbolic (Rachlin)", "Hyperbolic (Kirby)", "Hyperbolic (L&P)", "Arithmetic", "Similarity (Leland)", "Similarity+exponential", "Similarity+Mazur", "Similarity+Rachlin", "Similarity+Kirby", "Similarity+L&P", "Similarity+arithmetic")
pred_acc_subj_gain2$model <- factor(pred_acc_subj_gain2$model, labels = model_names_pred)	# rename model labels
pred_acc_subj_gain2$subject <- rep(unique(preds_subj_gain2$subject), total_num_models2)		# create column of subject numbers

pred_acc_means_gain2 <- aggregate(pred_acc_subj_gain2$pred_acc, by=list(pred_acc_subj_gain2$model), FUN = "mean", na.rm = TRUE)	# aggregate predictive accuracy by subject
names(pred_acc_means_gain2) <- c("model", "pred_acc")
pred_acc_mean_gain2 <- normDataWithin(data = pred_acc_subj_gain2, idvar = "subject", measurevar="pred_acc", na.rm = TRUE)	# calculate normalized data for within-subjects CIs
pred_acc_norm_gain2 <- summarySEwithin(data = pred_acc_mean_gain2, measurevar="pred_accNormed", withinvars = "model", idvar = "subject", na.rm = TRUE)	# calculate within-subjects CIs
pred_acc_norm_gain2$pred_acc <- pred_acc_means_gain2$pred_acc											# copy over actual means
pred_acc_norm_gain2$lci <- pred_acc_norm_gain2$pred_acc - pred_acc_norm_gain2$ci  # add CIs
pred_acc_norm_gain2$uci <- pred_acc_norm_gain2$pred_acc + pred_acc_norm_gain2$ci	# add CIs

## Test difference between Mazur and similarity+Mazur
hyper_sim_diff_subj_gain2 <- preds_subj_gain2$sim_h_m - preds_subj_gain2$hyper_m										# calculate difference between predictive accuracy of Mazur and similarity+Mazur
hyper_sim_diff_ci_gain2 <- ci(hyper_sim_diff_subj_gain2)																						# calculate mean, sd, CIs of difference
hyper_sim_diff_ci_gain2$ci <- hyper_sim_diff_ci_gain2$mean - hyper_sim_diff_ci_gain2$lower95ci			# find CI difference
hyper_sim_diff_ci_gain2$effect_size_d <-  hyper_sim_diff_ci_gain2$mean / hyper_sim_diff_ci_gain2$sd	# calculate Cohen's d effect size

#####
## Calculate predictive accuracy for loss data
#####

## Prepare data
validation_loss2 <- subset(validation2, gainloss == "loss" & subject %in% good_subjects_loss2) 	# extract loss data for subjects with complete training data
validation_loss2$small_amount <- -validation_loss2$small_amount																	# make small amounts negative
validation_loss2$large_amount <- -validation_loss2$large_amount																	# make large amounts negative

## Create data frame of questions and mean responses for loss data
validation_loss2_aggregation <- aggregate(validation_loss2$choice, by = list(validation_loss2$small_amount, validation_loss2$large_amount, validation_loss2$short_delay, validation_loss2$long_delay, validation_loss2$question), FUN = c("mean", "sd")) # calculate mean and sd for each question
validation_loss2_aggregation <- validation_loss2_aggregation[, -5]    							# remove question column
validation_loss2_aggregation <- validation_loss2_aggregation[, c(3, 4, 1, 2, 5, 6)] # reorder columns with delays first
names(validation_loss2_aggregation) <- c("Short delay", "Long delay", "Small amount", "Large amount", "Mean choice for LL", "Standard deviation")
validation_loss2_aggregation <- validation_loss2_aggregation[order(validation_loss2_aggregation$"Short delay", validation_loss2_aggregation$"Long delay", validation_loss2_aggregation$"Small amount", validation_loss2_aggregation$"Large amount"), ]																			# sort columns by delays and amounts

# Similarity data in similarity domain
validation_loss2$sim_pred <- ifelse(validation_loss2$sim_diff == 2, 1, ifelse(validation_loss2$sim_diff == 3, 0, NA)) 	# generate predictions when similarity can make a prediction (only one of two attributes is similar)
validation_loss2$sim_acc <- ifelse(!is.na(validation_loss2$sim_pred), ifelse(validation_loss2$choice == validation_loss2$sim_pred, 1, 0), NA) # assess whether prediction matches observed choice
validation_loss2 <- validation_loss2[order(validation_loss2$subject, validation_loss2$question),]												# reorder by subject and question
sim_acc_loss2 <- mean(validation_loss2$sim_acc, na.rm = TRUE) 																													# calculate overall mean predictive accuracy of similarity when it can make predictions
sim_acc_subject_loss2 <- aggregate(sim_acc ~ subject, data = validation_loss2, FUN = mean, na.rm = TRUE) 								# calculate subject-wise mean predictive accuracy of similarity when it can make predictions

## Predict binary choice data (generalization/validation)
preds_trials_loss2 <- data.frame(validation_loss2$subject, matrix(rep(NA, length(validation_loss2$subject) * (total_num_models2 - 1)), nrow = length(validation_loss2$subject)))  # create data frame of NAs for trial predictions
names(preds_trials_loss2) <- c("subject", "expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith")
preds_subj_loss2 <- data.frame(good_subjects_loss2, matrix(rep(NA, num_good_subjects_loss2 * (total_num_models2 - 1)), nrow = num_good_subjects_loss2))	# create data frame of NAs for subject predictions
names(preds_subj_loss2) <- names(preds_trials_loss2) 				# create data frame of NAs for subject predictions
num_questions2 <- length(unique(validation_loss2$question))	# find number of questions

# Exponential discounting model
validation_loss2$delta <- rep(parameters_loss2$delta, each = num_questions2)   	# create column of 46 consecutives instances of each subject's fitted parameter (46 questions have simliarity judgments)
exp_pred_loss2 <- ifelse(validation_loss2$small_amount * exp(-validation_loss2$delta * validation_loss2$short_delay) < validation_loss2$large_amount * exp(-validation_loss2$delta * validation_loss2$long_delay), 1, 
  ifelse(validation_loss2$small_amount * exp(-validation_loss2$delta * validation_loss2$short_delay) > validation_loss2$large_amount * exp(-validation_loss2$delta * validation_loss2$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_loss2$expo <- ifelse(exp_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$expo <- aggregate(preds_trials_loss2$expo, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model Mazur
validation_loss2$k_m <- rep(parameters_loss2$k_m, each = num_questions2)
hyp_m_pred_loss2 <- ifelse(validation_loss2$small_amount / (1 + validation_loss2$k_m * validation_loss2$short_delay) < validation_loss2$large_amount / (1 + validation_loss2$k_m * validation_loss2$long_delay), 1, ifelse(validation_loss2$small_amount / (1 + validation_loss2$k_m * validation_loss2$short_delay) > validation_loss2$large_amount / (1 + validation_loss2$k_m * validation_loss2$long_delay), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_loss2$hyper_m <- ifelse(hyp_m_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$hyper_m <- aggregate(preds_trials_loss2$hyper_m, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model Rachlin
validation_loss2$k_r <- rep(parameters_loss2$k_r, each = num_questions2)
validation_loss2$sigma <- rep(parameters_loss2$sigma, each = num_questions2)
hyp_r_pred_loss2 <- ifelse(validation_loss2$small_amount / (1 + validation_loss2$k_r * validation_loss2$short_delay ^ validation_loss2$sigma)  < validation_loss2$large_amount / (1 + validation_loss2$k_r * validation_loss2$long_delay ^ validation_loss2$sigma), 1, ifelse(validation_loss2$small_amount / (1 + validation_loss2$k_r * validation_loss2$short_delay ^ validation_loss2$sigma) > validation_loss2$large_amount / (1 + validation_loss2$k_r * validation_loss2$long_delay ^ validation_loss2$sigma), 0, NA))  # determine whether model predicts choosing LL (1) or SS (0)
preds_trials_loss2$hyper_r <- ifelse(hyp_r_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$hyper_r <- aggregate(preds_trials_loss2$hyper_r, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model Kirby
validation_loss2$k_k <- rep(parameters_loss2$k_k, each = num_questions2)
validation_loss2$mu <- rep(parameters_loss2$mu, each = num_questions2)
hyp_k_pred_loss2 <- ifelse(validation_loss2$small_amount / (1 + validation_loss2$k_k * abs(validation_loss2$small_amount) ^ validation_loss2$mu * validation_loss2$short_delay) < validation_loss2$large_amount / (1 + validation_loss2$k_k * abs(validation_loss2$large_amount) ^ validation_loss2$mu * validation_loss2$long_delay), 1, ifelse(validation_loss2$small_amount / (1 + validation_loss2$k_k * abs(validation_loss2$small_amount) ^ validation_loss2$mu * validation_loss2$short_delay) > validation_loss2$large_amount / (1 + validation_loss2$k_k * abs(validation_loss2$large_amount) ^ validation_loss2$mu * validation_loss2$long_delay), 0, NA))  # determine whether model predicts choosing LL (1) or SS (0)
preds_trials_loss2$hyper_k <- ifelse(hyp_k_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$hyper_k <- aggregate(preds_trials_loss2$hyper_k, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Hyperbolic discounting model L&P
validation_loss2$alpha <- rep(parameters_loss2$alpha, each = num_questions2)
validation_loss2$beta <- rep(parameters_loss2$beta, each = num_questions2)
hyp_lp_pred_loss2 <- ifelse(validation_loss2$small_amount / (1 + validation_loss2$alpha * validation_loss2$short_delay) ^ (validation_loss2$beta / validation_loss2$alpha)  < validation_loss2$large_amount / (1 + validation_loss2$alpha * validation_loss2$long_delay) ^ (validation_loss2$beta / validation_loss2$alpha), 1, ifelse(validation_loss2$small_amount / (1 + validation_loss2$alpha * validation_loss2$short_delay) ^ (validation_loss2$beta / validation_loss2$alpha) > validation_loss2$large_amount / (1 + validation_loss2$alpha * validation_loss2$long_delay) ^ (validation_loss2$beta / validation_loss2$alpha), 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_loss2$hyper_lp <- ifelse(hyp_lp_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$hyper_lp <- aggregate(preds_trials_loss2$hyper_lp, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Arithmetic model
validation_loss2$lambda <- rep(parameters_loss2$lambda, each = num_questions2)
arith_pred_loss2 <- ifelse(validation_loss2$small_amount - validation_loss2$lambda * validation_loss2$short_delay < validation_loss2$large_amount - validation_loss2$lambda * validation_loss2$long_delay, 1,   
  ifelse(validation_loss2$small_amount - validation_loss2$lambda - validation_loss2$short_delay > validation_loss2$large_amount - validation_loss2$lambda - validation_loss2$long_delay, 0, NA))	# determine whether model predicts choosing LL (1) or SS (0)
preds_trials_loss2$arith <- ifelse(arith_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$arith <- aggregate(preds_trials_loss2$arith, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (exponential)--first look for similarity differences, then use exponential discounting if similarity doesn't distinguish
sim_e_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, exp_pred_loss2)
preds_trials_loss2$sim_e <- ifelse(sim_e_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$sim_e <- aggregate(preds_trials_loss2$sim_e, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--Mazur)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hm_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, hyp_m_pred_loss2)
preds_trials_loss2$sim_h_m <- ifelse(sim_hm_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$sim_h_m <- aggregate(preds_trials_loss2$sim_h_m, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--Rachlin)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hr_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, hyp_r_pred_loss2)
preds_trials_loss2$sim_h_r <- ifelse(sim_hr_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$sim_h_r <- aggregate(preds_trials_loss2$sim_h_r, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--Kirby)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hk_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, hyp_k_pred_loss2)
preds_trials_loss2$sim_h_k <- ifelse(sim_hk_pred_loss2 == validation_loss2$choice, 1, 0)  # create column of predictive accuracy
preds_subj_loss2$sim_h_k <- aggregate(preds_trials_loss2$sim_h_k, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (hyperbolic--L&P)--first look for similarity differences, then use hyperbolic discounting if similarity doesn't distinguish
sim_hlp_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, hyp_lp_pred_loss2)
preds_trials_loss2$sim_h_lp <- ifelse(sim_hlp_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$sim_h_lp <- aggregate(preds_trials_loss2$sim_h_lp, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (Arithmetic)--first look for similarity differences, then use Arithmetic discounting if similarity doesn't distinguish
sim_a_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, arith_pred_loss2)
preds_trials_loss2$sim_arith <- ifelse(sim_a_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of predictive accuracy
preds_subj_loss2$sim_arith <- aggregate(preds_trials_loss2$sim_arith, by = list(preds_trials_loss2$subject), FUN = mean)[, 2]	# calculate mean predictive accuracy for each subject

# Similarity model (Leland)
sim_l_pred_loss2 <- ifelse(validation_loss2$sim_diff == 2 | validation_loss2$sim_diff == 3, validation_loss2$sim_pred, NA)	# create array differentiating between higher sim for amount, higher sim for delay, and equal sim

# Add Leland similarity model and other columns to preds_trials
preds_trials_loss2$sim_l <- ifelse(sim_l_pred_loss2 == validation_loss2$choice, 1, 0)	# create column of correct predictions for similarity
preds_trials_loss2$sim_diff <- validation_loss2$sim_diff															# append similarity difference values
preds_trials_loss2$sim_domain <- ifelse(is.na(validation_loss2$sim_delay), NA, ifelse(is.na(validation_loss2$sim_amt), NA, ifelse(validation_loss2$sim_diff == 1 | validation_loss2$sim_diff == 4, 0, 1)))	# create column with domain, non-domain, and unrated questions
preds_trials_loss2$sim_preds <- sim_l_pred_loss2																			# append similarity predictions
preds_trials_loss2$choice <- validation_loss2$choice																	# append observed choices

# Select questions within similarity domain
subjects_df_loss2 <- data.frame(subject = good_subjects_loss2)								# create data frame of good subjects
domain_loss2 <- subset(preds_trials_loss2, sim_domain == 1)										# select domain questions
domain_sim_loss2 <- aggregate(domain_loss2$sim_l, by = list(domain_loss2$subject), FUN = c("mean", "sum", "length"), na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(domain_sim_loss2) <- c("subject", "acc", "numcorr", "numquestions")
domain_sim_loss2 <-merge(subjects_df_loss2, domain_sim_loss2, all.x = TRUE) 	# merge with subject list to include subjects with no choices in similarity domain
domain_models_subj_loss2 <- aggregate(domain_loss2[, c(2:6, 12)], by = list(domain_loss2$subject), FUN = c("mean"), na.rm = TRUE)		# aggregate similarity predictive accuracy by subject
names(domain_models_subj_loss2)[1] <- "subject"			
domain_models_loss2 <- colMeans(domain_models_subj_loss2[, -1], na.rm = TRUE)	# calculate mean predictive accuracy for models in domain

# Select questions outside of similarity domain
bothsim_loss2 <- subset(preds_trials_loss2, sim_diff == 1)  											# select both similar questions
bothsim_sim_loss2 <- aggregate(bothsim_loss2$choice, by = list(bothsim_loss2$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothsim_sim_loss2) <- c("subject", "acc", "numLL", "numquestions")
bothsim_sim_loss2 <- merge(subjects_df_gain2, bothsim_sim_loss2, all.x = TRUE)		# merge subjects list with subjects with both
bothsim_sim_loss2$obs <- bothsim_sim_loss2$numLL / bothsim_sim_loss2$numquestions # calculate observed proportion choice for LL
bothsim_sim_loss2$deviation <- abs(bothsim_sim_loss2$obs - 0.5)										# calculate deviation from chance (0.5)
bothsim_sim_loss2$acc <- bothsim_sim_loss2$deviation / 0.5												# calculate predictive accuracy for non-domain questions
bothdis_loss2 <- subset(preds_trials_loss2, sim_diff == 4) 												# select both dissimilar questions
bothdis_sim_loss2 <- aggregate(bothdis_loss2$choice, by = list(bothdis_loss2$subject), FUN = c("mean", "sum", "length"),na.rm = TRUE)	# aggregate similarity predictive accuracy by subject
names(bothdis_sim_loss2) <- c("subject", "acc", "numLL", "numquestions")
bothdis_sim_loss2$obs <- bothdis_sim_loss2$numLL / bothdis_sim_loss2$numquestions # calculate observed proportion choice for LL
bothdis_sim_loss2$deviation <- abs(bothdis_sim_loss2$obs - 0.5)										# calculate deviation from chance (0.5)
bothdis_sim_loss2$acc <- bothdis_sim_loss2$deviation / 0.5												# calculate predictive accuracy for non-domain questions

# Create overall predictive accuracy for similarity
preds_subj_loss2$sim_l <- rep(NA, num_good_subjects_loss2)	# create vector of NAs for similarity (Leland) predictive accuracy
for(i in 1:num_good_subjects_loss2) {		# calculate predictive accuracy as weighted average of accuracy within similarity domain and deviation from chance outside of domain
  preds_subj_loss2$sim_l[i] <- (domain_sim_loss2$numquestions[i] * (domain_sim_loss2$numcorr[i] / domain_sim_loss2$numquestions[i]) + bothsim_sim_loss2$numquestions[i] * (1 - (abs((bothsim_sim_loss2$numLL[i] / bothsim_sim_loss2$numquestions[i]) - 0.5) / 0.5)) + bothdis_sim_loss2$numquestions[i] * (1 - (abs((bothdis_sim_loss2$numLL[i] / bothdis_sim_loss2$numquestions[i]) - 0.5) / 0.5))) / (domain_sim_loss2$numquestions[i] + bothsim_sim_loss2$numquestions[i] + bothdis_sim_loss2$numquestions[i])
}

parameters_medians_loss2 <- apply(parameters_loss2[, c(2:10)], 2, FUN = median)	# calculate mean parameter values

## Aggregate analysis
accuracy_all_loss2 <- colMeans(preds_subj_loss2[, 2:(total_num_models2 + 1)], na.rm = TRUE) * 100   # calculate mean predictive accuracies for all models in all questions
accuracy_all_loss2 <- accuracy_all_loss2[c(1:6, 13, 7:12)]

## Calculate CIs and effect sizes
pred_acc_subj_loss2 <- stack(preds_subj_loss2[c(2:14)])      		# reshape data to long/stacked format
names(pred_acc_subj_loss2) <- c("pred_acc", "model")
pred_acc_subj_loss2$model <- factor(pred_acc_subj_loss2$model, levels = c("expo", "hyper_m", "hyper_r", "hyper_k", "hyper_lp", "arith", "sim_l", "sim_e", "sim_h_m", "sim_h_r", "sim_h_k", "sim_h_lp", "sim_arith"))	# reorder model levels
model_names_pred <- c("Exponential", "Hyperbolic (Mazur)", "Hyperbolic (Rachlin)", "Hyperbolic (Kirby)", "Hyperbolic (L&P)", "Arithmetic", "Similarity (Leland)", "Similarity+exponential", "Similarity+Mazur", "Similarity+Rachlin", "Similarity+Kirby", "Similarity+L&P", "Similarity+arithmetic")
pred_acc_subj_loss2$model <- factor(pred_acc_subj_loss2$model, labels = model_names_pred)	# rename model labels
pred_acc_subj_loss2$subject <- rep(unique(preds_subj_loss2$subject), total_num_models2)		# create column of subject numbers

pred_acc_means_loss2 <- aggregate(pred_acc_subj_loss2$pred_acc, by=list(pred_acc_subj_loss2$model), FUN = "mean", na.rm = TRUE)	# aggregate predictive accuracy by subject
names(pred_acc_means_loss2) <- c("model", "pred_acc")
pred_acc_mean_loss2 <- normDataWithin(data = pred_acc_subj_loss2, idvar = "subject", measurevar="pred_acc", na.rm = TRUE)	# calculate normalized data for within-subjects CIs
pred_acc_norm_loss2 <- summarySEwithin(data = pred_acc_mean_loss2, measurevar="pred_accNormed", withinvars = "model", idvar = "subject", na.rm = TRUE)	# calculate within-subjects CIs
pred_acc_norm_loss2$pred_acc <- pred_acc_means_loss2$pred_acc											# copy over actual means
pred_acc_norm_loss2$lci <- pred_acc_norm_loss2$pred_acc - pred_acc_norm_loss2$ci  # add CIs
pred_acc_norm_loss2$uci <- pred_acc_norm_loss2$pred_acc + pred_acc_norm_loss2$ci	# add CIs
# pred_acc_norm_loss2 <- rbind(pred_acc_norm_loss2[1:3, ], rep(NA, 9), pred_acc_norm_loss2[4:9, ], rep(NA, 9), pred_acc_norm_loss2[10:11, ])

## Test difference between Mazur and similarity+Mazur
hyper_sim_diff_subj_loss2 <- preds_subj_loss2$sim_h_m - preds_subj_loss2$hyper_m										# calculate difference between predictive accuracy of Mazur and similarity+Mazur
hyper_sim_diff_ci_loss2 <- ci(hyper_sim_diff_subj_loss2)																						# calculate mean, sd, CIs of difference
hyper_sim_diff_ci_loss2$ci <- hyper_sim_diff_ci_loss2$mean - hyper_sim_diff_ci_loss2$lower95ci			# find CI difference
hyper_sim_diff_ci_loss2$effect_size_d <-  hyper_sim_diff_ci_loss2$mean / hyper_sim_diff_ci_loss2$sd	# calculate Cohen's d effect size

#####
## Compare gain and loss data
#####

## Create xtable of questions and mean responses for both gain and loss data
validation2_aggregation <- data.frame(validation_gain2_aggregation, validation_loss2_aggregation[, 5:6]) # merge gain and loss data
names(validation2_aggregation) <- c("Short delay", "Long delay", "Small amount", "Large amount", "Choice LL (gain)", "SD (gain)", "Choice LL (loss)", "SD (loss)")
validation2_xtable <- xtable(validation2_aggregation, caption = 'Questions and Mean Responses for Experiment 3', label = "tab:EXP3QUEST")	# create xtable for import into LaTeX

## Compare discount rates for gains and losses
# Calulate MLE for gain data using the good loss subjects (to compare directly with loss data)
train_gain_loss2 <- subset(validation2, subject %in% good_subjects_loss2 & gainloss == "gain")	# select gain data for good subjects (in loss condition)
expon_fit_group_gainloss2 <- mle2(choice ~ dbinom(prob = invlogit(expon(large_amount, long_delay, delta) - expon(small_amount, short_delay, delta)), size = 1), data = train_gain_loss2, start = list(delta = 0.05))  # calculate MLE for exponential discounting model
expon_par_group_gainloss2 <- coef(expon_fit_group_gainloss2)			# extract delta
expon_ci_group_gainloss2 <- confint(expon_fit_group_gainloss2)		# calculate delta CIs
expon_ci_diff_group_gainloss2 <- expon_par_group_gainloss2 - expon_ci_group_gainloss2[1]  	# calculate delta CI difference

## Compare similarity for gains and losses
sim_subj <- aggregate(similarity ~ gainloss * subject, sim_amt2, mean)	# calculate mean similarity jugdments over gain/loss condition and subjects
sim_gain <- subset(sim_subj, gainloss == "gain")$similarity							# select similarity judgements from gain data
sim_gain_mean <- ci(sim_gain)																						# calculate mean, sd, and CIs
sim_loss <- subset(sim_subj, gainloss == "loss")$similarity							# select similarity judgements from loss data
sim_loss_mean <- ci(sim_loss)																						# calculate mean, sd, and CIs
sim_gainloss_diff <- sim_gain - sim_loss																# calculate similarity judgment differences across gain/loss conditions
sim_gainloss_diff_mean <- ci(sim_gainloss_diff)													# calculate mean, sd, and CIs
sim_gainloss_diff_mean$ci <- sim_gainloss_diff_mean$mean - sim_gainloss_diff_mean$lower95ci				# calculate CI differences
sim_gainloss_diff_mean$effect_size_d <- sim_gainloss_diff_mean$mean / sim_gainloss_diff_mean$sd		# calculate Cohen's d effect size

#####
## Compare predictive accuracy in similarity vs. discounting models
#####

sim_model_diff1 <- pred_acc_norm1$pred_acc[8:13] - pred_acc_norm1$pred_acc[1:6]									# calculate difference in predictive accuracy between similarity and discounting models	for Expt 1
sim_model_diff_gain2 <- pred_acc_norm_gain2$pred_acc[8:13] - pred_acc_norm_gain2$pred_acc[1:6]	# calculate difference in predictive accuracy between similarity and discounting models	for Expt 2 gains condition

