###################################################
### stevens.2007.animcog.rcode.R
### Created by Jeffrey R. Stevens on 17 Aug 2010
### Summary: This script generates figures and calculates inferential statistics 
###  for the analysis of number discrimination in tamarins:
###  Stevens, J.R., Wood, J.N., & Hauser, M.D. 2007. Quantity trumps number: 
###   discrimination experiments in cotton-top tamarins (Saguinus oedipus) 
###   and common marmosets (Callithrix jacchus).  Animal Cognition 10:429-437.
### Instructions: Place this file and the data files (stevens.2007.animcog.data*.txt,
###  where * represents 1-3; please find description of data files below)
###  in the same directory. In R set the working directory to this directory. Type 
###  > source("stevens.2007.animcog.rcode.R")
###  This will run the script, adding all of the calculated variables to the
###  workspace and saving PDF versions of the figures in the figures directory.
###  **Note that this script was not used to generate the graphs or analyses in the
###    original paper.  It produces box plots instead of bar plots.**
### Uses: This script can be reproduced and modified for personal and scientific use.
###################################################


#################
## Experiment 1
#################

## Load data
exp1 <- read.csv("stevens.2007.animcog.data1.txt")					# load data for experiment 1
tam <- subset(exp1, species == "T")							# extract tamarin data
marm <- subset(exp1, species == "M")							# extract marmoset data

## Aggregate data by ratio
both <- aggregate(exp1$correct, by = list(exp1$ratio), FUN = mean, na.rm = T)		# aggregate correct numerical discrimination by numerical ratio for both species
tam.m <- aggregate(tam$correct, by = list(tam$ratio), FUN = mean, na.rm = T)		# aggregate correct numerical discrimination by numerical ratio for tamarins
marm.m <- aggregate(marm$correct, by = list(marm$ratio), FUN = mean, na.rm = T)		# aggregate correct numerical discrimination by numerical ratio for marmosets
ratio.means <- cbind(both, tam.m$x, marm.m$x)					# combine into a single data frame
names(ratio.means) <- c("ratio", "both", "tamarins", "marmosets")			# rename columns of data frame

## Generate Figure 2
pdf(file = "figure2.pdf", width = 6, height = 5)
exp1.plot <- xyplot(both + tamarins + marmosets ~ ratio, data = ratio.means, type = "l",
  ylim = c(0.48, 1.02), ylab = "Proportion choosing larger reward", xlab = "Ratio of larger to smaller reward",
  lwd = c(4, 2, 2), lty = c(1, 2, 3), col = c("black", "red", "blue"), cex = 1.25, aspect = "xy",
    key = list(corner = c(0.05, 0.95), padding.text = 3, cex = 1.25,
      text = list(c("Combined", "Tamarins", "Marmosets"), adj = 1), 
      lines = list(col = c("black", "red", "blue")), lwd = c(4, 2, 2), lty = c(1, 2, 3), pch = 1, type = c("l", "l", "l"), divide = 1)
)
plot(exp1.plot)
dev.off()

## Generate ANOVA for ratio effect on transformed proportion correct
exp1.subj <- aggregate(exp1$correct, by = list(exp1$ratio, exp1$species, exp1$subject),
  FUN = mean, na.rm = T)								# aggregate correct numerical discrimination by ratio and subject
names(exp1.subj) <- c("ratio", "species", "subject", "correct")
exp1.subj$tcorr <- asin(sqrt(exp1.subj$correct))					# create arcsine, square-root transformed proportion correct
exp1.aov <- summary(aov(tcorr ~ as.factor(ratio) * species + Error(subject), data = exp1.subj))	# generate repeated-measures ANOVA of ratio effect

## Generate ANOVA for trial number effect on transformed proportion correct
exp1.trial <- aggregate(exp1$correct, by = list(exp1$trial., exp1$species, exp1$subject),
  FUN = mean, na.rm = T)								# aggregate correct numerical discrimination by trial number and subject
names(exp1.trial) <- c("trial", "species", "subject", "correct")
exp1.trial$tcorr <- asin(sqrt(exp1.trial$correct))					# create arcsine, square-root transformed proportion correct
exp1.trial.aov <- summary(aov(tcorr ~ as.factor(trial) * species + Error(subject), data = exp1.trial))	# generate repeated-measures ANOVA of trial number

## Generate ANOVA for session effect on transformed proportion correct
exp1.session <- aggregate(exp1$correct, by = list(exp1$session, exp1$species, exp1$subject),
  FUN = mean, na.rm = T)								# aggregate correct numerical discrimination by session number and subject
names(exp1.session) <- c("session", "species", "subject", "correct")
exp1.session$tcorr <- asin(sqrt(exp1.session$correct))					# create arcsine, square-root transformed proportion correct
exp1.session.aov <- summary(aov(tcorr ~ as.factor(session) * species + Error(subject), data = exp1.session))	# generate repeated-measures ANOVA of session number

#################
## Experiment 2
#################

## Load data
exp2 <- read.csv("stevens.2007.animcog.data2.txt")

## Aggregate data by correlation condition and subject
exp2.subj <- aggregate(exp2$correct, by = list(exp2$sizecorrelation, exp2$subject),
  FUN = mean, na.rm = T)								# aggregate correct numerical discrimination by correlation condition and subject
names(exp2.subj) <- c("correlation", "subject", "correct")
exp2.subj$correlation <- relevel(exp2.subj$correlation, ref = "=")			# reorder levels to match Figure 4
exp2.subj$correlation <- relevel(exp2.subj$correlation, ref = "+")			# reorder levels to match Figure 4
levels(exp2.subj$correlation) <- c("Greater amount", "Equal amount", "Smaller amount")	# relabel levels to match Figure 4

## Generate Figure 4
pdf(file = "figure4.pdf", width = 6, height = 5)
exp2.plot <- bwplot(correct ~ correlation, data = exp2.subj, ylim=c(-0.05, 1.05),
  ylab = "Proportion choosing larger reward", xlab = "Correlation condition",
  panel = function(x, y, groups) {
    panel.bwplot(x, y, pch = "|", horizontal = F, coef = 0, col = "black", box.ratio = 2)	# plots boxplot with medians as bars
    mean.values <<- tapply(y, x, mean)								# calculates means
    panel.points(mean.values, pch = 18, cex = 1.5, col = "black")				# plots means as diamonds
  }
)
plot(exp2.plot)
dev.off()

## Generate ANOVA for correlation condition effect on proportion correct
exp2.subj2 <- aggregate(exp2$correct, by = list(exp2$sizecorrelation, exp2$numcond, exp2$density, exp2$subject),
  FUN = mean, na.rm = T)									# aggregate by correlation, numerical condition, density, and subject
names(exp2.subj2) <- c("correlation", "numcond", "density", "subject", "correct")
exp2.subj2$tcorr <- asin(sqrt(exp2.subj2$correct))						# create arcsine, square-root transformed proportion correct
exp2.aov <- summary(aov(tcorr ~ numcond * density * correlation + 
  Error(subject/(numcond * density * correlation)), data = exp2.subj2))				# generate repeated-measures ANOVA of correlation, condition, and density

#################
## Experiment 3
#################

## Load data
exp3 <- read.csv("stevens.2007.animcog.data3.txt")

## Aggregate data by condition and subject
exp3.subj <- aggregate(exp3$correct, by = list(exp3$cond, exp3$subject), FUN = mean)
names(exp3.subj) <- c("cond", "subject", "correct")						# aggregate by condition and subject
exp3.subj$cond <- relevel(exp3.subj$cond, ref = "+L")						# reorder levels to match Figure 5
exp3.subj$cond <- relevel(exp3.subj$cond, ref = "=L")						# reorder levels to match Figure 5
levels(exp3.subj$cond) <- c("Same/Asymmetric", "Different/Asymmetric", "Different/Symmetric")	# relabel levels to match Figure 5

## Generate Figure 5
pdf(file = "figure5.pdf", width = 6, height = 5)
exp3.plot <- bwplot(correct ~ cond, data = exp3.subj, ylim=c(-0.05, 1.05),
  ylab = "Proportion choosing larger reward", xlab = "Condition",
  panel = function(x, y, groups) {
    panel.bwplot(x, y, pch = "|", horizontal = F, coef = 0, col = "black", box.ratio = 2)	# plots boxplot with medians as bars
    mean.values <<- tapply(y, x, mean)								# calculates means
    panel.points(mean.values, pch = 18, cex = 1.5, col = "black")				# plots means as diamonds
  }
)
plot(exp3.plot)
dev.off()

## Generate t-test and ANOVA
exp3.ttest <- t.test(subset(exp3.subj$correct, exp3.subj$cond == "Same/Asymmetric"), mu = 0.5)	# generate t-test comparing same/asymmetric condition to 0.5
exp3.subj$tcorr <- asin(sqrt(exp3.subj$correct))						# create arcsine, square-root transformed proportion correct
exp3.aov <- summary(aov(tcorr ~ cond + Error(subject), data = exp3.subj))			# generate repeated-measures ANOVA of condition effect

###################################################
### Description of data files
###################################################
#################
## stevens.2007.animcog.data1.txt
## This file contains data for experiment 1
#################
## subject = name of subject
## species = species (T=tamarin, M=marmoset)
## session = session number
## trial# = trial number
## L = number of pellets on left tray
## R = number of pellets on right tray
## choice = chosen amount of pellets
## side = chosen side
## larger = larger amount of pellets
## smaller = smaller amount of pellets
## ratio = ratio of larger to smaller amount of pellets
## dist = numerical distance between amount of pellets
## mag = smaller amount of pellets
## correct = choice for larger amount of pellets (1 = choosing larger, 0 = choosing smaller)
#################


#################
## stevens.2007.animcog.data2.txt
## This file contains data for experiment 2
#################
## subject = name of subject
## session = session number
## trial# = trial number
## numcond = number condition (1v2 = 1 vs. 2 pellets, 4v8 = 4 vs. 8 pellets)
## sizecorrelation = correlation between number and amount (+ larger number has greater amount, = larger number has equal amount, - larger number has smaller amount)
## sidelarger = side which has larger amount (R = right, L = left)
## density = whether density is the same/symmetric (B) or different/asymmetric (L) between options
## sidechosen = chosen side
## correctchoice = binary description of whether they chose the larger amount (Y = yes, N = no)
## correct = choice for larger amount of pellets (1 = choosing larger, 0 = choosing smaller)
#################

#################
## stevens.2007.animcog.data3.txt
## This file contains data for experiment 3
#################
## subject = name of subject
## session = session number
## trial# = trial number
## size = pellet size (= same, + different)
## density = whether density is the same/symmetric (B) or different/asymmetric (L) between options
## cond = condition, combination of size and density (=L same/asymmetric, +L different/asymmetric, +B different/symmetric)
## side = side which has larger amount (R = right, L = left)
## choice = chosen side
## correct = choice for larger amount of pellets (1 = choosing larger, 0 = choosing smaller)
#################

