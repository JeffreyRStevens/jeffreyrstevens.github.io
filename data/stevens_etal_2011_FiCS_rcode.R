###################################################
### stevens_rcode.R
### Created by Jeffrey R. Stevens on 1 Nov 2010 (jstevens@mpib-berlin.mpg.de)
### Summary: This script calculates descriptive statistics and generates figures
### 	for the analysis of the cooperative memory experiment, agent-based
###	simulation and game theoretical analysis for the paper entitled
###	"Forgetting constrains the emergence of cooperative decision strategies"
###	(Stevens, Volstorf, Schooler, & Rieskamp).
### Instructions: Place this file and the data files (stevens_memory_data.csv and
###	stevens_simulation_data.csv) in the same directory.  Create a folder 
###	called "figures". In R set the working directory to this directory.  Type 
### 	> source("stevens_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PDF versions of the figures in the figures directory.
### Description of the data columns for stevens_memory_data.csv:
### 	Subject - name of subject
###	Age - age of subject
###	correctresponses - total number of correct responses
###	numquests - total number of questions received
###	SessionDate - date of session
###	SessionTime - time of session
###	Sex - gender of participant
###	TotalPayment - total amount of money received (in euros)
###	Trial - the trial number
###	CorrectAnswer - the correct key that the participant should press (p = not cooperate and q = cooperate)
###	GroupNum - the number for the current group (or replicate)
###	GroupSize - the number of partners in the group
###	GroupTrials - the trial number within a group
###	OldState.ACC - binary value for whether participant was correct on this trial (1 = correct, 0 = incorrect)
###	PartnerAction - partner's action chosen in current trial
###	PartnerName - name for partner
###	ImageFile - name of image file used for partner
###	Round - the number for the current round
###	RoundTrials - the trial number within a round
###	SubjectType - the number representing which of 9 conditions the participant received
###	AccumulatedPayment - accumulated amount of money received at that trial
###	TotalGroups - total number of groups experienced
###	TotalRounds - total number of rounds experienced
###	TotalTrials - total number of trials experienced
###	TrialsPerGroup - number of trials experienced per group
###	InterveningEvents - number of intervening events since last interaction with this partner
### Description of data columns for stevens_simulation_data.csv
###	Error - error rate
###	gTFT - proportion of simulations won by GTFT strategists
###	WSLS - proportion of simulations won by WSLS strategists
###	RANDOM - proportion of simulations won by RAND strategists
###	GRIM - proportion of simulations won by GRIM strategists
###	cTFT - proportion of simulations won by CTFT strategists
###	TFT - proportion of simulations won by TFT strategists
###	TF2T - proportion of simulations won by TF2T strategists
###	AllC - proportion of simulations won by ALLC strategists
###	AllD - proportion of simulations won by ALLD strategists
###	Coop - proportion of cooperative interactions in final generation
### Permission: This script can be reproduced and modified for personal and scientific use.
###################################################

###################################################
### Load libraries, data, and R version
###################################################
library(lattice)        						# required for plots
library(epicalc)							# required for aggregatation
library(Hmisc)								# required for plotting error bars
library(foreach)							# required for foreach function
all.data <- read.csv("data/stevens_memory_data.csv")			# cooperative memory experimental data
sim.all <- read.csv("data/stevens_simulation_data.csv")			# agent-based simulation data
ver <- getRversion()							# get R version

###################################################
### Cooperative memory experiment
###################################################
############
## Aggregate data, calculate descriptive statistics
############
all.data$Subject <- factor(all.data$Subject)							# change Subject to factor
all.data$Correct <- as.numeric(all.data$OldState.ACC) - 1					# rescore ACC to Correct
all.data$Incorrect <- ifelse(all.data$Correct == 1, 0, ifelse(all.data$Correct == 0, 1, NA))

# Aggregrate data by subject
subj.m <- aggregate(x = all.data$Incorrect * 100, by = list(partners = all.data$GroupSize, 
  rounds = all.data$TotalRounds, subject = all.data$Subject), FUN = mean)	# calculate means per subject
subj.m$partners <- as.factor(subj.m$partners)					# change parnters to factor
subj.m$rounds <- as.factor(subj.m$rounds)					# change rounds to factor
names(subj.m) = c("partners", "rounds",  "subject", "error")			# rename subj.m columns

# Aggregate data by number of partners
subj5 <- subset(subj.m, partners == 5)						# extract rows in which partners = 5
subj10 <- subset(subj.m, partners == 10)					# extract rows in which partners = 10
subj15 <- subset(subj.m, partners == 15)					# extract rows in which partners = 15
m5 <- mean(subj5$error)								# calculate meanfor each group
m10 <- mean(subj10$error)
m15 <- mean(subj15$error)
sd5 <- sd(subj5$error)								# calculate standard deviation for each group
sd10 <- sd(subj10$error)
sd15 <- sd(subj15$error)
n5 <- length(subj5$error)							# calculate sample size for each group
n10 <- length(subj10$error)
n15 <- length(subj15$error)
ci5 <- qt(0.975,df = n5 - 1) * sd5 / sqrt(n5)					# calculate 95% confidence intervals for each group
ci10 <- qt(0.975,df = n10 - 1) * sd10 / sqrt(n10)
ci15 <- qt(0.975,df = n15 - 1) * sd15 / sqrt(n15)

# Aggregate data by sex
sex <- aggregate(all.data$Incorrect * 100, by = list(sex = all.data$Sex),
  FUN = c("mean", "sd", "count"))					# calculate mean, sd, N for each sex
male <- sex$mean[sex$sex == "male"]
female <- sex$mean[sex$sex == "female"]
male.sd <- sex$sd[sex$sex == "male"]
female.sd <- sex$sd[sex$sex == "female"]
male.n <- sex$count[sex$sex == "male"]
female.n <- sex$count[sex$sex == "female"]
male.ci <- qt(0.975, df = male.n - 1) * male.sd / sqrt(male.n)		# calculate 95% CI for each sex
female.ci <- qt(0.975, df = female.n - 1) * female.sd / sqrt(female.n)

# Aggregate data by whether correct choice was cooperate or defect
CorD <- aggregate(all.data$Incorrect * 100, by = list(action = all.data$CorrectAnswer),
  FUN = c("mean", "sd", "count"))					# calculate mean, sd, N for correct choice
coop <- CorD$mean[CorD$action == "q"]					# q = cooperate and p = defect
defect <- CorD$mean[CorD$action == "p"]
coop.sd <- CorD$sd[CorD$action == "q"]
defect.sd <- CorD$sd[CorD$action == "p"]
coop.n <- CorD$count[CorD$action == "q"]
defect.n <- CorD$count[CorD$action == "p"]
coop.ci <- qt(0.975, df = coop.n - 1) * coop.sd / sqrt(coop.n)		# calculate 95% CI for correct choice
defect.ci <- qt(0.975, df = defect.n - 1) * defect.sd / sqrt(defect.n)

# Generate forgetting function data
group5 <- subset(all.data, GroupSize == 5)			# extract rows in which GroupSize = 5
group10 <- subset(all.data, GroupSize == 10)			# extract rows in which GroupSize = 10
group15 <- subset(all.data, GroupSize == 15)			# extract rows in which GroupSize = 15
forget5 <- aggregate(100 * group5$Incorrect, by = list(events = group5$InterveningEvents), 
  FUN = c("mean", "sd", "count"), na.rm = TRUE)			# calculate means, sd, N by number of intervening events for group size 5 
names(forget5) <- c("events", "error", "sd", "N")
forget5$se <- forget5$sd/sqrt(forget5$N)			# calculate se for each number of intervening events
forget10 <- aggregate(100 * group10$Incorrect, by = list(events = group10$InterveningEvents),
  FUN = c("mean", "sd", "count"))				# calculate means, sd, N by number of intervening events for group size 10
names(forget10) <- c("events", "error", "sd", "N")
forget10$se <- forget10$sd/sqrt(forget10$N)			# calculate se for each number of intervening events
forget15 <- aggregate(100 * group15$Incorrect, by = list(events = group15$InterveningEvents),
  FUN = c("mean", "sd", "count"))			# calculate means, sd, N by number of intervening events for group size 15
names(forget15) <- c("events", "error", "sd", "N")
forget15$se <- forget15$sd/sqrt(forget15$N)			# calculate se for each number of intervening events
forget <- rbind(forget5, forget10, forget15)			# combine group size 5, 10, and 15
forget$partners <- as.factor(c(rep(5, length(forget5$events)), rep(10, length(forget10$events)),
  rep(15, length(forget15$events)))) 				# create column of group sizes
forget1015 <- subset(forget, partners != 5)			# extract rows for group size 10 or 15

# Conduct nonlinear regression to generate forgetting functions
forget1015.nls <- nls(100 - error ~ a * (1 + events) ^ b, data = forget1015,
  start = list(a = 100, b = -0.07))				# run nonlinear regression on group size 10 and 15 data
forget1015.coef <- coef(forget1015.nls)				# extract regression coefficients
forget5.nls <- nls(100 - error ~ a * (1 + events) ^ b, data = forget5,
  start = list(a = 100, b = -0.07))				# run nonlinear regression on group size 5 data
forget5.coef <- coef(forget5.nls)				# extract regression coefficients

# Aggregate data by number of rounds
round5 <- aggregate(100 * group5$Incorrect, by = list(round = group5$Round),
  FUN = c("mean", "sd", "count"), na.rm = TRUE)						# calculate mean, sd, N for round 5
names(round5) <- c("round", "error", "sd", "N")
round5$se <- round5$sd / sqrt(round5$N)							# calculate se for round 5
round10 <- aggregate(100 * group10$Incorrect, by = list(round = group10$Round),
  FUN = c("mean", "sd", "count"), na.rm = TRUE)						# calculate mean, sd, N for round 10
names(round10) <- c("round", "error", "sd", "N")
round10$se <- round10$sd / sqrt(round10$N)						# calculate se for round 10
round15 <- aggregate(100 * group15$Incorrect, by = list(round = group15$Round),
  FUN = c("mean", "sd", "count"), na.rm = TRUE)						# calculate mean, sd, N for round 15
names(round15) <- c("round", "error", "sd", "N")
round15$se <- round15$sd / sqrt(round15$N)						# calculate se for round 15
round <- rbind(round5, round10, round15)						# combine round numbers 5, 10, 15
round$partners <- as.factor(c(rep(5, length(round5$round)), rep(10, length(round10$round)),
  rep(15, length(round15$round))))							# create column of round numbers

# Aggregate data by names of partners and correlate name attractiveness with error rate
names.s <- aggregate(x = 100 * all.data$Incorrect, by = list(all.data$PartnerName, all.data$Subject), FUN = "mean", na.rm = T)
names(names.s) <- c("names", "subject", "error")
names.s$names <- names.s$names[, drop=T]
names.sm <- aggregate(x = names.s$error, by = list(names.s$names),FUN = c("mean", "sd", "length"), na.rm=T)
names(names.sm) <- c("names", "error", "sd", "N")
names.sm$sem <- names.sm$sd / sqrt(names.sm$N)
names.sm.plot <- Dotplot(reorder(names, error) ~ Cbind(error, error - sem, error + sem), data = names.sm, xlim=c(0,27))

names.m <- aggregate(x = 100 * all.data$Incorrect, by = list(all.data$PartnerName),
  FUN = mean, na.rm = T)								# calculate means per name
names(names.m) <-  c("names", "error")							# assign names to columns
names.m$attract <- c(4.13, 3.51, NA, NA, 3.89, NA, 3.92, NA, NA, NA, NA, 3.92, 3.94, 3.86, 4.09, NA, NA, NA, 4.23, 4.21, NA, NA, 
  3.44, 4.04, 3.64, 3.99, NA, 3.58, NA, NA, NA, NA, NA, 3.53, 3.35, 4.13, NA, NA, 2.93, NA)	# input name attractiveness
names.cor <- cor.test(names.m$error, names.m$attract, na.rm = TRUE)			# calculate correlation between error and attractiveness
names.r <- names.cor$estimate								# extract correlation coefficient
names(names.r) <- NULL
names.p <- names.cor$p.value								# extract p-value

# Aggregate data by images of partners and correlate image attractiveness with error rate
images <- aggregate(x = 100 * all.data$Incorrect, by = list(names = all.data$ImageFile),
  FUN = mean, na.rm = T)								# calculate means per name
names(images) = c("image", "error")							# assign names to columns
images$attract <- c(1.21, 1.08, 1.83, 0.75, 0.63, 1.08, 1.08, 0.58, 0.92, 2.13, 1.21, 1.54, 0.21, 1.71, 0.83, 0.92, 0.38, 1.71, 0.42, 
  0.96, 1.42, 1.25, 2.00, 1.54, 0.83, 1.25, 0.92, 0.92, 1.79, 1.96, 0.58, 0.79, 0.33, 1.04, 0.96, 1.50, 1.38, 1.25, 1.21, 2.04)	# input image attractiveness
images.cor <- cor.test(images$error, images$attract)					# calculate correlation between error and attractiveness
images.r <- images.cor$estimate								# extract correlation coefficient
names(images.r) <- NULL
images.p <- images.cor$p.value								# extract p-value

############
## Generate plots
############

# Set plot parameters
cb.col <- c(rgb(0, 45, 70, max = 100), rgb(0, 60, 50, max = 100), rgb(90, 60, 0, max = 100), 
  rgb(80, 40, 0, max = 100), rgb(80, 60, 70, max = 100), "black")		# create vector of color-blind-safe colors
trellis.par.set(box.rectangle = list(fill = c(cb.col[1], cb.col[3], cb.col[2]), col = "black"), clip=list(panel="off", strip="off"),
  box.umbrella = list(col = "black", lty = 1), par.ylab.text = list(cex = 1.5), par.xlab.text = list(cex = 1.5), 
  axis.text = list(cex = 1.15), strip.border = list(col = 0), layout.heights = list(strip = 2))

# Create box plot of group size and number of rounds effect on memory (Figure 2)
pdf(file = "../figures/stevens_fig2.pdf", width = 7, height = 5.5)
bplot <- bwplot(error ~ rounds | partners, data = subj.m, horizontal = F, ylab = "Memory error rate (%)", xlab = "Number of interactions",
  layout = c(3, 1), strip = strip.custom(factor.levels = c("5 Partners", "10 Partners", "15 Partners"), bg = 0,
  par.strip.text = list(cex = 1.35)),
  panel = function(x, y, groups) {
    panel.bwplot(x, y, pch = "|", horizontal = F, coef = 0, col = "black", box.ratio = 2)	# plots boxplot with medians as bars
    mean.values <<- tapply(y, x, mean)								# calculates means
    panel.points(mean.values, pch = 18, cex = 1.5, col = "black")				# plots means as diamonds
  }
)
print(bplot)
dev.off()

# Create line plot of intervening events and group size effect on memory (Figure 3)
pdf(file = "../figures/stevens_fig3.pdf", width = 7, height = 5.5)
forget.plot <- xYplot(Cbind(error, error - se, error + se) ~ events, data = forget, groups = partners, type = "b", 
  aspect = 0.75, col = cb.col, pch = c(17, 15, 16), label.curves = FALSE, cex = 1.2,
  ylim = c(-2, 40), ylab = "Memory error rate (%)", xlab = "Number of intervening events",
  key = list(corner = c(0.9, 0.05), padding.text = 3, 
    text = list(c("5 partners", "10 partners", "15 partners", "Best fit for 5 partners", "Best fit for 10 and 15 partners"), adj = 1), 
    lines = list(col = c(cb.col[1], cb.col[2], cb.col[3], 1, 1)), lwd = 2, lty = c(1, 1, 1, 2, 1),
    pch = c(17, 15, 16), type = c("b", "b", "b", "l", "l"), divide = 1),
  panel = function(...){
    panel.xYplot(...)
    panel.curve(100 - forget1015.coef[1] * (1 + x) ^ forget1015.coef[2], from = 0, type = "l", lwd = 2)	# add forgetting function for group size = 10, 15
    panel.curve(100 - forget5.coef[1] * (1 + x) ^ forget5.coef[2], from = 0, type = "l", lwd = 2, lty = 2)	# add forgetting function for group size = 5
  }
)
print(forget.plot)
dev.off()

# Create line plot of round number and group size effect on memory (Figure 4)
pdf(file = "../figures/stevens_fig4.pdf", width = 7, height = 5.5)
round.plot <- xYplot(Cbind(error, error - se, error + se) ~ round, data = round, groups = partners, type = "b", 
  aspect = 0.75, col = cb.col, pch = c(17, 15, 16), label.curves = FALSE, cex = 1.2,
  ylim = c(-2, 40), ylab = "Memory error rate (%)", xlab = "Round number",
  key = list(corner = c(0.9, 0.95), padding.text = 3, 
    text = list(c("5 partners", "10 partners", "15 partners"), adj = 1), 
    lines = list(col = cb.col[1:3]), lwd = 2, pch = c(17, 15, 16), type = "b", divide = 1)
)
print(round.plot)
dev.off()


###################################################
### Simulation
###################################################
############
## Rearrange data
############
sim.all$Error <- 100 * sim.all$Error						# convert error rates to percentages
names(sim.all) <- c("Error", "GTFT", "WSLS", "RAND", "GRIM", "CTFT", "TFT", "TF2T", "ALLC", "ALLD", "COOP")
sim.wide <- sim.all[c(3, 5, 6, 7, 10, 11)]					# extract strategies with overall success > 0.05% (ALLD, CTFT, GRIM, TFT, WSLS)
sim.df <- stack(sim.wide)							# stack columns to generate a single column
sim.df$Error <- rep(0:50, 6) 							# create column of error rates
names(sim.df) <- c("Wins", "Strategy", "Error")

############
## Generate plot
############

# Set plot parameters
sim.col <- c(rgb(0, 45, 70, max = 100), "black", rgb(0, 60, 50, max = 100), rgb(90, 60, 0, max = 100), 
  rgb(80, 40, 0, max = 100), rgb(80, 60, 70, max = 100))		# create vector of strategy colors
sim.names <- c("ALLD", "CTFT", "GRIM", "TFT", "WSLS", "Cooperation")		# create vector of strategy names

# Create line plot of error rate effect on strategies' success (Figure 6)
pdf(file = "../figures/stevens_fig5.pdf", width = 7, height = 5.5)
simplot <- xyplot(Wins ~ Error, data = sim.df, groups = Strategy, type = "l",
  col = sim.col, lwd = c(3, 4, 3, 3, 3, 3), lty = c(1, 1, 2, 3, 4, 5),
  xlab = "Memory error rate (%)", ylab = "Mean percentage of replicates won",
  key = list(corner = c(0.9, 0.5), padding.text = 3, text = list(sim.names, adj = 1), 
    lines = list(col = c(sim.col[1], sim.col[3], sim.col[4], sim.col[5], sim.col[6],
    sim.col[2])), lwd = c(3, 3, 3, 3, 3, 4), lty = c(1, 2, 3, 4, 5, 1)),
  panel = function(...){
    panel.rect(9.6, -7, 24, 107, col = "grey85", border = NULL)			# add shaded box for empirical error rates
    panel.xyplot(...)
  }
)
print(simplot)
dev.off()


###################################################
### Game theoretical analysis
###################################################
############
## Calculate memory-1 strategy payoffs
############

# Determine each memory-1 strategy's probability of cooperating following T, R, P, and S for every other strategy
allc <- matrix(c(1,1,1,1,1,1,1,1,1,1,  			#ALLC vs ALLC (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,1,1,1,1,1,  				#ALLD vs ALLC
	1,0,1,0,0,1,1,1,1,1,  				#GRIM vs ALLC
	1,1,1,0.33,0.33,1,1,1,1,1,  			#gTFT vs ALLC
	0.5,0.5,0.5,0.5,0.5,1,1,1,1,1,  		#RAND vs ALLC
	1,1,1,0,0,1,1,1,1,1,  				#TFT vs ALLC
	1,0,1,1,0,1,1,1,1,1),7,byrow=T) 		#WSLS vs ALLC
alld <- matrix(c(1,1,1,1,1,0,0,0,0,0,  			#ALLC vs ALLD (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,0,0,0,0,0,  				#ALLD vs ALLD
	1,0,1,0,0,0,0,0,0,0,  				#GRIM vs ALLD
	1,1,1,0.33,0.33,0,0,0,0,0,  			#gTFT vs ALLD
	0.5,0.5,0.5,0.5,0.5,0,0,0,0,0, 		 	#RAND vs ALLD
	1,1,1,0,0,0,0,0,0,0,  				#TFT vs ALLD
	1,0,1,1,0,0,0,0,0,0),7,byrow=T)  		#WSLS vs ALLD
grim <- matrix(c(1,1,1,1,1,1,0,1,0,0,  			#ALLC vs GRIM (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,1,0,1,0,0,  				#ALLD vs GRIM
	1,0,1,0,0,1,0,1,0,0,  				#GRIM vs GRIM
	1,1,1,0.33,0.33,1,0,1,0,0,  			#gTFT vs GRIM
	0.5,0.5,0.5,0.5,0.5,1,0,1,0,0,  		#RAND vs GRIM
	1,1,1,0,0,1,0,1,0,0,  				#TFT vs GRIM
	1,0,1,1,0,1,0,1,0,0),7,byrow=T) 		#WSLS vs GRIM
gtft <- matrix(c(1,1,1,1,1,1,1,1,0.33,0.33,  		#ALLC vs gTFT (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,1,1,1,0.33,0.33,  			#ALLD vs gTFT
	1,0,1,0,0,1,1,1,0.33,0.33,  			#GRIM vs gTFT
	1,1,1,0.33,0.33,1,1,1,0.33,0.33,  		#gTFT vs gTFT
	0.5,0.5,0.5,0.5,0.5,1,1,1,0.33,0.33,  		#RAND vs gTFT
	1,1,1,0,0,1,1,1,0.33,0.33,  			#TFT vs gTFT
	1,0,1,1,0,1,1,1,0.33,0.33),7,byrow=T)		#WSLS vs gTFT
rand <- matrix(c(1,1,1,1,1,0.5,0.5,0.5,0.5,0.5, 	#ALLC vs RAND (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,0.5,0.5,0.5,0.5,0.5,  		#ALLD vs RAND
	1,0,1,0,0,0.5,0.5,0.5,0.5,0.5,  		#GRIM vs RAND
	1,1,1,0.33,0.33,0.5,0.5,0.5,0.5,0.5,  		#gTFT vs RAND
	0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,	#RAND vs RAND
	1,1,1,0,0,0.5,0.5,0.5,0.5,0.5,  		#TFT vs RAND
	1,0,1,1,0,0.5,0.5,0.5,0.5,0.5),7,byrow=T)	#WSLS vs RAND
tft <- matrix(c(1,1,1,1,1,1,1,1,0,0,  			#ALLC vs TFT (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,1,1,1,0,0,  				#ALLD vs TFT
	1,0,1,0,0,1,1,1,0,0,  				#GRIM vs TFT
	1,1,1,0.33,0.33,1,1,1,0,0,  			#gTFT vs TFT
	0.5,0.5,0.5,0.5,0.5,1,1,1,0,0,  		#RAND vs TFT
	1,1,1,0,0,1,1,1,0,0,  				#TFT vs TFT
	1,0,1,1,0,1,1,1,0,0),7,byrow=T)  		#WSLS vs TFT
wsls <- matrix(c(1,1,1,1,1,1,0,1,1,0,  			#ALLC vs WSLS (c1,t1,r1,p1,s1,c2,t2,r2,p2,s2)
	0,0,0,0,0,1,0,1,1,0,  				#ALLD vs WSLS
	1,0,1,0,0,1,0,1,1,0,  				#GRIM vs TFT
	1,1,1,0.33,0.33,1,0,1,1,0,  			#gTFT vs WSLS
	0.5,0.5,0.5,0.5,0.5,1,0,1,1,0,  		#RAND vs WSLS
	1,1,1,0,0,1,0,1,1,0,  				#TFT vs WSLS
	1,0,1,1,0,1,0,1,1,0),7,byrow=T)  		#WSLS vs WSLS
strat <- list(allc = allc, alld = alld, grim = grim, gtft = gtft, rand = rand, tft = tft, wsls = wsls)	# create list including all strategy probabilities
T <- 5						# define payoff matrix
R <- 3
P <- 1
S <- 0
alpha <- 0.9 					# probability of future interaction
epsilon2 <- c(0:50)/100				# vector of error rates
mem1.payoffs <- 0				# initiate variable

# Generate the expected payoff for the vector of error rates
payoff.calc <- function(epsilon2) {			
  if (tt1 == rr1 && rr1 == pp1 && pp1 == ss1) {		# if the prob(coop) of the invading strategy is the same (ALLC, ALLD, RAND)...
    epsilon1 <- 0					# ...the strategy has no error rate
    if (tt2 == rr2 && rr2 == pp2 && pp2 == ss2) {	# if the prob(coop) of the resident strategy is the same (ALLC, ALLD, RAND)...
      epsilon2 <- 0					# ...the strategy has no error rate
    }
  }
  else {	
    epsilon1 <- epsilon2 				# otherwise, the error rates are the same for the two strategies
    if (tt2 == rr2 && rr2 == pp2 && pp2 == ss2) {	# but if the prob(coop) of the resident strategy is the same (ALLC, ALLD, RAND)...
      epsilon2 <- 0					# ...the strategy has no error rate
    }
  }
  t1 <- tt1 * (1 - epsilon1) + (1 - tt1) * epsilon1	# assign error probabilities to prob(coop)
  t2 <- tt2 * (1 - epsilon2) + (1 - tt2) * epsilon2
  r1 <- rr1 * (1 - epsilon1) + (1 - rr1) * epsilon1
  r2 <- rr2 * (1 - epsilon2) + (1 - rr2) * epsilon2
  p1 <- pp1 * (1 - epsilon1) + (1 - pp1) * epsilon1
  p2 <- pp2 * (1 - epsilon2) + (1 - pp2) * epsilon2
  s1 <- ss1 * (1 - epsilon1) + (1 - ss1) * epsilon1
  s2 <- ss2 * (1 - epsilon2) + (1 - ss2) * epsilon2
  M <- matrix(c(s2 * (1 - t1), (1 - r1) * r2, (1 - p1) * p2, (1 - s1) * t2,	# generate transition matrix
    s2 * t1, r1 * r2, p1 * p2, s1 * t2, (1 - s2) * (1 - t1),
    (1 - r1) * (1 - r2), (1 - p1) * (1 - p2), (1 - s1) * (1 - t2),
    (1 - s2) * t1, r1 * (1 - r2), p1 * (1 - p2), s1 * (1 - t2)), 4, byrow = T)
  y <- c((1 - c1) * c2, c1 * c2,(1 - c1)*(1 - c2), (1 - c2) * c1);		# generate vector of initial cooperation
  v <- c(T, R, P, S);								# initiate payoff matrix
  iden <- diag(4);								# generate identity matrix
  s <- solve(iden - alpha * M);							# find geometric series
  E <<- v %*% s %*% y;								# calculate expected benefit
}

# Calculate expected payoff for all strategies
strat.payoff.calc <- function(strat) {			
  c1 <<- strat[1]					# assign initial prob(coop) from strat matrix
  c2 <<- strat[6]
  tt1 <<- strat[2]					# assign prop(coop) from strat matrix
  rr1 <<- strat[3]
  pp1 <<- strat[4]
  ss1 <<- strat[5]
  tt2 <<- strat[7]
  rr2 <<- strat[8]
  pp2 <<- strat[9]
  ss2 <<- strat[10]
  EE <<- apply(as.matrix(epsilon2), 1, payoff.calc)	# apply payoff.calc to all invading strategies
}

# Apply strat.payoff.calc to all resident strategies
for (i in 1:length(strat)) {				
  payoffs <- apply(strat[[i]], 1, strat.payoff.calc)
  mem1.payoffs <- rbind(mem1.payoffs, payoffs)		# create payoff table
}

dim.payoffs <- dim(mem1.payoffs)			# Get dimensions of payoff table
mem1.payoffs <- mem1.payoffs[2:dim.payoffs[1],]		# Remove initiating row
mem1.payoffs <- data.frame(mem1.payoffs)		# Convert to data frame
mem1.payoffs$err.rate <- rep(c(0:50),7)			# Create column of error rates
mem1.payoffs$against <- c(rep("allc", 51), rep("alld", 51), rep("grim", 51), 	#Create column of strategy names
  rep("gtft", 51), rep("rand", 51), rep("tft", 51), rep("wsls", 51))
names(mem1.payoffs) <- c("ALLC", "ALLD", "GRIM", "GTFT", "RAND", "TFT", "WSLS", "err.rate", "against")	# Create vector of strategy names

############
## Simulate memory-2 strategy payoffs
############

# Input memory-2 strategy data from files
path <- "data/other_payoffs/"					# set path to memory-2 strategy simulation data
files <- list.files(path)					# get file names from current working directory
filenames <- paste(path, files, sep = "")
ct.pay <- 0
rows <- (c(1:51) * 3) - 1					# set number of rows in an error rate block
L <- lapply(filenames, function(x) {				# for each file in filenames...
  X <<- readLines(con = x)					# read in data lines
  payrows <<- X[rows]
  payoffs <<- as.numeric(substr(payrows,17,23))			# find number of replicates for this file
  ct.pay <<- c(ct.pay, payoffs)
})
ct.pay <- ct.pay[2:length(ct.pay)]				# cut initial row of zeros
ct.pay <- data.frame(ct.pay * 10)				# convert to data frame and to total payoffs (by multiplying by 10 interactions)
ct.pay$err.rate <- c(rep(0:50, 32))				# create error rate vector
ct.pay$against <- c(rep("allc", 51), rep("allc", 51), rep("alld", 51), rep("alld", 51), rep("ctft", 51), rep("ctft", 51), 	# create vector of 'strategies against'
  rep("ctft", 51), rep("ctft", 51), rep("ctft", 51), rep("ctft", 51), rep("ctft", 51), rep("ctft", 51), rep("ctft", 51), 
  rep("grim", 51), rep("grim", 51), rep("gtft", 51), rep("gtft", 51), rep("rand", 51), rep("rand", 51), rep("tf2t", 51), 
  rep("tf2t", 51), rep("tf2t", 51), rep("tf2t", 51), rep("tf2t", 51), rep("tf2t", 51), rep("tf2t", 51), rep("tf2t", 51), 
  rep("tf2t", 51), rep("tft", 51), rep("tft", 51), rep("wsls", 51), rep("wsls", 51))
ct.pay$strat1 <- c(rep("CTFT", 51), rep("TF2T", 51), rep("CTFT", 51), rep("TF2T", 51), rep("ALLC", 51), rep("ALLD", 51), 	# create vector of 'strategies for'
  rep("CTFT", 51), rep("GRIM", 51), rep("GTFT", 51), rep("RAND", 51), rep("TF2T", 51), rep("TFT", 51), rep("WSLS", 51), 
  rep("CTFT", 51), rep("TF2T", 51), rep("CTFT", 51), rep("TF2T", 51), rep("CTFT", 51), rep("TF2T", 51), rep("ALLC", 51), 
  rep("ALLD", 51), rep("CTFT", 51), rep("GRIM", 51), rep("GTFT", 51), rep("RAND", 51), rep("TF2T", 51), rep("TFT", 51), 
  rep("WSLS", 51), rep("CTFT", 51), rep("TF2T", 51), rep("CTFT", 51), rep("TF2T", 51))
CTFT.payoffs <- data.frame(ct.pay[(ct.pay$strat1 == "CTFT"), ])		# extract data for CTFT for (column of data to be added to all.payoffs)
names(CTFT.payoffs) <- c("CTFT", "err.rate", "against", "strat1")
TF2T.payoffs <- data.frame(ct.pay[(ct.pay$strat1 == "TF2T"), ])		# extract data for TF2T for (column of data to be added to all.payoffs)
names(TF2T.payoffs) <- c("TF2T", "err.rate", "against", "strat1")

# Create matrix of payoffs with 'for strategies' as columns and 'against strategies' as rows (to be appended to mem1.payoffs to generate all.payoffs)
ctft.vec <- ct.pay[(ct.pay$against == "ctft" & ct.pay$strat1 != "CTFT" & ct.pay$strat1 != "TF2T"), 1]	# extract data for 'ctft against'
ctft.mat <- matrix(ctft.vec, nrow = 51)									# rearrange into matrix
ctft.df <- data.frame(ctft.mat)										# convert to data frame
names(ctft.df) <- c("ALLC", "ALLD", "GRIM", "GTFT", "RAND", "TFT", "WSLS")	
ctft.df$err.rate <- c(0:50)										# create vector of error rates
ctft.df$against <- c(rep("ctft", 51))									# create vector of 'strategies against'
tf2t.vec <- ct.pay[(ct.pay$against == "tf2t" & ct.pay$strat1 != "CTFT" & ct.pay$strat1 != "TF2T"), 1]	# extract data for 'tf2t against'
tf2t.mat <- matrix(tf2t.vec, nrow = 51)									# rearrange into matrix
tf2t.df <- data.frame(tf2t.mat)										# convert to data frame
names(tf2t.df) <- c("ALLC", "ALLD", "GRIM", "GTFT", "RAND", "TFT", "WSLS")
tf2t.df$err.rate <- c(0:50)										# create vector of error rates
tf2t.df$against <- c(rep("tf2t", 51))									# create vector of 'strategies against'
mem2.payoffs <- rbind(ctft.df, tf2t.df)									# append ctft.df with tf2t.df

mem12.payoffs <- rbind(mem1.payoffs, mem2.payoffs)							# append mem1.payoffs with mem2.payoffs
strat.ord <- order(mem12.payoffs$against)								# create order by 'strategies against'
all.payoffs <- mem12.payoffs[strat.ord, ]								# create all.payoffs ordered by 'strategies against'
all.payoffs$CTFT <- CTFT.payoffs$CTFT									# add CTFT column
all.payoffs$TF2T <- TF2T.payoffs$TF2T									# add TF2T column
all.payoffs <- all.payoffs[, c(1, 2, 10, 3, 4, 5, 11, 6, 7, 8, 9)]					# reorder by 'strategies for'

############
## Generate plot
############

# Set plot parameters
strat.names <- c("ALLC", "ALLD", "CTFT", "GRIM", "GTFT", "RAND", "TF2T", "TFT", "WSLS")			# generate list of 'strategy against' names
strat.col <- c("lightgreen", rgb(0, 114, 178, max = 255), rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255), 	# Create vector of strategy colors to match simulation colors
  rgb(86, 180, 233, max = 255), rgb(240, 228, 66, max = 255), rgb(213, 94, 0, max = 255), "black", rgb(204, 121, 167, max = 255))


# Create line plot of payoffs for each strategy against each other strategy (Figure 7)
pdf(file = "../figures/stevens_fig6.pdf", width = 7, height = 5.5)
simplot <- xyplot(ALLC + ALLD + CTFT + GRIM + GTFT + RAND + TF2T + TFT + WSLS ~ err.rate | against, data = all.payoffs, type = "l", 
  layout = c(3,3), as.table = T, aspect = 0.9,  col = strat.col, lwd = 1.5, lty = c(1:9), 
  ylab = "Payoffs", xlab = "Memory error rate (%)", between = list(y = 0.5),
  scales = list(y = list(relation = "free", cex = 0.8), x = list(cex = 0.8)), 
  key = list(x = 1, y = 0.5, corner = c(0, 0), text = list(strat.names), lines = list(col = strat.col, lty = c(1:9), lwd = 2)),
  strip = strip.custom(factor.levels = strat.names), par.settings = list(strip.border = list(col = 1), layout.heights = list(strip = 1))
)
print(simplot)
dev.off()



