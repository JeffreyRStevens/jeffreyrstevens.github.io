###################################################
### stevens.2010.IJCP.rcode.R
### Created by Jeffrey R. Stevens on 16 April 2010 (jstevens@mpib-berlin.mpg.de)
### Summary: This script calculates descriptive statistics and generates figures
### 	for the analysis of delayed gratification data in bonobos:
###	Stevens, J.R., Rosati, A.G., Heilbronner, S.R., & Schmücking, N. forthcoming. 
###	 Waiting for grapes: expectancy and delayed gratification in bonobos.  
###	 International Journal of Comparative Psychology.
### Instructions: Place this file and the data file (stevens.2010.IJCP.data.txt)
### 	in the same directory.  In R set the working directory to this directory.  Type 
### 	> source("stevens.2010.IJCP.rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PDF versions of the figures in the figures directory.
### Data file: Description of the data columns:
### 	subject - name of subject
### 	date - date of testing
### 	phase - phase number (1, 2, 3)
### 	cond - experimenter condition (control, high reliability, low reliability)
### 	replicate - session number within each phase and condition (phases 1 and 2
###		have one session each, phase 3 has three sessions)
### 	exptr_order - code for order of experimenter condition
### 	session - code for session number
### 	test_session - code for test session (high or low reliability) number
### 	trial - trial number
### 	stop - number of grapes placed before stopping (main dependent variable)
### 	interrupted - code for interruptions: "no_interruption" if the trial did not 
###		have interruptions, otherwise the grape placement number at which
###		the trial was interrupted.
### Uses: This script can be reproduced and modified for personal and scientific use.
###################################################


###################################################
### Load bonobo data
###################################################
data.all <- read.csv("stevens.2010.IJCP.data.txt")
library(plyr)       							#required for ddply
library(lattice)        						#required for bwplot, xyplot


###################################################
### Calculates descriptive statistics
###################################################
## Calculates means, medians, std. dev, and 95% CI for each subject for each condition and phase
subj.means <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase), subj = data.all$subject), mean) 
subj.medians <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase), subj = data.all$subject), median) 
subj.sd <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase), subj = data.all$subject), sd) 
subj.N <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase), subj = data.all$subject), length) 
subj.ci95 <- qt(0.975, df = subj.N$x - 1) * subj.sd$x / sqrt(subj.N$x)	#calculate 95% confidence intervals
subj.data<-cbind(subj.medians, subj.means$x, subj.sd$x, subj.ci95, subj.N$x)
names(subj.data) <- c("cond", "phase", "subj", "median", "mean", "sd", "ci95", "N")

## Calculates overall means, medians, std. dev, and 95% CI for each condition and phase
overall.means <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase)), mean) 
overall.medians <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase)), median) 
overall.sd <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase)), sd) 
overall.N <- aggregate(data.all$stop, list(cond = data.all$cond, phase = as.factor(data.all$phase)), length) 
overall.ci95 <- qt(0.975, df = overall.N$x - 1) * overall.sd$x / sqrt(overall.N$x)	#calculate 95% confidence intervals
overall.data<-cbind(overall.medians, overall.means$x, overall.sd$x, overall.ci95, overall.N$x)
names(overall.data) <- c("cond", "phase", "median", "mean", "sd", "ci95", "N")

## Calculates descriptive statistics on data by condition only for phase 3 
ph3.all <- data.all[data.all$phase == "Phase 3", ] 				#extract phase 3 data
ph3 <- (ddply(ph3.all, c("cond","replicate"), function(df)
    return(c(stop = mean(df$stop), SD = sd(df$stop), N = length(df$stop))))) 	#calculate mean, SD, and N
ci95 <- qt(0.975, df = ph3$N - 1) * ph3$SD / sqrt(ph3$N) 			#calculate 95% confidence intervals
ph3.desc <- cbind(ph3, ci95)
mph1 <- round(ph3.desc$stop[ph3.desc$cond == "Low" & ph3.desc$replicate == 1], 1) 
mph2 <- round(ph3.desc$stop[ph3.desc$cond == "Low" & ph3.desc$replicate == 2], 1) 
mph3 <- round(ph3.desc$stop[ph3.desc$cond == "Low" & ph3.desc$replicate == 3], 1) 
ciph1 <- round(ph3.desc$ci95[ph3.desc$cond == "Low" & ph3.desc$replicate == 1], 1) 
ciph2 <- round(ph3.desc$ci95[ph3.desc$cond == "Low" & ph3.desc$replicate == 2], 1) 
ciph3 <- round(ph3.desc$ci95[ph3.desc$cond == "Low" & ph3.desc$replicate == 3], 1) 

## Calculates expected values for phase 3
stop8 <- 19     							#assign the number of actual interruptions for each step
stop9 <- 11
stop10 <- 18
pstop8 <- stop8 / (stop8 + stop9 + stop10)   				#calculate probability of no interruption at each step
pstop9 <- stop9/(stop9 + stop10)
pstop10 <- 0
ev <- c(1:10)           						#initiate array of payoffs with no interruption
ev[8] <- (1 - pstop8) * 8 + pstop8 * ev[7 - 5]      			#calculate expected value
ev8 <- round(ev[8], 2)
ev[9] <- (1 - pstop8) * (1 - pstop9) * 9 + (1 - pstop8) * pstop9 * ev[8 - 5] + pstop8 * ev[7 - 5]
ev9 <- round(ev[9], 2)
ev[10] <- 4     							#expected value when interrupted at step 10


###################################################
### Generate Figure 2
###################################################
pdf(file = "fig2.pdf", width = 6, height = 5)
fig2 <- bwplot(stop ~ cond | phase, data = data.all, layout = c(3,1), pch = "|",	#plot median as line
  ylab = "Number of grapes before stopping", ylim = c(-1:11),
  panel = function(x, y, groups) {
    panel.bwplot(x,y , pch = "|", horizontal = F, coef = 0)     			#plot boxplot
    mean.values <<- tapply(y, x, mean)          					#calculates means
    panel.points(mean.values, pch = 18, cex = 1.5, col = "black")       		#plot means as diamonds
  })
trellis.par.set(box.umbrella = list(lty = 1, col = "black"), box.rectangle = list(col = "black")
)
print(fig2)
dev.off()

###################################################
### Generate Figure 3
###################################################
data.stop <- aggregate(x = data.all$stop, by = list(subject = data.all$subject, phase = as.factor(data.all$phase), cond = data.all$cond, 
  trial = data.all$trial), mean)
my.col <- c(rgb(0, 45, 70, max = 100), rgb(90, 60, 0, max = 100), rgb(0, 60, 50, max = 100))	#assign line colors
my.pch <- c(15, 17, 19)										#assign symbols
pdf(file = "fig3.pdf", width = 6, height = 5)
fig3 <- xyplot(x ~ trial | phase * subject,
  data = data.stop, groups = cond, type = "b", aspect = 1, pch = my.pch, col = my.col,
  xlab = "Trial Number", ylab = "Number of grapes before stopping",
  ylim = c(-1:11),
  key = list(space = "right", text = list(levels(data.stop$cond)), points = list(pch = my.pch, col = my.col)), 
  scales = list(alternating = 1, x = (list(limits = 0:9, tck = c(0.5,0),
  at = 1:8)), y = list(limits = 0:11, at = 1:10, tck = c(0.5,0.5))),
  panel = function(...) {
    panel.number = panel.number()
    if (panel.number == 3 | panel.number == 6) {    			#for panels 3 and 6
      panel.abline(h = 7, lty = 2)          				#draw line at y=7
   }
    panel.xyplot(...)                   				#plot lines
  }
)
print(fig3)
dev.off()

