###################################################
### stevens_etal_rcode.R
### Created by Jeffrey R. Stevens on 2014-04-02 (jeffrey.r.stevens@gmail.com),
###  finalized on 2015-09-09
### Summary: This script analyzes and produces plots for caching and operant
### 	delay choice data pinyon jays.
### Instructions: Place this file and the data files in the main directory.
### 	Create a folder called "figures". Set the R working directory to the 
### 	main directory.  At the R command prompt, type
### 	> source("stevens_etal_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving figures in the figures directory.
### Uses: This script can be reproduced and modified for personal and scientific use.
### Data files:
###  stevens_operant_data.csv
### 	experiment - experiment number (1, 2, or 3 [pellet/pine nut comparison])
### 	subject - subject ID
### 	date - date (YYYY-MM-DD)
### 	choiceLL - choice for smaller, sooner option (0) or larger, later option (1)
### 	rt - reaction time in milliseconds
### 	delay - length of large delay in seconds(
###   condition - experimental condition: for experiment 2, short-day (10:14 h) or long-day (10:14 h)
###     condition; for experiment 3, pellets or pine nuts
###  stevens_cache_data.csv
### 	experiment - experiment number (1, 2, or 3 [pellet/pine nut comparison])
### 	date - date (YYYY-MM-DD)
### 	subject - subject ID
### 	cached - number items cached in tray at end of session
###   condition - experimental condition: for experiment 2, short-day (10:14 h) or long-day (10:14 h)
###     condition; for experiment 3, pellets or pine nuts; for experiment 4, prefed or not prefed
###################################################

###################
## Load libraries and create functions
###################
rm(list=ls())  # clear all variables

library(epicalc)      # needed for ci
library(foreach)      # needed for foreach loops
library(lattice)      # needed for lattice graphics
library(latticeExtra) # needed for layer
library(xtable)       # needed for xtable

# create function analogous to %in% that searches for items not in a vector
"%notin%" <- function(x, table) {
  match(x, table, nomatch = 0) == 0
}

###################
## Experiment 1: Correlating caching and operant delay choice data
###################
########
# Operant data
########
# Load and clean data
operant_data <- read.csv("stevens_operant_data.csv")	# read in operant delay choice data
operant_data1 <- subset(operant_data, experiment == 1)  	# subset operant data for experiment 1 with complete sessions
operant_data1$date <- as.Date(operant_data1$date)					# convert date column to date format
operant_data1$subject <- as.factor(operant_data1$subject)	# convert subject column to factor

subjects <- unique(operant_data1$subject) 								# create vector of subject ids
num_subjects <- length(subjects)  												# calculate number of subjects

## Aggregate trials to one row per session
# Prepare data frame and counter
session_delay1 <- data.frame(matrix(rep(NA, 4), ncol = 4)) 	# initiate data frame for session_delay
names(session_delay1) <- c("subject", "date", "delay", "session_counter")
session_counter1 <- 1  	# initiate session_counter

# Extract only one row per session
for(i in 1:length(operant_data1[, 1])) { 					# for each row of data
  if(i > 1 & i < length(operant_data1[, 1])) { 		 # if it is neither the first or last row
    if(operant_data1$subject[i] != operant_data1$subject[i + 1]) {  # if this is not the last row for this subject
      current_session1 <- cbind(operant_data1[i, c(2:3, 6)], session_counter1)  # add session counter column
      names(session_delay1) <- names(current_session1)  						# sync column names
      session_delay1 <- rbind(current_session1, session_delay1)  		# combine session means with current session data
      session_counter1 <- 1  																				# initiate session counter
    } else {  		# if the last row for this subject
      if(operant_data1$date[i] != operant_data1$date[i + 1]) {  		# if the date does not change
        current_session1 <- cbind(operant_data1[i, c(2:3, 6)], session_counter1)  # add session counter column
        names(session_delay1) <- names(current_session1)
        session_delay1 <- rbind(current_session1, session_delay1)  	# combine current session data to all data
        session_counter1 <- session_counter1 + 1  									# increment session counter
      }
    }
  }
  if(i == length(operant_data1[, 1])) { 		# if the last row
    current_session1 <- cbind(operant_data1[i, c(2:3, 6)], session_counter1)  # add session counter column
    names(session_delay1) <- names(current_session1)
    session_delay1 <- rbind(current_session1, session_delay1)  			# combine current session data to all data
  }
}

# Clean and reorder data
session_delay1 <- session_delay1[-length(session_delay1[, 1]), ] 		# remove last row of NAs
session_delay1 <- session_delay1[order(session_delay1$subject, session_delay1$session_counter), ] # reorder based on subject and session_counter

# Calculate means
mean_rt1 <- aggregate(rt ~ subject * date, data = operant_data1, FUN = "mean")		# calculate mean reaction time per subject per session
session_means1 <- merge(session_delay1, mean_rt1)										# merge dely and reaction time data

## Generate columns of moving averages of last five and previous five sessions
# Prepare vectors and counter
session_means1$mean_diff <- session_means1$last_five <- session_means1$first_five <- rep(NA, length(session_means1[, 1])) # initiate columns of sessions_means
subject_counter1 <- 1  # initiate subject counter (counts numbers of sessions for a subject)

# Calculate moving averages of last five and previous five sessions
for(i in 1:length(session_means1[, 1])){ 			# for each row in session_means
  if(i != 1 & i != length(session_means1[, 1])) { 		# if it is not the first or last line
    if(session_means1[i, 1] != session_means1[i-1, 1]) { 		# if the subject changes
      subject_counter1 <- 1  # initiate subject_counter
    } else {  # if subject stays the same
      if(subject_counter1 >= 5) {  	# if subject counter (session number) is greater than or equal to 5 (at least five sessions have occurred)
        last_five1 <- mean(session_means1[(i - 4):i, ]$delay) 	# calculate mean delay in last five sessions
        session_means1[i, ]$last_five <- last_five1							# assign mean delay of last five sessions
        subject_counter1 <- subject_counter1 + 1								# increment counter
        if(subject_counter1 >= 10) { 		# if there are at least 10 rows of data for this subject  (at least ten sessions have occurred)
          first_five1 <- mean(session_means1[(i - 9):(i - 5), ]$delay) 	# calculate mean delay in previous five sessions
          #last_five1 <- mean(session_means1[(i - 4):i, ]$delay)
          if(first_five1 < last_five1) {  	# if first_five is less than last_five
            mean_diff1 <- first_five1 / last_five1							# calculate proportional difference between last and previous five
          } else {  												# if last_five is less than first_five
            mean_diff1 <- last_five1 / first_five1							# calculate proportional difference between last and previous five
          }
          session_means1[i, ]$first_five <- first_five1					# assign mean delay of previous five sessions
          #session_means1[i, ]$last_five <- last_five1
          session_means1[i, ]$mean_diff <- mean_diff1						# assign proportional difference
        }
      } else {		# if subject counter (session number) is less than 5
        subject_counter1 <- subject_counter1 + 1								# increment counter
      }
    }
  } else { 		# if it is the first or last line
    if(subject_counter1 >= 10) { 	# if there are at least 10 rows of data for this subject
      first_five1 <- mean(session_means1[(i - 9):(i - 5), ]$delay) 	# calculate mean delay in last five sessions
      last_five1 <- mean(session_means1[(i - 4):i, ]$delay) 				# calculate mean delay in previous five sessions
      if(first_five1 < last_five1) {  	# if first_five is less than last_five
        mean_diff1 <- first_five1 / last_five1											# calculate proportional difference between last and previous five
      } else {													# if last_five is less than first_five
        mean_diff1 <- last_five1 / first_five1											# calculate proportional difference between last and previous five
      }
      session_means1[i, ]$first_five <- first_five1									# assign mean delay of last five sessions
      session_means1[i, ]$last_five <- last_five1										# assign mean delay of previous five sessions
      session_means1[i, ]$mean_diff <- mean_diff1										# assign proportional difference
      subject_counter1 <- subject_counter1 + 1											# increment counter
    } else {		# if there are not 10 rows of data for this subject
      subject_counter1 <- subject_counter1 + 1											# increment counter
    }
  }
}

## Find indifference points for each subject
#Prepare data frames and counter
session_stop1 <- data.frame(matrix(rep(NA, 8), ncol = 8))  # initiate dataframe
names(session_stop1) <- c("subject", "date", "delay", "session", "rt", "first_five", "last_five", "mean_diff")  # name columns
indiff_means1 <- data.frame(subjects, rep(NA, num_subjects), rep(NA, num_subjects), rep(NA, num_subjects))  # initiate dataframe
names(indiff_means1) <- c("subject", "start_date", "end_date", "delay_indiff")  # name columns
subject_counter1 <- 1 # initiate subject counter

#	Calculate indifference points and create data frame of indifference points for each subject
foreach(curr_subj = subjects) %do% {  # for each subject
  subj_subset1 <- subset(session_means1, subject == curr_subj)  # subset this subject's data
  subject_counter1 <- subject_counter1 + 1  	# increment subject counter
  subj_length1 <- length(subj_subset1[, 1]) 	# calculate number of rows
  stop_counter1 <- 1 					# initiate stop counter
  stop_flag1 <- 0  						# initiate stop flag
  for(i in 1:subj_length1) { 	# for each row in the subject's data
    subj_thusfar1 <- subj_subset1[1:i, ]  # add rows to list of rows included thus far
    if(!is.na(subj_subset1[i, ]$mean_diff)) {  # if mean_diff is not NA
      # test whether stop criteria have been met
      if(min(subj_thusfar1$delay, na.rm = TRUE) %notin% subj_thusfar1[(i-4):i, ]$delay & # if not minimum delay so far
           subj_thusfar1[i, ]$last_five != max(subj_thusfar1$last_five, na.rm = TRUE) & # if mean of last five is not highest mean of five
           subj_thusfar1[i, ]$mean_diff >= 0.9 &  # if ratios of last and previous five delays is greater than or equal to 90%
           stop_flag1 == 0) {  # if the stop flag has not been triggered already
          stop_counter1 <- i    # assign stop counter to current row number
        stop_flag1 <- 1    # raise stop flag
      }
    }
  }
  if(stop_flag1 == 0) {  # if stop flag has not been raised
    stop_counter1 <- i   # assign stop counter to current row number
  }
  subject_complete1 <- subj_subset1[1:stop_counter1, ] 			# subset subject data until stopping point
  names(session_stop1) <- names(subject_complete1)
  session_stop1 <- rbind(session_stop1, subject_complete1) 	# append subject data until stopping point
  indiff_means1[which(indiff_means1$subject == curr_subj), ]$delay_indiff <- session_stop1[length(session_stop1[, 1]), ]$last_five # find mean delay at indifference
  indiff_means1[which(indiff_means1$subject == curr_subj), ]$start_date <- session_stop1[(length(session_stop1[, 1]) - 10), ]$date # find start date for indifference
  indiff_means1[which(indiff_means1$subject == curr_subj), ]$end_date <- session_stop1[length(session_stop1[, 1]), ]$date # find stop date for indifference
}

# Clean data
session_stop1 <- session_stop1[-1, ]  # remove NA
session_stop1$date <- as.Date(session_stop1$date, origin = "1970-01-01")
indiff_means1$start_date <- as.Date(indiff_means1$start_date, origin = "1970-01-01")
indiff_means1$end_date <- as.Date(indiff_means1$end_date, origin = "1970-01-01")

# Calculate max and min sessions to indifference points
min_indiff_sessions1 <- min(aggregate(session_counter1 ~ subject, data = session_stop1, FUN = "max")$session_counter)	# calculate minimum number of sessions to indifference
max_indiff_sessions1 <- max(session_stop1$session_counter1)		# calculate maximum number of sessions to indifference

## Plot titration indifference point data for each subject
stopped_plot1 <- xyplot(delay ~ session_counter1, groups = subject, data = session_stop1,    # plot delay to large at indifference by date for each subject
  type = "b", pch = c(0:7, 15:18),
  xlab = "Session", ylab = "Delay to large reward (s)",
  strip = strip.custom(factor.levels = levels(operant_data1$subject), par.strip.text = list(cex = 1.25)),
  par.settings = list(axis.text = list(cex = 1.15), par.xlab.text = list(cex = 2), par.ylab.text = list(cex = 2), layout.heights = list(strip = 1.25))
)
cairo_pdf(file = "figures/expt1_subj_indiff_points.pdf", width = 8, height = 6)
plot(stopped_plot1)
dev.off()

########
# Caching data
########
## Load and clean data
cache_data <- read.csv("stevens_cache_data.csv")		# read in all caching data
cache_data1 <- subset(cache_data, experiment == 1)			# subset out experiment 1 data
cache_data1$date <- as.Date(cache_data1$date)						# convert date column to date format
cache_data1$subject <- factor(cache_data1$subject)			# convert subject column to factor

## Calculate mean caches per subject
# Prepare data frames
cache_stopped1 <- cache_indiff1 <- rep(NA, 4)

# Subset out cache data during indifference point calculation
foreach(curr_subj = subjects) %do% {  # for each subject
  subj_subset1 <- subset(cache_data1, subject == curr_subj)   # subset this subject's data
  cache_end1 <- indiff_means1$end_date[which(indiff_means1$subject == curr_subj)]  # find this subject's stop date for indifference point in the operant task
  cache_stopped1 <- rbind(cache_stopped1, subset(subj_subset1, date <= cache_end1)[, 2:5])		# append this subject's data through the indifference point date
}

# Clean and calculate means
cache_stopped1 <- cache_stopped1[-1, ]	# remove empty data row
cache_means_stopped1 <- aggregate(cbind(cached, eaten) ~ subject, data = cache_stopped1, FUN = "mean")		# calculate each subject's mean number of caches

########
# Combined operant and caching data
########
## Merge data frames and write to file
all_data1 <- merge(indiff_means1, cache_means_stopped1)		# merge operant delay choice and caching data
all_data1$subject <- as.factor(all_data1$subject)					# convert subject column to factor

## Calculate means and 95% CIs for operant delay choice with caching data
operant_summary1 <- ci(all_data1$delay_indiff)  # calculate 95% CIs
operant_mean1 <- operant_summary1$mean     # extract mean
operant_ci1 <- operant_mean1 - operant_summary1$lower95ci  # extract 95% CI
cache_summary1 <- ci(all_data1$cached)  # calculate 95% CIs
cache_mean1 <- cache_summary1$mean     # extract mean
cache_ci1 <- cache_mean1 - cache_summary1$lower95ci  # extract 95% CI
consume_summary1 <- ci(all_data1$eaten)  # calculate 95% CIs
consume_mean1 <- consume_summary1$mean     # extract mean
consume_ci1 <- consume_mean1 - consume_summary1$lower95ci  # extract 95% CI

## Correlate operant delay choice with caching data
# Calculate correlation statistics
indiff_cache_cor1 <- cor.test(all_data1$cached, all_data1$delay_indiff)		# calculate correlation between operant delay and number of caches
indiff_cache_cor_r1 <- round(indiff_cache_cor1$estimate, 2)		# extract correlation coefficient
indiff_cache_cor_p1 <- round(indiff_cache_cor1$p.value, 2)		# extract correlation coefficient

## Plot correlation between operant delay and number of caches
cache_indiff_cor_plot1 <- xyplot(delay_indiff ~ cached, data = all_data1,
  xlab = "Mean number of caches", ylab = "Mean delay to large reward (s)",
  par.settings = list(axis.text = list(cex = 2), par.xlab.text = list(cex = 2.5), par.ylab.text = list(cex = 2.5),
		layout.widths = list(axis.right = 0), layout.heights = list(axis.top = 0)),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, type = c("p", "r"), lwd = 2, cex = 2, ...)
    panel.text(x = 8, y = 25, labels = paste("r =", indiff_cache_cor_r1), cex = 2)		# add correlation coefficient text
  }
)
# Print to PDF
cairo_pdf(file = "figures/expt1_choice_cache_xy.pdf", width = 8, height = 7)
print(cache_indiff_cor_plot1)
dev.off()

# Create xtable of operant and caching task data
expt1_table <- all_data1[, c(1, 4:5)]
names(expt1_table) <- c("Subject", "Indifference point", "Items cached")
expt1_xtable <- xtable(expt1_table, digits = 1, caption = 'Mean Responses for Experiment 1: Correlating Caching and Operant Delay Choice', label = "tab:EXP1MEANS")	# create xtable for import into LaTeX

###################
## Experiment 2: Using photoperiod to manipulate caching and operant delay choice
###################
########
# Operant data
########
## Load and prepare data
operant_data2_all <- subset(operant_data, experiment == 2)  # subset experiment 2 operant delay choice data
operant_data2_all$date <- as.Date(operant_data2_all$date)		# convert date column to date
operant_data2_stable <- subset(operant_data2_all, date >= "2013-11-18")  #  subset data for which the date is after the photoperiod stabilized (2013-11-18)
operant_data2_stable$trial <- c(1, rep(NA, length(operant_data2_stable[, 1]) - 1))  # initialize column of trial numbers

## Create column of total free choice trials for each subject
trial_counter <- 2  # initialize trial counter
for(i in 1:length(operant_data2_stable[, 1])) {  				# for each row of data
  if(i > 1) { # if not the first row
    if(operant_data2_stable$subject[i] == operant_data2_stable$subject[i - 1]) {  # if the subject for this row is the same as the previous row
      operant_data2_stable$trial[i] <- trial_counter   	# write the trial number
      trial_counter <- trial_counter + 1  							# increment the trial counter
    } else {  # if the subject changes between rows
      trial_counter <- 1  															# re-initialize the trial counter for this subject
      operant_data2_stable$trial[i] <- trial_counter   	# write the trial number
      trial_counter <- trial_counter + 1 								# increment the trial counter
    }
  }
}

## Prepare data from first 240 trials
operant_data2 <- subset(operant_data2_stable, trial < 240)		# subset data from first 240 trials
operant_data2$subject <- as.factor(operant_data2$subject)			# convert subject column to factor
operant_data2 <- operant_data2[order(operant_data2$subject, operant_data2$date), ]  # reorder by subject and date

## Calculate mean choice for larger later per subject per session
operant_session_means2 <- aggregate(choiceLL ~ subject * date * condition, data = operant_data2, FUN = "mean", na.rm = TRUE)	# aggregate mean choiceLL over subjects and session
names(operant_session_means2) <- c("subject", "date", "condition", "choiceLL")			# rename columns
operant_end_dates2 <- aggregate(date ~ subject, data = operant_data2, FUN = "max")	# find latest date (end date) for each subject
operant_subj_means2 <- aggregate(choiceLL ~ subject * condition, data = operant_session_means2, FUN = "mean")		# calculate choice means per subject per session
operant_subj_means2$choiceLL <- operant_subj_means2$choiceLL * 100									# convert choice proportions to percentages

## Calculate mean choice for larger later per condition
# Subset data per condition
operant_means_cond_short2 <- subset(operant_subj_means2, condition == "Short day")	# subset subject in short-day photoperiod condition (10:14 h)
operant_means_cond_long2 <- subset(operant_subj_means2, condition == "Long day")		# subset subject in long-day photoperiod condition (14:10 h)

# Calculate mean choice, CIs, and effect sizes per condition
operant_means2 <- aggregate(choiceLL ~ condition, data = operant_subj_means2, FUN = "mean")		# aggregate mean choice per condition
names(operant_means2) <- c("condition", "choiceLL")
operant_means_cond2 <- rbind(ci(subset(operant_subj_means2, condition == "Short day")$choiceLL), ci(subset(operant_subj_means2, condition == "Long day")$choiceLL))		# calculate 95% CIs for each condition
operant_means_cond2$ci <- operant_means_cond2$mean - operant_means_cond2$lower95ci										# calculate the CI range
operant_means_cond2$condition <- as.factor(c("Short day", "Long day"))																					# add condition as factor
operant_means_cond_diff2 <- max(operant_means_cond2$mean) - min(operant_means_cond2$mean)							# calculate the difference between means
operant_means_cond_t2 <- t.test(operant_means_cond_short2$choiceLL, operant_means_cond_long2$choiceLL)	# calculate t-test on mean difference
operant_means_cond_diff_ci2 <- operant_means_cond_diff2 - operant_means_cond_t2$conf.int[1]						# calculate CI for mean difference
operant_means_cond_g2 <- operant_means_cond_diff2 / operant_means_cond2$sd[2]													# calculate Glass's g effect size for mean difference

## Plot choice as a function of photoperiod
operant_light <- stripplot(choiceLL ~ condition, data = operant_subj_means2,
  xlab = "Photoperiod", ylab = "Mean percent choice for larger option", subscripts = TRUE, aspect = 1.4,
  par.settings = list(axis.text = list(cex = 1.6), par.xlab.text = list(cex = 2.2),
    par.ylab.text = list(cex = 2.2)),
  panel = function(x, y, ...) {
    panel.stripplot(x, y, cex = 1.75, ...)
    mean.values <<- tapply(y, x, mean, na.rm=T)          				# calculates mean choice percentage per photoperiod
    panel.points(mean.values, pch = 18, col = "black", cex = 2)	# plot mean choice percentages
  }
)
addCI <- layer_(panel.arrows(condition, lower95ci[subscripts], condition, upper95ci[subscripts], col = "black",  unit = "native", angle = 0, code = 3, lwd = 2), data = operant_means_cond2)  # create layer with within-subject CIs

# Print to PDF
cairo_pdf(file = "figures/expt2_operant_light.pdf", width = 5.75, height = 7)
plot(operant_light + addCI)
dev.off()

########
# Caching data
########
## Load and prepare data
cache_data2_all <- subset(cache_data, experiment == 2)  		# subset experiment 2 caching data
cache_data2_all$date <- as.Date(cache_data2_all$date)				# convert date column to date format
cache_data2_all$subject <- factor(cache_data2_all$subject)	# convert subject column to factor

## Subset caching data within appropriate period (after photoperiod stabilization and before 241st operant trial)
cache_data2 <- cache_data2_all[1, ]		# initiate data frame
foreach(subj = subjects) %do% {				# for each subject
  subj_data <- subset(cache_data2_all, subject == subj & date > "2013-11-18" & date <= operant_end_dates2$date[which(operant_end_dates2$subject == subj)])		# extract data after photoperiod stabiliation (2013-11-18) and up until 240 operant trial
  cache_data2 <- rbind(cache_data2, subj_data)		# append subject data to data frame
}
cache_data2 <- cache_data2[-1, ]		# remove empty row

## Calculate mean number of caches per subject
cache_subj_means2 <- aggregate(cbind(cached, eaten) ~ subject * condition, data = cache_data2, FUN = "mean")	# aggregate mean number of caches per subject
cache_subj_means2$condition <- as.factor(cache_subj_means2$condition)			# convert condition column to factor

## Calculate mean number of caches per subject
# Subset data per condition
cache_subj_means_short2 <- subset(cache_subj_means2, condition == "Short day")	# subset subject in short-day photoperiod condition (10:14 h)
cache_subj_means_long2 <- subset(cache_subj_means2, condition == "Long day")		# subset subject in long-day photoperiod condition (14:10 h)

# Calculate mean caches, CIs, and effect sizes per condition
cache_means_cond2 <- rbind(ci(cache_subj_means_short2$cached), ci(cache_subj_means_long2$cached))		# aggregate mean caches per condition
cache_means_cond2$ci <- cache_means_cond2$mean - cache_means_cond2$lower95ci											# calculate the CI range
cache_means_cond2$condition <- as.factor(c("Short day", "Long day"))																				# add condition as factor
cache_means_cond_diff2 <- max(cache_means_cond2$mean) - min(cache_means_cond2$mean)								# calculate the difference between means
cache_means_cond_t2 <- t.test(cache_subj_means_short2$cached, cache_subj_means_long2$cached, var.equal = FALSE)	# calculate t-test on mean difference
cache_means_cond_diff_ci2 <- cache_means_cond_diff2 - cache_means_cond_t2$conf.int[1]							# calculate CI for mean difference
cache_means_cond_g2 <-  cache_means_cond_diff2 / cache_means_cond2$sd[2]													# calculate Glass's g effect size for mean difference

## Plot choice as a function of photoperiod
cache_light <- stripplot(cached ~ condition, data = cache_subj_means2,
  xlab = "Photoperiod", ylab = "Mean number of caches", ylim = c(-0.5, 18.5),
  subscripts = TRUE, aspect = 1.4,
  par.settings = list(axis.text = list(cex = 1.6), par.xlab.text = list(cex = 2.2),
    par.ylab.text = list(cex = 2.2)),
  panel = function(x, y, ...) {
    panel.stripplot(x, y, cex = 1.75, ...)
    mean.values <<- tapply(y, x, mean, na.rm=T)      		        # calculates mean choice percentage per photoperiod
    panel.points(mean.values, pch = 18, col = "black", cex = 2)	# plot mean choice percentages
  }
)
addCI <- layer_(panel.arrows(condition, lower95ci[subscripts], condition, upper95ci[subscripts], col = "black",  unit = "native", angle = 0, code = 3, lwd = 2), data = cache_means_cond2)  # create layer with within-subject CIs

# Print to PDF
cairo_pdf(file = "figures/expt2_cache_light.pdf", width = 5.75, height = 7)
plot(cache_light + addCI)
dev.off()

# Create xtable of operant and caching task data
expt2_table <- merge(operant_subj_means2, cache_subj_means2)[1:4]
names(expt2_table) <- c("Subject", "Condition", "Percent LL", "Items cached")
expt2_xtable <- xtable(expt2_table, digits = 1, caption = 'Mean Responses for Experiment 2: Manipulating Caching Effects on Operant Delay Choice', label = "tab:EXP2MEANS")	# create xtable for import into LaTeX

## Compare caching in experiments 1 and 2
# Prepare data
cache12 <- merge(all_data1[, c(1, 5)], cache_subj_means2, by = "subject")		# merge caching data from experiments 1 and 2
names(cache12) <- c("subject", "cache_expt1", "condition", "cache_expt2")

# Correlate caching in experiments 1 and 2
cache_consistency_cor <- cor.test(cache12$cache_expt1, cache12$cache_expt2, paired = TRUE)	# correlate caching in experiments 1 and 2
cache_consistency_cor_r <- round(cache_consistency_cor$estimate, 2)  				# extract correlation coefficient
cache_consistency_cor_p <- round(cache_consistency_cor$p.value, 2)   				# extract correlation p-value

## Plot correlation between number of caches across experiments 1 and 2
cache_consistency_cor_plot <- xyplot(cache_expt2 ~ cache_expt1, data = cache12,
  xlab = "Mean caches in Experiment 1", ylab = "Mean caches in Experiment 2",
  par.settings = list(axis.text = list(cex = 2), par.xlab.text = list(cex = 2.5), par.ylab.text = list(cex = 2.5)),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, type = c("p", "r"), lwd = 2, cex = 2, ...)
    panel.text(x = 7.4, y = 17.2, labels = paste("r =", cache_consistency_cor_r), cex = 2)	# add correlation coefficient text
  }
)
# Print to PDF
cairo_pdf(file = "figures/expt1_expt2_cache_xy.pdf", width = 7, height = 7)
print(cache_consistency_cor_plot)
dev.off()

###################
## Experiment 3: Comparing pellet and pine nut caching
###################
## Load and prepare data
cache_data3 <- subset(cache_data, experiment == 3)  # subset experiment 3 caching data
cache_data3$date <- as.Date(cache_data3$date)				# convert date column to date format
cache_data3$subject <- factor(cache_data3$subject)	# convert subject column to factor

## Calculate mean number of caches per subject per condition
mean_pellet_pinenuts_caching <- aggregate(cbind(cached, eaten) ~ subject * condition, data = cache_data3, FUN = "mean")	# calculate mean number of caches per subject per condition
pinenut_caching <- subset(mean_pellet_pinenuts_caching, condition == "pine nuts")	# subset pine nut condition
pellet_caching <- subset(mean_pellet_pinenuts_caching, condition == "pellets")		# subset pellet condition

# Correlate caching pellets and pine nuts
pellet_pinenuts_caching_cor <- cor.test(pinenut_caching$cached, pellet_caching$cached)	# correlate caches of pine nuts and pellets
pellet_pinenuts_caching_cor_r <- sprintf("%.2f", pellet_pinenuts_caching_cor$estimate)	# extract correlation coefficient
pellet_pinenuts_caching_cor_p <- sprintf("%.2f", pellet_pinenuts_caching_cor$p.value)		# extract correlation p-value

## Plot correlation between caching pellets and pine nuts
pellet_pinenut_plot <- xyplot(jitter(pinenut_caching$cached, factor = 1) ~ jitter(pellet_caching$cached, factor = 0.2),
  xlab = "Pellets cached per session", ylab = "Pine nuts cached per session",
  par.settings = list(axis.text = list(cex = 1.8), par.xlab.text = list(cex = 2.5), par.ylab.text = list(cex = 2.5)),
  panel = function(x, y, ...) {
    panel.xyplot(x, y, type = c("p", "r"), lwd = 2, cex = 2, ...)
    panel.text(x = 1.07, y = 22, labels = paste("r =", pellet_pinenuts_caching_cor_r), cex = 2)	# add correlation coefficient text
  }
)
# Print to PDF
cairo_pdf(file = "figures/pellet_pinenut_xy.pdf", width = 7, height = 7)
print(pellet_pinenut_plot)
dev.off()

# Create xtable of pine nut and pellet caching task data
expt3_table <- merge(pinenut_caching, pellet_caching, by = "subject")[, c(1, 3, 6)]
names(expt3_table) <- c("Subject", "Pine nuts cached", "Pellets cached")
expt3_xtable <- xtable(expt3_table, digits = 1, caption = 'Mean Responses for Experiment 3: Caching Pine Nuts vs. Pellets', label = "tab:EXP3MEANS")	# create xtable for import into LaTeX

###################
## Experiment 4: Comparing prefed and not prefed pine nuts
###################
## Load and prepare data
cache_data4 <- subset(cache_data, experiment == 4)  # subset experiment 4 caching data
cache_subj_means4 <- aggregate(cached ~ subject * condition, cache_data4, mean)  # aggregate cached pine nuts by subject and condition (prefed or not)
prefed <- subset(cache_subj_means4, condition == "Pre-fed")        # subset prefed data
notprefed <- subset(cache_subj_means4, condition == "Not pre-fed") # subset nonprefed data

## Calculate means, confidence intervals, and difference
prefed_means_cond4 <- rbind(ci(subset(cache_subj_means4, condition == "Pre-fed")$cached), ci(subset(cache_subj_means4, condition == "Not pre-fed")$cached))		# calculate 95% CIs for each condition
prefed_means_cond4$ci <- prefed_means_cond4$mean - prefed_means_cond4$lower95ci							# calculate the CI range
prefed_means_cond4$condition <- as.factor(c("Pre-fed", "Not pre-fed"))											# add condition as factor
prefed_means_cond_diff4 <- max(prefed_means_cond4$mean) - min(prefed_means_cond4$mean)			# calculate the difference between means
prefed_means_cond_t4 <- t.test(subset(cache_subj_means4, condition == "Not pre-fed")$cached, subset(cache_subj_means4, condition == "Pre-fed")$cached, paired = TRUE)	# calculate t-test on mean difference
prefed_means_cond_diff_ci4 <- prefed_means_cond_diff4 - prefed_means_cond_t4$conf.int[1]		# calculate CI for mean difference
prefed_means_cond_g4 <- prefed_means_cond_diff4 / prefed_means_cond4$sd[2]									# calculate Glass's g effect size for mean difference

# Create xtable of prefed and not prefed caching task data
expt4_table <- merge(prefed, notprefed, by = "subject")[, c(1, 3, 5)]
names(expt4_table) <- c("Subject", "Pre-fed cached", "Not pre-fed cached")
expt4_xtable <- xtable(expt4_table, digits = 1, caption = 'Mean Responses for Experiment 4: Caching When Pre-fed vs. Not Pre-fed', label = "tab:EXP4MEANS")	# create xtable for import into LaTeX

