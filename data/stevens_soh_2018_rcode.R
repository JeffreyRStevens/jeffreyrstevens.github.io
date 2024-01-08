###################################################
### stevens_soh_2018_rcode.R
### Created by Jeffrey R. Stevens on 16 Feb 2015 (jeffrey.r.stevens@gmail.com),
###  finalized on 16 Oct 2017
### Summary: This script investigates which mathematical operations of small
###   and large values predict similarity judgments.
### Instructions: Place this file and the data files in the main directory.
### 	Create a folder called "figures" and create a subfolder within "figures"
###   called "trees". Set the R working directory to the main directory.  
### 	At the R command prompt, type
### 	> source("stevens_soh_2018_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving figures in the figures directory. If packages do not
###		load properly, install them with install.packages("package_name").
### Data files:
###  stevens_soh_2018_data.csv
###   dataset - Number for data set (1 or 2)
###   subject - Participant ID number with data set as first number
###   age - Participant age
###   gender - Participant gender
###   small - Small value
###   large - Large value
###   type - Judgment type (Amount = amount judgment; Delay = delay judgment)
###   similarity - Similarity judgment (0 = dissimilar; 1 = similar)
###   response_time - Response time in seconds from presentation of value pairs to choice
###   difference - Numerical difference: L-S
###   ratio - Numerical ratio: S/L
###   mean_ratio - Mean ratio: S/((S+L)/2)
###   log_ratio - Log ratio: log(S/L)
###   relative_diff - Relative/proportional difference: (L−S)/L
###   disparity - Disparity ratio: (L−S)/((S+L)/2)
###   salience - Salience: (L−S)/(S+L)
###   discrim - Discriminability: log(L/(L−S))
###   logistic - Logistic function: 1/(1+e^(L−S))
### License: This script is released under the Creative Commons 
###   Attribution-NonCommercial-ShareAlike 4.0 International license (CC BY-NC-SA 4.0). 
###   You may share and adapt this content with attribution, for non-commercial purposes 
###   if you ShareAlike (distribute any contributions under the same license).
###################################################

###############################
### Load libraries and clear variables
###############################
suppressPackageStartupMessages(library(BayesFactor))  
# needed for Bayesian statistics (ttestBF)
# Morey, R. D., & Rouder, J. N. (2015). BayesFactor: Computation of Bayes Factors for Common Designs. (R package version 0.9.12-2). Retrieved from https://CRAN.R-project.org/package=BayesFactor
suppressPackageStartupMessages(library(car))          
# needed for Recode
# Fox, J., & Weisberg, S. (2011). An R Companion to Applied Regression (Second edition). Thousand Oaks, CA: Sage. Retrieved from http://socserv.socsci.mcmaster.ca/jfox/Books/Companion
suppressPackageStartupMessages(library(cowplot))          
# needed for plot_grid
# Wilke, C.O. (2016). cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'. R package version 0.8.0.
#https://CRAN.R-project.org/package=cowplot
suppressPackageStartupMessages(library(plyr))        
# needed for data manipulation (ddply, rename)
# Wickham, H. (2011). The split-apply-combine strategy for data analysis. Journal of Statistical Software, 40(1), 1–29.
suppressPackageStartupMessages(library(dplyr))        
# needed for data manipulation (bind_rows, group_by, summarize)
# Wickham, H., & Francois, R. (2016). dplyr: A Grammar of Data Manipulation. (R package version 0.7.4). Retrieved from https://CRAN.R-project.org/package=dplyr
suppressPackageStartupMessages(library(foreach))      
# needed for foreach
# Revolution Analytics, & Weston, S. (2015). foreach: Provides Foreach Looping Construct for R. (R package version 1.4.3). Retrieved from http://CRAN.R-project.org/package=foreach
suppressPackageStartupMessages(library(ggplot2))      
# needed for ggplot
# Wickham, H. (2009). ggplot2: Elegant Graphics for Data Analysis. New York: Springer.
suppressPackageStartupMessages(library(lattice))         
# needed for lattice plots
# Sarkar, D. (2008). Lattice: Multivariate Data Visualization with R. New York: Springer. Retrieved from http://lmdvr.r-forge.r-project.org.
suppressPackageStartupMessages(library(lme4))         
# needed for multilevel modeling (glmer)
# Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software, 67(1), 1–48. https://doi.org/10.18637/jss.v067.i01
suppressPackageStartupMessages(library(MBESS))        
# needed for confidence limits for mean differences (ci.smd)
# Kelley, K. (2017). MBESS: The MBESS R Package. (R package version 4.2.0). Retrieved from https://CRAN.R-project.org/package=MBESS
suppressPackageStartupMessages(library(papaja))        
# needed for formatting R Markdown document in APA style
# Aust, F and Barth, M. (2017). papaja: Create APA manuscripts with R Markdown. (R package version 0.1.0.9492). Retrieved from https://github.com/crsh/papaja
suppressPackageStartupMessages(library(rpart))        
# needed for CART algorithm (rpart)
# Therneau, T., Atkinson, B., & Ripley, B. (2015). rpart: Recursive Partitioning and Regression Trees. (R package version 4.1-11). Retrieved from https://CRAN.R-project.org/package=rpart
suppressPackageStartupMessages(library(rpart.plot))        
# needed to plot rpart trees
# Milborrow, S. (2017). rpart.plot: Plot “rpart” Models: An Enhanced Version of “plot.rpart.” (R package version 2.1.2). Retrieved from https://CRAN.R-project.org/package=rpart.plot
suppressPackageStartupMessages(library(tidyr))        
# needed for tidying up data (gather, spread)
# Wickham, H. (2017). tidyr: Easily Tidy Data with “spread()” and “gather()” Functions. (R package version 0.6.1). Retrieved from https://CRAN.R-project.org/package=tidyr

## Clear all variables
rm(list=ls())         

###############################
### Define functions
###############################
################
# Find items that do not match the vector (the inverse of %in%)
################
"%notin%" <- function(x, table) {  # create function analogous to %in% that searches for items not in a vector
  match(x, table, nomatch = 0) == 0
}

################
# Calculate modified z-score
################
z_score <- function(data, cutoff) {
  MAD <- median(abs(data - median(data)))         # calculate median absolute deviation
  mod_z <- (0.6745 * (data - median(data))) / MAD # calculate modified z score
  data[mod_z > cutoff] <- NA   # replace z scores > cutoff with NA
  return(data)  # return new vector
}

################
# Convert GLMM BIC values to Bayes factor
################
bic_bf <- function(null, alternative) {
  new_bf <- exp((null - alternative) / 2) # convert BICs to Bayes factor
  names(new_bf) <- NULL   # remove BIC label
  return(new_bf)  # return Bayes factor
}

################
# Calculate classification tree
################
tree_model <- function(dataset) {             # runs rpart's CART algorithm for all attributes
  rpart(similarity ~ small + large + difference + ratio + mean_ratio + log_ratio + relative_diff + salience + disparity + discrim + logistic, data = dataset, method = "class", control = rpart.control(maxdepth = 3, minbucket = 2, cp = 0))
}

################
# Calculate fit for given classification tree
################
cart_fit <- function(dataset) {
  cart_tree <<- tree_model(dataset) # run CART algorithm on data set
  cart_predictions <<- predict(cart_tree, type = "class")       # calculate predicted responses
  data_no_nas <- subset(dataset, !is.na(similarity))            # remove cases with NA for similarity judgment
  cart_acc <<- mean(cart_predictions == data_no_nas$similarity) # calculate mean accuracy
}

################
# Conduct regression model
################
regression_model <- function(dataset) {             # runs logistic regression for all attributes
  suppressWarnings(glm(similarity ~ small + large + difference + ratio + mean_ratio + log_ratio + relative_diff + salience + disparity + discrim + logistic, data = dataset, family = binomial(link='logit'), control = list(maxit = 500)))
}

################
# Calculate fit for logistic regression
################
regression_fit <- function(dataset) {
  regression <<- regression_model(dataset)  # run regression model on data set
  regression_predictions <<- predict(regression, type = "response")   # calculate predicted responses
  binary_predictions <<- ifelse(regression_predictions >= 0.5, 1, 0)  # make responses binary
  data_no_nas <- subset(dataset, !is.na(similarity))                  # remove cases with NA as similarity judgment
  reg_acc <<- mean(binary_predictions == data_no_nas$similarity)      # calculate mean accuracy
}

################
# Calculate cross-validated predictive accuracy for tree and regression models
################
crossValidateModels <- function(dataset, N) {
  n <- nrow(dataset)  # find number of rows
  K <- 2              # assign number of folds
  cut <- n %/% K      # find number of rows in first split
  tree_acc <- reg_acc <- numeric(0) # initialize variables
  # Cross validate for N iterations
  for(i in 1:N) { # for each iteration
    random_nums <- runif(n)             # randomly draw n samples from a uniform distribution between 0 and 1
    ranking <- rank(random_nums)        # rank samples
    block <- (ranking - 1) %/% cut + 1  # assign rows to folds
    block <- as.factor(block)           # convert to factor
    tree.acc <- reg.acc <- numeric(0)   # initialize variables
    # Apply to each fold
    for(k in 1:K) { # for each fold
      training_data <- dataset[block != k, ]  # assign training sample
      test_data <- dataset[block == k, ]      # assign test sample
      # Calculate predictive accuracy for tree
      if(sum(as.numeric(as.character(training_data$similarity)), na.rm = TRUE) > 0) { # if at least one judgment is 'similar'
        tree <<- tree_model(training_data)  # run CART algorithm on training data
        tree_predictions <- predict(tree, newdata = test_data, type = "class")  # calculate predicted responses
        acc <- mean(test_data$similarity == tree_predictions, na.rm = TRUE)     # compare predicted and actual responses
      } else { # if no judgments are 'similar'
        acc <- NA   # assign accuracy as NA
      }
      tree.acc <- c(tree.acc, acc)  # append tree accuracy to vector of accuracies
      # Calculate predictive accuracy for regression
      regression <- regression_model(training_data)   # run regression on training data
      regression_predictions <- suppressWarnings(predict(regression, newdata = test_data, type = "response")) # calculate predicted responses
      binary_predictions <<- ifelse(regression_predictions >= 0.5, 1, 0)    # make responses binary
      test_data$reg_acc <- test_data$similarity == binary_predictions       # compare predicted and actual responses
      acc <- mean(test_data$similarity == binary_predictions, na.rm = TRUE) # calculate mean predictive accuracy
      reg.acc <- c(reg.acc, acc)  # append regression accuracy to vector of accuracies
    }
    tree_acc <- c(tree_acc, mean(tree.acc, na.rm = TRUE))   # append tree accuracy to vector of accuracies
    reg_acc <- c(reg_acc, mean(reg.acc, na.rm = TRUE))      # append regression accuracy to vector of accuracies
  }
  tree_accuracy <<- mean(tree_acc, na.rm = TRUE)      # calculate mean predictive accuracy for trees
  regression_accuracy <<- mean(reg_acc, na.rm = TRUE) # calculate mean predictive accuracy for regression
}

###############################
### Define plot theme
###############################
## General plots
plot_theme <- trellis.par.get()
plot_theme$box.dot$pch <- "|"             # define boxplot median symbol
plot_theme$box.rectangle$col <- "black"   # define boxplot rectangle color
plot_theme$box.rectangle$lwd <- 3         # define boxplot rectangle line width
plot_theme$box.rectangle$fill <- "grey90" # define boxplot rectangle fill color
plot_theme$box.umbrella$lty <- 1          # define boxplot whisker line type
plot_theme$box.umbrella$lwd <- 2          # define boxplot whisker line width
plot_theme$box.umbrella$col <- "black"    # define boxplot whisker color
plot_theme$par.xlab.text$cex <- 4         # define x-axis label font size
plot_theme$par.ylab.text$cex <- 4         # define y-axis label font size
plot_theme$axis.text$cex <- 3             # define axis value lable font size
plot_theme$layout.heights$strip <- 3.5    # define strip height
plot_theme$strip.background$col <- c("grey80", "grey90")         # define strip background color

## Plots with two strips
plot_theme2strips <- plot_theme   # copy plot theme
plot_theme2strips$layout.heights$strip <- 2 # define strip height


###############################
### Load and prepare data and calculate descriptive statistics
###############################
## Load and prepare data
data_all <- read.csv("stevens_soh_2018_data.csv")  # read in all data
data_all <- subset(data_all, small != large & (large - small < 80))   # remove attention check questions with identical or very different small/large values
data1 <- subset(data_all, dataset == 1)         # ABC web panel
data2 <- subset(data_all, dataset == 2)         # UNL students--collected for this study
data2$response_time <- z_score(data2$response_time, 3)  # replace outliers with Z scores > 3 with NA
num_subjects1 <- length(unique(data1$subject))  # calculate number of subjects for data set 1
num_subjects2 <- length(unique(data2$subject))  # calculate number of subjects for data set 2
num_subjects_all <- length(unique(data_all$subject))    # calculate number of subjects for all data

## Descriptive stats for amount judgments for data set 1
amount_data1 <- subset(data1, type == "Amount") # subset amount judgments for data set 1
amount_data1_subject <- aggregate(similarity ~ subject, data = amount_data1, FUN = "mean", na.rm = TRUE)  # calculate mean proportion of value pairs judged as similar for each subject
amount_data1_mean <- mean(amount_data1_subject$similarity)  #  calculate overall mean proportion of value pairs judged as similar
amount_data1_sd <- sd(amount_data1_subject$similarity)      #  calculate overall standard deviatin of proportion of value pairs judged as similar
amount_data1_min <- min(amount_data1_subject$similarity)    #  calculate overall minimum proportion of value pairs judged as similar
amount_data1_max <- max(amount_data1_subject$similarity)    #  calculate overall maximum proportion of value pairs judged as similar

## Descriptive stats for amount judgments for data set 2
amount_data2 <- subset(data2, type == "Amount") # subset amount judgments for data set 2
amount_data2_subject <- aggregate(similarity ~ subject, data = amount_data2, FUN = "mean", na.rm = TRUE)# calculate mean proportion of value pairs judged as similar for each subject
amount_data2_mean <- mean(amount_data2_subject$similarity)  #  calculate overall mean proportion of value pairs judged as similar
amount_data2_sd <- sd(amount_data2_subject$similarity)      #  calculate overall standard deviatin of proportion of value pairs judged as similar
amount_data2_min <- min(amount_data2_subject$similarity)    #  calculate overall minimum proportion of value pairs judged as similar
amount_data2_max <- max(amount_data2_subject$similarity)    #  calculate overall maximum proportion of value pairs judged as similar

## Descriptive stats for delay judgments for data set 1
delay_data1 <- subset(data1, type == "Delay")   # subset delay judgments for data set 1
delay_data1_subject <- aggregate(similarity ~ subject, data = delay_data1, FUN = "mean", na.rm = TRUE)# calculate mean proportion of value pairs judged as similar for each subject
delay_data1_mean <- mean(delay_data1_subject$similarity)    #  calculate overall mean proportion of value pairs judged as similar
delay_data1_sd <- sd(delay_data1_subject$similarity)        #  calculate overall standard deviatin of proportion of value pairs judged as similar
delay_data1_min <- min(delay_data1_subject$similarity)      #  calculate overall minimum proportion of value pairs judged as similar
delay_data1_max <- max(delay_data1_subject$similarity)      #  calculate overall maximum proportion of value pairs judged as similar

## Descriptive stats for delay judgments for data set 
delay_data2 <- subset(data2, type == "Delay")   # subset delay judgments for data set 2
delay_data2_subject <- aggregate(similarity ~ subject, data = delay_data2, FUN = "mean", na.rm = TRUE)# calculate mean proportion of value pairs judged as similar for each subject
delay_data2_mean <- mean(delay_data2_subject$similarity)    #  calculate overall mean proportion of value pairs judged as similar
delay_data2_sd <- sd(delay_data2_subject$similarity)        #  calculate overall standard deviatin of proportion of value pairs judged as similar
delay_data2_min <- min(delay_data2_subject$similarity)      #  calculate overall minimum proportion of value pairs judged as similar
delay_data2_max <- max(delay_data2_subject$similarity)      #  calculate overall maximum proportion of value pairs judged as similar

## For tree analysis, remove subjects with inconsistent and/or extreme choice
data_tree1 <- subset(data1, subject %notin% c("1-12", "1-17", "1-27", "1-28", "1-33", "1-34", "1-35", "1-36", "1-39", "1-43", "1-49", "1-54", "1-57", "1-60", "1-63"))  # remove inconsistent/bad subjects from data set 1
data_tree2 <- subset(data2, subject %notin% c("2-7", "2-19", "2-30", "2-31", "2-106", "2-111", "2-112", "2-118", "2-124", "2-126", "2-129", "2-136", "2-139", "2-140", "2-142", "2-150"))      # remove inconsistent/bad subjects from data set 1
data_all_tree <- bind_rows(data_tree1, data_tree2)        # combine data sets 1 and 2
num_subjects_tree1 <- length(unique(data_tree1$subject))  # calculate number of subjects with good data for data set 1
num_subjects_tree2 <- length(unique(data_tree2$subject))  # calculate number of subjects with good data for data set 2
subjects <- unique(data_all_tree$subject)                 # find all subject numbers
num_subjects <- length(subjects)                          # calculate number of subjects with good data for all data

###############################
### Plot example classification plots
###############################
## Participant with difference as only predictor
# Prepare data
example_subject_diff <- subset(data2, subject == "2-3" & type == "Delay")
example_subject_diff$letters <- Recode(example_subject_diff$similarity, "0='D';1='S'")
example_subject_diff_rpart <- tree_model(example_subject_diff)

# Plot similarity by difference and ratio
example_subject_diff_plot <- ggplot(example_subject_diff, aes(x = ratio, y = difference)) +
  geom_text(aes(label = letters, color = letters), size = 12) +  # add letters for similar/dissimilar
  scale_color_manual(values = c("#0072B2", "#009E73")) +  # set color of letters
  geom_hline(yintercept = 5.5) +   # add horizontal line
  labs(x = "Ratio", y = "Difference") +   # label axes
  theme_classic() +   # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), legend.position = "none")  # adjust font sizes

## Participant with ratio as only predictor
# Prepare data
example_subject_ratio <- subset(data2, subject == "2-6" & type == "Delay")
example_subject_ratio$letters <- Recode(example_subject_ratio$similarity, "0='D';1='S'")
example_subject_ratio_rpart <- tree_model(example_subject_ratio)

# Plot similarity by difference and ratio
example_subject_ratio_plot <- ggplot(example_subject_ratio, aes(x = ratio, y = difference)) +
  geom_text(aes(label = letters, color = letters), size = 12) +  # add letters for similar/dissimilar
  scale_color_manual(values = c("#0072B2", "#009E73")) +  # set color of letters
  geom_vline(xintercept = 0.45) +   # add horizontal line
  labs(x = "Ratio", y = "Difference") +   # label axes
  theme_classic() +   # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), legend.position = "none")  # adjust font sizes

## Participant with difference and ratio as predictors
# Prepare data
example_subject_diff_ratio <- subset(data2, subject == "2-21" & type == "Delay")
example_subject_diff_ratio$letters <- Recode(example_subject_diff_ratio$similarity, "0='D';1='S'")
example_subject_diff_ratio_rpart <- tree_model(example_subject_diff_ratio)

# Plot similarity by difference and ratio
example_subject_diff_ratio_plot <- ggplot(example_subject_diff_ratio, aes(x = ratio, y = difference)) +
  geom_text(aes(label = letters, color = letters), size = 12) +  # add letters for similar/dissimilar
  scale_color_manual(values = c("#0072B2", "#009E73")) +  # set color of letters
  geom_hline(yintercept = 3.5) +   # add horizontal line
  geom_segment(aes(x = 0.71, y = 3.5, xend = 0.71, yend = 10)) +  # add vertical line
  labs(x = "Ratio", y = "Difference") +   # label axes
  theme_classic() +   # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), legend.position = "none")  # adjust font sizes

png(file = "figures/example_subject_plot.png", width = 2500, height = 750) # open device
plot(plot_grid(example_subject_diff_plot, NULL, example_subject_ratio_plot, NULL, example_subject_diff_ratio_plot, ncol = 5, rel_widths = c(1, 0.075, 1, 0.075, 1)))
dev.off() # close device

###############################
### Descriptive analysis of ratio and difference effect on similarity
###############################
## Prepare data
data2$question <- paste(data2$small, data2$large, sep = "_")  # create column for question
sim_ratio_diff_type2 <- data2 %>% group_by(difference, ratio, question, type) %>% summarize(similarity = mean(similarity))  # calculate mean similarity judgments by difference, ratio, and type
sim_ratio_diff_type2 <- subset(sim_ratio_diff_type2, ratio > 0) # subset data with nonzero ratio
sim_ratio_diff_subject <- data2 %>% group_by(difference, ratio, question, type, subject) %>% summarize(similarity = mean(similarity))  # calculate mean similarity judgments by difference, ratio, and type for each subject

## Test main effects and interactions of difference, ratio, and type on similarity judgments
# Calculate GLMMs for full model and for removing each main effect and interaction
sim_ratio_diff_type_glmer <- suppressWarnings(glmer(similarity ~ difference * ratio + type + (1 | subject), data = sim_ratio_diff_subject, family = binomial(link='logit'))) # calculate GLMM for full model
sim_ratio_diff_type_glmer_full <- summary(sim_ratio_diff_type_glmer) # summarize GLMM
sim_ratio_diff_type_glmer_summary <- data.frame(parameter = c("difference", "ratio", "type", "differencexratio"), estimate = sim_ratio_diff_type_glmer@beta[2:5], lci = c(-1.101, 0.511, 0.681, 0.403), uci = c(-0.914, 1.692, 0.970, 0.660))  # create data frame of beta estimates and CIs
sim_ratio_diff_type_glmer_nodiff <- suppressWarnings(summary(glmer(similarity ~ ratio + type + (1 | subject), data = sim_ratio_diff_subject, family = binomial(link='logit'))))  # calculate GLMM for model with ratio and type
sim_ratio_diff_type_glmer_noratio <- suppressWarnings(summary(glmer(similarity ~ difference + type + (1 | subject), data = sim_ratio_diff_subject, family = binomial(link='logit'))))  # calculate GLMM for difference and type
sim_ratio_diff_type_glmer_notype <- suppressWarnings(summary(glmer(similarity ~ difference * ratio + (1 | subject), data = sim_ratio_diff_subject, family = binomial(link='logit'))))  # calculate GLMM for difference and ratio
sim_ratio_diff_type_glmer_nointeraction <- suppressWarnings(summary(glmer(similarity ~ difference + ratio + type + (1 | subject), data = sim_ratio_diff_subject, family = binomial(link='logit'))))  # calculate GLMM for difference, ratio, type, but no interaction

# Convert BICs to Bayes factor between full and null models
sim_ratio_diff_nodiff_bf <- bic_bf(sim_ratio_diff_type_glmer_nodiff$AICtab[2], sim_ratio_diff_type_glmer_full$AICtab[2]) # ratio and type
sim_ratio_diff_noratio_bf <- bic_bf(sim_ratio_diff_type_glmer_noratio$AICtab[2], sim_ratio_diff_type_glmer_full$AICtab[2]) # difference and type
sim_ratio_diff_notype_bf <- bic_bf(sim_ratio_diff_type_glmer_notype$AICtab[2], sim_ratio_diff_type_glmer_full$AICtab[2]) # difference and ratio
sim_ratio_diff_nointeraction_bf <- bic_bf(sim_ratio_diff_type_glmer_nointeraction$AICtab[2], sim_ratio_diff_type_glmer_full$AICtab[2]) # no interation

## Generate similiarity by ratio and difference scatter plot
sim_ratio_diff_plot2 <- xyplot(similarity ~ difference | factor(ratio, labels = c("Ratio = 0.5", "Ratio = 0.667", "Ratio = 0.75", "Ratio = 0.8", "Ratio = 0.9")) * type, data = subset(sim_ratio_diff_type2, ratio %in% c(0.5, 0.6667, 0.75, 0.8, 0.9)),
  cex = 2.5, pch = 16, type = "b", col = "grey40", lwd = 4,
  xlab = "Difference", ylab = "Proportion judged similar", ylim = c(-0.05, 1.05),
  scales = list(x = list(relation = "free"), y = list(alternating = 1)),
  strip = strip.custom(par.strip.text = list(cex = 2)),
  between = list(x = 0.5),
  par.settings = plot_theme2strips
)
png(file = "figures/similarity_diff_ratio2.png", height = 800, width = 1200)  # open plot device
plot(sim_ratio_diff_plot2) # plot figure
dev.off() #  close plot device

###############################
### Similarity effect on decision response time
###############################
## Prepare data
response_time_sim <- data2 %>% group_by(question, type, dataset) %>% summarize(similarity = mean(similarity), response_time = median(response_time, na.rm = TRUE)) # calculate mean similarity judgment for each question, type, and data set

## Generate response time by similiarity scatter plot
response_time_similarity_plot <- ggplot(response_time_sim, aes(x = similarity, y = response_time)) +
  geom_point(size = 7) +
  facet_wrap(~ type) +
  labs(x = "Proportion judged similar", y = "Response time (s)") +
  theme_classic() +
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), strip.text.x = element_text(size = 45), strip.background = element_rect(fill="grey80"))  # adjust label font sizes
png(file = "figures/response_time_similarity.png", height = 650, width = 1100)  # open plot device
plot(response_time_similarity_plot) # plot figure
dev.off() #  close plot device

###############################
### Decision tree and regression analysis of attribute effects on similarity
###############################
#######
# Calculate trees and regressions
#######
## Prepare data
data_all_tree$similarity <- as.factor(data_all_tree$similarity)         # convert similarity to factor
data_all_tree$log_ratio <- Recode(data_all_tree$log_ratio, "'-Inf'=NA") # recode -Inf to NA

iterations <- 100   # set number of iterations
subject_count <- 0  # initialize subject counter
cart_subject <- data.frame(matrix(rep(NA, num_subjects * 20), ncol = 20)) # initiate data frame
names(cart_subject) <- c("subject", "overall_tree_accuracy", "overall_regression_accuracy", "overall_tree_fit", "overall_regression_fit", "overall_attribute", "amount_tree_accuracy", "amount_regression_accuracy", "amount_tree_fit", "amount_regression_fit", "amount_attribute", "amount_value", "amount_nodes", "delay_tree_accuracy", "delay_regression_accuracy", "delay_tree_fit", "delay_regression_fit", "delay_attribute", "delay_value", "delay_nodes")  # rename columns
cart_subject$subject <- subjects  # create column of subjects
samples_per_subject <- aggregate(as.numeric(as.character(similarity)) ~ subject, data = data_all_tree, FUN = sum) # sum number of similarity judgments per subject
names(samples_per_subject) <- c("subject", "sum")   # rename columns
small_sample_subjects <- as.character(subset(samples_per_subject, sum < 15)$subject)  # find subjects with fewer than 15 similarity judgments

## Conduct tree and regression analyses for each subject
foreach(current_subject = subjects) %do% {  # for each subject
  current_subject_data <- subset(data_all_tree, subject == current_subject)   # subset data from current subject
  current_subject_amount <- subset(current_subject_data, type == "Amount")    # subset amount data from current subject
  current_subject_delay <- subset(current_subject_data, type == "Delay")      # subset amount data from current subject
  subject_count <- subject_count + 1  # increment subject counter

  # Cross validate tree and regression models for all data
  crossValidateModels(current_subject_data, iterations)          # run tree analysis on current subject's data
  cart_subject[subject_count, c("overall_tree_accuracy")] <- tree_accuracy  # assign tree results to data frame row
  cart_subject[subject_count, c("overall_regression_accuracy")] <- regression_accuracy  # assign tree results to data

  # Fit tree and regression models for all data
  cart_fit(current_subject_data)   # run CART on all data
  cart_subject[subject_count, c("overall_tree_fit")] <- cart_acc  # assign tree results to data frame row
  if(!is.null(cart_tree$splits)) {  # if a tree can be generated
    cart_subject$overall_attribute[subject_count] <- as.character(cart_tree$frame[1, 1]) # extract the most important attribute
  } else {  # if a tree cannot be generated
    cart_subject$overall_attribute[subject_count] <- NA   # assign NA
  }
  regression_fit(current_subject_data)   # run regression fit on all data
  cart_subject[subject_count, c("overall_regression_fit")] <- reg_acc  # assign regression results to data frame row

  # Fit tree and regression models for amount data
  cart_fit(current_subject_amount)  # run CART on amount data
  cart_amount <- cart_tree   # assign amount data to cart_amount
  if(!is.null(cart_amount$splits)) {  # if a tree can be generated
    crossValidateModels(current_subject_amount, iterations)          # run tree analysis on current subject's data
    cart_subject[subject_count, c("amount_tree_accuracy")] <- tree_accuracy  # assign tree results to data frame row
    cart_subject[subject_count, c("amount_regression_accuracy")] <- regression_accuracy  # assign regression results to data
    cart_fit(current_subject_amount)   # run CART for amount
    cart_subject[subject_count, c("amount_tree_fit")] <- cart_acc  # assign tree results to data frame row
    regression_fit(current_subject_amount)   # run regression for amount
    cart_subject[subject_count, c("amount_regression_fit")] <- reg_acc  # assign regression results to data frame row
    cart_subject$amount_attribute[subject_count] <- as.character(cart_amount$frame[1, 1]) # extract the most important attribute
    cart_subject$amount_value[subject_count] <- cart_amount$splits[[1, 4]]   # calculate number of incorrect instances
    cart_subject$amount_nodes[subject_count] <-  length(which(cart_amount$frame[, 1] != "<leaf>"))  # extract number of nodes
  }

  # Fit tree and regression models for delay data
  if(current_subject %notin% small_sample_subjects) { # if the current subject is not in the list of small-sample subjects
    cart_fit(current_subject_delay)  # run CART on delay data
    cart_delay <- cart_tree   # assign delay data to cart_delay
    if(!is.null(cart_delay$splits)) {  # if a tree can be generated
      crossValidateModels(current_subject_delay, iterations)          # run tree analysis on current subject's data
      cart_subject[subject_count, c("delay_tree_accuracy")] <- tree_accuracy  # assign tree results to data frame row
      cart_subject[subject_count, c("delay_regression_accuracy")] <- regression_accuracy  # assign regression results to data
      cart_fit(current_subject_delay)   # run CART for delay
      cart_subject[subject_count, c("delay_tree_fit")] <- cart_acc  # assign tree results to data frame row
      regression_fit(current_subject_delay)   # run regression for delay
      cart_subject[subject_count, c("delay_regression_fit")] <- reg_acc  # assign regression results to data frame row
      cart_subject$delay_attribute[subject_count] <- as.character(cart_delay$frame[1, 1]) # extract the most important attribute
      cart_subject$delay_value[subject_count] <- cart_delay$splits[[1, 4]]   # calculate number of incorrect instances
      cart_subject$delay_nodes[subject_count] <-  length(which(cart_delay$frame[, 1] != "<leaf>"))  # extract number of nodes
    }
  }
}
accuracy <- colMeans(cart_subject[, c("overall_tree_accuracy", "overall_regression_accuracy", "overall_tree_fit", "overall_regression_fit", "amount_tree_accuracy", "amount_regression_accuracy", "amount_tree_fit", "amount_regression_fit", "delay_tree_accuracy", "delay_regression_accuracy", "delay_tree_fit", "delay_regression_fit")], na.rm = TRUE)   # calculate mean predictive accuracy for amount and delay

#######
# Find best attributes
#######
## Extract best attributes
# All amount data
amount_attributes <- table(cart_subject$amount_attribute) # count number of subjects for each attribute
amount_attributes[5] <- 0               # add element for large
names(amount_attributes)[5] <- "large"  # name large element
amount_attributes_percent <- amount_attributes / sum(amount_attributes) * 100 # calculate percentage for each attribute for amount judgments

# All delay data
delay_attributes <- table(cart_subject$delay_attribute) # count number of subjects for each attribute
delay_attributes_percent <- delay_attributes / sum(delay_attributes) * 100  # calculate percentage for each attribute for delay judgments

# Data set 1 amount data
amount_attributes1 <- table(subset(cart_subject, substr(subject, 1, 1) == 1)$amount_attribute) # count number of subjects for each attribute
amount_attributes1[5] <- 0               # add element for large
names(amount_attributes1)[5] <- "large"  # name large element

# Data set 2 amount data
amount_attributes2 <- table(subset(cart_subject, substr(subject, 1, 1) == 2)$amount_attribute) # count number of subjects for each attribute
amount_attributes2[5] <- 0               # add element for large
names(amount_attributes2)[5] <- "large"  # name large element

# Data set 1 delay data
delay_attributes1 <- table(subset(cart_subject, substr(subject, 1, 1) == 1)$delay_attribute) # count number of subjects for each attribute

# Data set 2 delay data
delay_attributes2 <- table(subset(cart_subject, substr(subject, 1, 1) == 2)$delay_attribute) # count number of subjects for each attribute

# Create data frame of best attribute tallies
best_attributes_table <- data.frame(c(1, 1, 2, 2, "All", "All"), c("Amount", "Delay", "Amount", "Delay", "Amount", "Delay"), c(amount_attributes1['large'], delay_attributes1['large'], amount_attributes2['large'], delay_attributes2['large'], paste(sprintf("%0.1f", amount_attributes_percent['large']), "%", sep = ""), paste(sprintf("%0.1f", delay_attributes_percent['large']), "%", sep = "")), c(amount_attributes1['difference'], delay_attributes1['difference'], amount_attributes2['difference'], delay_attributes2['difference'], paste(sprintf("%0.1f", amount_attributes_percent['difference']), "%", sep = ""), paste(sprintf("%0.1f", delay_attributes_percent['difference']), "%", sep = "")), c(amount_attributes1['ratio'], delay_attributes1['ratio'], amount_attributes2['ratio'], delay_attributes2['ratio'], paste(sprintf("%0.1f", amount_attributes_percent['ratio']), "%", sep = ""), paste(sprintf("%0.1f", delay_attributes_percent['ratio']), "%", sep = "")), c(amount_attributes1['relative_diff'], delay_attributes1['relative_diff'], amount_attributes2['relative_diff'], delay_attributes2['relative_diff'], paste(sprintf("%0.1f", amount_attributes_percent['relative_diff']), "%", sep = ""), paste(sprintf("%0.1f", delay_attributes_percent['relative_diff']), "%", sep = "")), c(amount_attributes1['logistic'], delay_attributes1['logistic'], amount_attributes2['logistic'], delay_attributes2['logistic'], paste(sprintf("%0.1f", amount_attributes_percent['logistic']), "%", sep = ""), paste(sprintf("%0.1f", delay_attributes_percent['logistic']), "%", sep = ""))) # create data frame of best attributes
names(best_attributes_table) <- c("Data Set", "Judgment Type", "Large", "Difference", "Ratio", "Relative Difference", "Logistic") # rename columns

#######
# Calculate fitted and predictive accuracy for tree and regression models
#######
## Prepare data
subject_wide <- cart_subject[, c("subject", "overall_tree_accuracy", "overall_regression_accuracy", "overall_tree_fit", "overall_regression_fit", "amount_tree_accuracy", "amount_regression_accuracy", "amount_tree_fit", "amount_regression_fit", "delay_tree_accuracy", "delay_regression_accuracy", "delay_tree_fit", "delay_regression_fit")] # extra relevant columns
subject_stacked <- subject_wide %>% gather(model_type, accuracy, -subject)  # convert from wide to narrow shape
subject_stacked$accuracy <- subject_stacked$accuracy * 100  # convert to percentage
subject_stacked$model_type <- factor(subject_stacked$model_type)  # convert to factor

## Calculate means and within-subjects 95% CIs
subject_summary <- subject_stacked %>% group_by(model_type) %>% summarize(N = length(accuracy), mean_accuracy = mean(accuracy, na.rm = TRUE), sd = sd(accuracy, na.rm = TRUE)) # calculate N, mean, and sd for each model type
subject_summary$ci <- wsci(subject_stacked, id = "subject", dv = "accuracy", factors = "model_type")$accuracy # calculate within-subject 95% CI
subject_summary$lci <- subject_summary$mean_accuracy - subject_summary$ci  	# add lower CIs
subject_summary$uci <- subject_summary$mean_accuracy + subject_summary$ci		# add upper CIs

# Create data frame of model accuracies
model_accuracy_table <- data.frame(rep(c("Amount", "Delay"), each = 2), rep(c("Regression", "Tree"), times = 2), c(paste(sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "amount_regression_accuracy"), "mean_accuracy"]), " [", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "amount_regression_accuracy"), "lci"]), ", ", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "amount_regression_accuracy"), "uci"]), "]", sep = ""), paste(sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "amount_tree_accuracy"), "mean_accuracy"]), " [", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "amount_tree_accuracy"), "lci"]), ", ", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "amount_tree_accuracy"), "uci"]), "]", sep = ""), paste(sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "delay_regression_accuracy"), "mean_accuracy"]), " [", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "delay_regression_accuracy"), "lci"]), ", ", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "delay_regression_accuracy"), "uci"]), "]", sep = ""), paste(sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "delay_tree_accuracy"), "mean_accuracy"]), " [", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "delay_tree_accuracy"), "lci"]), ", ", sprintf("%0.1f", subject_summary[which(subject_summary$model_type == "delay_tree_accuracy"), "uci"]), "]", sep = "")))  # create data frame of model accuracies
names(model_accuracy_table) <- c("Judgment Type", "Model", "Mean Accuracy") # rename columns

## Model comparison differences in predictive accuracy
# Subset data based on model type (regression or tree)
amount_regression_accuracy <- subset(subject_stacked, model_type == "amount_regression_accuracy")  # amount regression accuracy
amount_tree_accuracy <- subset(subject_stacked, model_type == "amount_tree_accuracy")  # amount tree accuracy
delay_regression_accuracy <- subset(subject_stacked, model_type == "delay_regression_accuracy" & !is.na(accuracy))  # delay regression accuracy
delay_tree_accuracy <- subset(subject_stacked, model_type == "delay_tree_accuracy" & !is.na(accuracy))  # delay tree accuracy

# Calculate mean difference and 95% CI
amount_model_diff <- mean(amount_tree_accuracy$accuracy - amount_regression_accuracy$accuracy)  # calculate difference
amount_model_sd <- max(c(subject_summary$sd[subject_summary$model_type == "amount_tree_accuracy"], subject_summary$sd[subject_summary$model_type == "amount_regression_accuracy"]))  # calculate standard deviation
amount_model_ci95 <- ci(amount_tree_accuracy$accuracy - amount_regression_accuracy$accuracy) # calculate 95% CI of mean difference
delay_model_diff <- mean(delay_tree_accuracy$accuracy - delay_regression_accuracy$accuracy)  # calculate difference
delay_model_sd <- max(c(subject_summary$sd[subject_summary$model_type == "delay_tree_accuracy"], subject_summary$sd[subject_summary$model_type == "delay_regression_accuracy"]))  # calculate standard deviation
delay_model_ci95 <- ci(delay_tree_accuracy$accuracy - delay_regression_accuracy$accuracy) # calculate 95% CI of mean difference

# Calculate Cohen's d
amount_model_d <- amount_model_diff / amount_model_sd   # calcuate Cohen's d
delay_model_d <- delay_model_diff / delay_model_sd   # calcuate Cohen's d

# Calculate Bayesian t-test comparing tree and regression predictive accuracy
amount_model_comparison_bf <- ttestBF(amount_regression_accuracy$accuracy, amount_tree_accuracy$accuracy, paired = TRUE) # amount predicitve accuracy t-test
delay_model_comparison_bf <- ttestBF(delay_regression_accuracy$accuracy, delay_tree_accuracy$accuracy, paired = TRUE) # delay predictive accuracy t-test

###############################
### Calculate effect of number of predicted number of node levels on judgment response time
###############################
## Create table of nodes
number_amount_nodes <- table(subset(cart_subject, substr(subject, 1, 1) == "2")$amount_nodes)  # find number of subjects with 1-4 amount nodes
number_amount_nodes[5] <- 0   # add 0 participants for 5 nodes
number_delay_nodes <- table(subset(cart_subject, substr(subject, 1, 1) == "2")$delay_nodes)  # find number of subjects with 1-4 delaynodes
node_table <- data.frame(c("Amount", "Delay"), rbind(number_amount_nodes, number_delay_nodes), row.names = NULL)
names(node_table) <- c("Judgment Type", "1 Node", "2 Nodes", "3 Nodes", "4 Nodes", "5 Nodes")

## Prepare data
# Amount
amount_node_subjects <- subset(cart_subject, amount_nodes > 1 & substr(subject, 1, 1) == "2")  # find subjects with more than 1 node for amount
amount_node_subjects$subject <- as.character(amount_node_subjects$subject)  # convert to character string
amount_nodes <- subset(data_all_tree, subject %in% unique(amount_node_subjects$subject)) # subset data of subjects with more than 1 node for amount
amount_nodes_subjects <- amount_nodes[1, ]  # create data frame from first row
amount_nodes_subjects$node_level <- amount_nodes_subjects$node_num <- NA  # create columns for node level and number of nodes

# Delay
delay_node_subjects <- subset(cart_subject, delay_nodes > 1 & substr(subject, 1, 1) == "2")  # find subjects with more than 1 node for delay
delay_node_subjects$subject <- as.character(delay_node_subjects$subject)  # convert to character string
delay_nodes <- subset(data_all_tree, subject %in% unique(delay_node_subjects$subject)) # subset data of subjects with more than 1 node for delay
delay_nodes_subjects <- delay_nodes[1, ]  # create data frame from first row
delay_nodes_subjects$node_level <- delay_nodes_subjects$node_num <- NA  # create columns for node level and number of nodes

## Calculate node levels for amount data
foreach(current_subject = amount_node_subjects$subject) %do% { # for each subject
  current_subject_amount <- subset(amount_nodes, subject == current_subject & type == "Amount") # subset current subject's amount data
  cart_fit(current_subject_amount)  # run cart on subject data to get nodes
  current_subject_amount$node_num <- as.numeric(row.names(cart_tree$frame)[cart_tree$where]) # get the predicted node number from the tree for each value pair
  current_subject_amount$node_level <- Recode(current_subject_amount$node_num, '1:3=1;4:7=2;8:hi=3')  # recode the node number to the node level 1, 2, or 3+
  amount_nodes_subjects <- rbind(amount_nodes_subjects, current_subject_amount) # append data to data frame
  # Plot participant trees
  # png(filename = paste("figures/trees/tree_amount_", current_subject, ".png", sep = ""), width = 200, height = 200)  # open plot device
  # prp(cart_tree)  # plot tree
  # dev.off() # close plot device
}
amount_nodes_subjects <- amount_nodes_subjects[-1, ]  # remove initial row
amount_nodes_subjects$node_level <- as.numeric(as.character(amount_nodes_subjects$node_level))  # convert to numeric

## Calculate node levels for delay data
foreach(current_subject = delay_node_subjects$subject) %do% { # for each subject
  current_subject_delay <- subset(delay_nodes, subject == current_subject & type == "Delay") # subset current subject's delay data
  cart_fit(current_subject_delay)  # run cart on subject data to get nodes
  current_subject_delay$node_num <- as.numeric(row.names(cart_tree$frame)[cart_tree$where]) # get the predicted node number from the tree for each value pair
  current_subject_delay$node_level <- Recode(current_subject_delay$node_num, '1:3=1;4:7=2;8:hi=3')  # recode the node number to the node level 1, 2, or 3
  delay_nodes_subjects <- rbind(delay_nodes_subjects, current_subject_delay) # append data to data frame
  # Plot participant trees
  # png(filename = paste("figures/trees/tree_delay_", current_subject, ".png", sep = ""), width = 200, height = 200)  # open plot device
  # prp(cart_tree)  # plot tree
  # dev.off() # close plot device
}
delay_nodes_subjects <- delay_nodes_subjects[-1, ]  # remove initial row
delay_nodes_subjects$node_level <- as.numeric(as.character(delay_nodes_subjects$node_level))  # convert to numeric

## Aggregate and combine data
amount_node_subject <- amount_nodes_subjects %>% group_by(subject, node_level, type) %>% summarize(response_time = median(response_time, na.rm = TRUE))  # calculate median response time per subject, node, and type
amount_node_subject$subject <- factor(amount_node_subject$subject)  # remove extraneous subjects
remove_amount_subjects <- table(amount_node_subject$subject)[table(amount_node_subject$subject) == 1] # find subjects not using fast-and-frugal trees
amount_node_subject <- subset(amount_node_subject, subject %notin% names(remove_amount_subjects)) # remove subjects not using fast-and-frugal trees

delay_node_subject <- delay_nodes_subjects %>% group_by(subject, node_level, type) %>% summarize(response_time = median(response_time, na.rm = TRUE))  # calculate median response time per subject, node, and type
delay_node_subject$subject <- factor(delay_node_subject$subject)  # remove extraneous subjects
remove_delay_subjects <- table(delay_node_subject$subject)[table(delay_node_subject$subject) == 1] # find subjects not using fast-and-frugal trees
delay_node_subject <- subset(delay_node_subject, subject %notin% names(remove_delay_subjects)) # remove subjects not using fast-and-frugal trees

node_data <- rbind(amount_node_subject, delay_node_subject)  # append amount and delay data


## Model effects of node level and judgment type on response time
# Conduct linear models with random effects for subject
node_lmer <- lmer(response_time ~ node_level * type + (1 | subject), node_data) # calculate linear mixed model
node_lmer_ci <- suppressMessages(confint(node_lmer, method = "profile"))  # calculate profile likelihood confidence intervals
node_lmer_summary <- data.frame(parameter = c("node", "type", "node*type"), estimate = node_lmer@beta[2:4], lci = node_lmer_ci[4:6, 1], uci = node_lmer_ci[4:6, 2])  # create data frame of beta estimates and CIs
node_lmer_null <- lmer(response_time ~ (1 | subject), node_data) # calculate null linear mixed model
node_lmer_node <- lmer(response_time ~ node_level + (1 | subject), node_data) # calculate linear mixed model for node level
node_lmer_type <- lmer(response_time ~ type + (1 | subject), node_data) # calculate linear mixed model for type
node_lmer_nointeraction <- lmer(response_time ~ node_level + type + (1 | subject), node_data) # calculate linear mixed model for no interaction

# Convert BICs to Bayes factors between full and null models
node_bf_node <- bic_bf(BIC(node_lmer_null), BIC(node_lmer_node))  # Bayes factor for effect of node
node_bf_type <- bic_bf(BIC(node_lmer_null), BIC(node_lmer_type))  # Bayes factor for effect of type
node_bf_node_type <- bic_bf(BIC(node_lmer_nointeraction), BIC(node_lmer))  # Bayes factor for interaction

## Generate plot for node data
# Prepare data
node_type_summary <- node_data %>% group_by(node_level, type) %>% summarize(N = length(accuracy), response_time = mean(response_time, na.rm = TRUE), sd = sd(accuracy, na.rm = TRUE)) # calculate N, mean, and sd for each node
node_type_summary$ci <- wsci(node_data, id = "subject", dv = "response_time", factors = c("node_level", "type"))$response_time # calculate within-subject ci
node_type_summary$lci <- node_type_summary$response_time - node_type_summary$ci  	# add lower CIs
node_type_summary$uci <- node_type_summary$response_time + node_type_summary$ci		# add upper CIs

# Plot boxplot of response time for node number and judgment type
response_time_ggplot <- ggplot(node_data, aes(x = as.factor(node_level), y = response_time)) +  # plot response time as a function of node level
  geom_boxplot(coef=5) +  # plot boxplot
  facet_wrap(~ type) +  # facet by judgment type
  geom_point(data = node_type_summary, aes(x = node_level, y = response_time), size = 6) +  # plot mean response time
  geom_linerange(data = node_type_summary, aes(ymin = lci, ymax = uci), size = 2) + # plot 95% within-subjects CIs
  labs(x = "Node level", y = "Response time (s)") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title=element_text(size=45), axis.text=element_text(size=30), strip.text.x = element_text(size = 45), strip.background = element_rect(fill="grey80"))  # adjust label font sizes
png(file = "figures/tree_response_time.png", height = 700, width = 1000)  # open plot device
plot(response_time_ggplot) # plot figure
dev.off() #  close plot device

###############################
### Demographic effects on similarity
###############################
## Prepare data
data1_subject <- summarize(group_by(data1, subject, age, gender), similarity = mean(similarity, na.rm = TRUE))  # summarize similarity by subject, age, and gender
data2_subject <- summarize(group_by(data2, subject, age, gender), similarity = mean(similarity)) # summarize similarity by subject, age, and gender

###############################
### Supplementary material
###############################
################
## Descriptive effects of attributes on similarity
################
## Data set 1 aggregations
sim_diff_type1 <- aggregate(similarity ~ type * difference, data = data1, FUN = mean)       # calculate mean similarity by type and difference
sim_ratio_type1 <- aggregate(similarity ~ type * ratio, data = data1, FUN = mean)           # calculate mean similarity by type and ratio
sim_mean_ratio_type1 <- aggregate(similarity ~ type * mean_ratio, data = data1, FUN = mean) # calculate mean similarity by type and mean ratio
sim_log_ratio_type1 <- aggregate(similarity ~ type * log_ratio, data = data1, FUN = mean)   # calculate mean similarity by type and log(ratio)
sim_relative_diff_type1 <- aggregate(similarity ~ type * relative_diff, data = data1, FUN = mean)     # calculate mean similarity by type and relative difference
sim_disparity_type1 <- aggregate(similarity ~ type * disparity, data = data1, FUN = mean)   # calculate mean similarity by type and disparity
sim_salience_type1 <- aggregate(similarity ~ type * salience, data = data1, FUN = mean)     # calculate mean similarity by type and salience
sim_discrim_type1 <- aggregate(similarity ~ type * discrim, data = data1, FUN = mean)       # calculate mean similarity by type and discriminability
sim_logistic_type1 <- aggregate(similarity ~ type * logistic, data = data1, FUN = mean)     # calculate mean similarity by type and logistic

## Data set 2 aggregations
sim_diff_type2 <- aggregate(similarity ~ type * difference, data = data2, FUN = mean)       # calculate mean similarity by type and difference
sim_ratio_type2 <- aggregate(similarity ~ type * ratio, data = data2, FUN = mean)           # calculate mean similarity by type and ratio
sim_mean_ratio_type2 <- aggregate(similarity ~ type * mean_ratio, data = data2, FUN = mean) # calculate mean similarity by type and mean ratio
sim_log_ratio_type2 <- aggregate(similarity ~ type * log_ratio, data = data2, FUN = mean)   # calculate mean similarity by type and log(ratio)
sim_relative_diff_type2 <- aggregate(similarity ~ type * relative_diff, data = data2, FUN = mean)     # calculate mean similarity by type and relative difference
sim_disparity_type2 <- aggregate(similarity ~ type * disparity, data = data2, FUN = mean)   # calculate mean similarity by type and disparity
sim_salience_type2 <- aggregate(similarity ~ type * salience, data = data2, FUN = mean)     # calculate mean similarity by type and salience
sim_discrim_type2 <- aggregate(similarity ~ type * discrim, data = data2, FUN = mean)       # calculate mean similarity by type and discriminability
sim_logistic_type2 <- aggregate(similarity ~ type * logistic, data = data2, FUN = mean)     # calculate mean similarity by type and logistic

## Plot similarity judgments by difference
# Prepare data
sim_diff_type <- dplyr::bind_rows(sim_diff_type1, sim_diff_type2) # combine data sets 1 and 2
sim_diff_type$dataset <- as.factor(c(rep("Data set 1", length(sim_diff_type1$similarity)), rep("Data set 2", length(sim_diff_type2$similarity)))) # add data set column

# Plot similarity as a function of type and difference
sim_diff_type_plot <- xyplot(similarity ~ difference | type, groups = dataset, data = subset(sim_diff_type, difference < 80),
  xlab = "Difference", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.995, 0.88), text = list(levels(sim_diff_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by ratio
# Prepare data
sim_ratio_type <- dplyr::bind_rows(sim_ratio_type1, sim_ratio_type2) # combine data sets 1 and 2
sim_ratio_type$dataset <- as.factor(c(rep("Data set 1", length(sim_ratio_type1$similarity)), rep("Data set 2", length(sim_ratio_type2$similarity)))) # add data set column

# Plot similarity as a function of type and ratio
sim_ratio_type_plot <- xyplot(similarity ~ ratio | type, groups = dataset, data = sim_ratio_type,
  xlab = "Ratio", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", layout = c(2, 1), aspect = 1, lwd = 2,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.005, 0.88), text = list(levels(sim_ratio_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by mean ratio (small / mean)
# Prepare data
sim_mean_ratio_type <- dplyr::bind_rows(sim_mean_ratio_type1, sim_mean_ratio_type2) # combine data sets 1 and 2
sim_mean_ratio_type$dataset <- as.factor(c(rep("Data set 1", length(sim_mean_ratio_type1$similarity)), rep("Data set 2", length(sim_mean_ratio_type2$similarity)))) # add data set column

# Plot similarity as a function of type and mean ratio
sim_mean_ratio_type_plot <- xyplot(similarity ~ mean_ratio | type, groups = dataset, data = subset(sim_mean_ratio_type),
  xlab = "Mean ratio", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.005, 0.88), text = list(levels(sim_mean_ratio_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by log ratio
# Prepare data
sim_log_ratio_type <- dplyr::bind_rows(sim_log_ratio_type1, sim_log_ratio_type2) # combine data sets 1 and 2
sim_log_ratio_type$dataset <- as.factor(c(rep("Data set 1", length(sim_log_ratio_type1$similarity)), rep("Data set 2", length(sim_log_ratio_type2$similarity)))) # add data set column

# Plot similarity as a function of type and log ratio
sim_log_ratio_type_plot <- xyplot(similarity ~ log_ratio | type, groups = dataset, data = sim_log_ratio_type,
  xlab = "Log ratio", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.005, 0.88), text = list(levels(sim_log_ratio_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by relative difference
# Prepare data
sim_relative_diff_type <- dplyr::bind_rows(sim_relative_diff_type1, sim_relative_diff_type2) # combine data sets 1 and 2
sim_relative_diff_type$dataset <- as.factor(c(rep("Data set 1", length(sim_relative_diff_type1$similarity)), rep("Data set 2", length(sim_relative_diff_type2$similarity))))

# Plot similarity as a function of type and relative difference
sim_relative_diff_type_plot <- xyplot(similarity ~ relative_diff | type, groups = dataset, data = sim_relative_diff_type,
  xlab = "Relative difference", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.995, 0.88), text = list(levels(sim_relative_diff_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by disparity ratio
# Prepare data
sim_disparity_type <- dplyr::bind_rows(sim_disparity_type1, sim_disparity_type2) # combine data sets 1 and 2
sim_disparity_type$dataset <- as.factor(c(rep("Data set 1", length(sim_disparity_type1$similarity)), rep("Data set 2", length(sim_disparity_type2$similarity))))

# Plot similarity as a function of type and disparity ratio
sim_disparity_type_plot <- xyplot(similarity ~ disparity | type, groups = dataset, data = sim_disparity_type,
  xlab = "Disparity ratio", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.995, 0.88), text = list(levels(sim_disparity_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by salience
# Prepare data
sim_salience_type <- dplyr::bind_rows(sim_salience_type1, sim_salience_type2) # combine data sets 1 and 2
sim_salience_type$dataset <- as.factor(c(rep("Data set 1", length(sim_salience_type1$similarity)), rep("Data set 2", length(sim_salience_type2$similarity))))

# Plot similarity as a function of type and salience
sim_salience_type_plot <- xyplot(similarity ~ salience | type, groups = dataset, data = sim_salience_type,
  xlab = "Salience", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.995, 0.88), text = list(levels(sim_salience_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by discriminability
# Prepare data
sim_discrim_type <- dplyr::bind_rows(sim_discrim_type1, sim_discrim_type2) # combine data sets 1 and 2
sim_discrim_type$dataset <- as.factor(c(rep("Data set 1", length(sim_discrim_type1$similarity)), rep("Data set 2", length(sim_discrim_type2$similarity))))

# Plot similarity as a function of type and discriminability
sim_discrim_type_plot <- xyplot(similarity ~ discrim | type, groups = dataset, data = sim_discrim_type,
  xlab = "Discriminability", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.68, 0.88), text = list(levels(sim_discrim_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

## Plot similarity judgments by logistic
# Prepare data
sim_logistic_type <- dplyr::bind_rows(sim_logistic_type1, sim_logistic_type2) # combine data sets 1 and 2
sim_logistic_type$dataset <- as.factor(c(rep("Data set 1", length(sim_logistic_type1$similarity)), rep("Data set 2", length(sim_logistic_type2$similarity)))) # add data set column

# Plot similarity as a function of type and logistic
sim_logistic_type_plot <- xyplot(similarity ~ logistic | type, groups = dataset, data = sim_logistic_type,
  xlab = "Logistic", ylab = "Proportion judged similar", cex = 2.5, pch = c(2, 16), col = "grey40", lwd = 2, layout = c(2, 1), aspect = 1,
  scales = list(x = list(relation = "free")),
  strip = strip.custom(factor.levels = c("Amount", "Delay"), par.strip.text = list(cex = 3)),
  par.settings = plot_theme,
  key = list(corner = c(0.995, 0.03), text = list(levels(sim_logistic_type$dataset)), lines = list(col = "grey40", lwd = 4), pch = c(2, 16), type = "p", divide = 1, cex = 3)
)

png(file = "figures/similarity_descriptives.png", height = 2000, width = 3000)  # open plot device
plot(sim_diff_type_plot, split = c(1, 1, 3, 3), more = TRUE)
plot(sim_ratio_type_plot, split = c(2, 1, 3, 3), more = TRUE)
plot(sim_mean_ratio_type_plot, split = c(3, 1, 3, 3), more = TRUE)
plot(sim_log_ratio_type_plot, split = c(1, 2, 3, 3), more = TRUE)
plot(sim_relative_diff_type_plot, split = c(2, 2, 3, 3), more = TRUE)
plot(sim_disparity_type_plot, split = c(3, 2, 3, 3), more = TRUE)
plot(sim_salience_type_plot, split = c(1, 3, 3, 3), more = TRUE)
plot(sim_discrim_type_plot, split = c(2, 3, 3, 3), more = TRUE)
plot(sim_logistic_type_plot, split = c(3, 3, 3, 3), more = TRUE)
dev.off() #  close plot device

###############################
### Prepare tables of stimuli for supplementary materials
###############################
## Experiment 1
# Amount
amount_data1_stimuli_raw <- summarize(group_by(amount_data1, small, large), similarity = round(mean(similarity, na.rm = TRUE), 2))  # calculate mean similarity per question
amount_data1_stimuli <- amount_data1_stimuli_raw[order(amount_data1_stimuli_raw$small, amount_data1_stimuli_raw$large), ] # order by small then large value

# Delay
delay_data1_stimuli_raw <- summarize(group_by(delay_data1, small, large), similarity = round(mean(similarity, na.rm = TRUE), 2))  # calculate mean similarity per question
delay_data1_stimuli <- delay_data1_stimuli_raw[order(delay_data1_stimuli_raw$small, delay_data1_stimuli_raw$large), ] # order by small then large value

## Experiment 2
# Amount
amount_data2_stimuli_raw <- summarize(group_by(amount_data2, small, large), similarity = round(mean(similarity, na.rm = TRUE), 2))  # calculate mean similarity per question
amount_data2_stimuli <- amount_data2_stimuli_raw[order(amount_data2_stimuli_raw$small, amount_data2_stimuli_raw$large), ] # order by small then large value

# Delay
delay_data2_stimuli_raw <- summarize(group_by(delay_data2, small, large), similarity = round(mean(similarity, na.rm = TRUE), 2))  # calculate mean similarity per question
delay_data2_stimuli <- delay_data2_stimuli_raw[order(delay_data2_stimuli_raw$small, delay_data2_stimuli_raw$large), ] # order by small then large value
names(amount_data1_stimuli) <- names(delay_data1_stimuli) <- names(amount_data2_stimuli) <- names(delay_data2_stimuli) <- c("Small Value", "Large Value", "Mean Similarity Rating") # rename columns for all dataframes
