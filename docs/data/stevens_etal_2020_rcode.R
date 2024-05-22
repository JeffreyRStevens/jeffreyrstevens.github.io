## ---
##
## Script name: stevens_etal_2020_rcode.R
##
## Purpose of script: This script investigates which machine learning algorithms best capture similarity judgments.
##
## Authors: Dr. Jeffrey R. Stevens (jeffrey.r.stevens@gmail.com) and Alexis Saltzman  (alexis.saltzman@gmail.com)
##
## Date Created: 2018-09-05
##
## Date Finalized: 2020-06-01
##
## License: This script is released under the Creative Commons 
##   Attribution-NonCommercial-ShareAlike 4.0 International license (CC BY-NC-SA 4.0). 
##   You may share and adapt this content with attribution, for non-commercial purposes 
##   if you ShareAlike (distribute any contributions under the same license).
##
## ---
##
## Notes:
## Instructions: Place this file and the data files in the main directory.
## 	Create a folder called "figures". Set the R working directory to the main directory.  
##  At the R command prompt, type
## 	> source("stevens_etal_2020_rcode.R")
## 	This will run the script, adding all of the calculated variables to the workspace and 
##  saving figures in the figures directory. If packages do not load properly, install them 
##  with install.packages("package_name").
## Data files:
##  stevens_etal_2020_data_raw.csv (raw data file from Stevens & Soh (2018))
##   dataset - Number for data set (1 or 2)
##   subject - Participant ID number with data set as first number
##   age - Participant age
##   gender - Participant gender
##   small - Small value
##   large - Large value
##   type - Judgment type (Amount = amount judgment; Delay = delay judgment)
##   similarity - Similarity judgment (0 = dissimilar; 1 = similar)
##   response_time - Response time in seconds from presentation of value pairs to choice
##   difference - Numerical difference: L-S
##   ratio - Numerical ratio: S/L
##   mean_ratio - Mean ratio: S/((S+L)/2)
##   log_ratio - Log ratio: log(S/L)
##   relative_diff - Relative/proportional difference: (L−S)/L
##   disparity - Disparity ratio: (L−S)/((S+L)/2)
##   salience - Salience: (L−S)/(S+L)
##   discrim - Discriminability: log(L/(L−S))
##   logistic - Logistic function: 1/(1+e^(L−S))
##  stevens_etal_2020_all_data_accuracy.csv (accuracy, precision, and recall data calculated from the script)
##   data_file - Identifier for data type (Amount or Delay) and set (1 or 2)
##   algorithm - Algorithm used
##   sample_size - Number of instances used in training set
##   train_test - Flag for whether results are for training set or testing set
##   subject - Participant ID (includes data set and ID number)
##   accuracy - Mean accuracy
##   precision - Mean precision
##   recall - Mean recall
##   order - Order in which training set was drawn (Random or Sequential)
##  stevens_etal_2020_data_importance.csv (predictor importance data calculated from this script)
##   data_file - Identifier for data type (Amount or Delay) and set (1 or 2)
##   algorithm - Algorithm
##   subject - Participant ID (includes data set and ID number)
##   predictor - Predictor
##   importance - Predictor importance
##
## ---

# Load libraries ---------------------------------------------------------------

# General and plotting packages
library(foreach)      # needed to use foreach function
library(GGally)       # needed for pairwise correlation plots
library(papaja)       # needed for within-subjects confidence intervals
library(patchwork)    # needed to arrange multiple ggplots
library(tidytext)     # needed to reorder independently for each facet
library(tidyverse)    # needed for tidyverse functions
# Machine learning packages
library(C50)          # needed for C5.0 algorithm
library(caret)        # needed for training and testing machine learning algorithms
library(e1071)        # needed for SVM algorithm
library(kernlab)      # needed for algorithm analyses
library(naivebayes)   # needed for Naive Bayes algorithm
library(nnet)         # needed for Neural Network algorithm
library(randomForest) # needed for Random Forest algorithm
library(rpart)        # needed for CART algorithm

# Define functions ----------------------------------------------------------

###
## Use methods from the caret package to train desired algorithm and test predictions
## against the testing data on desired data set with desired % data split
###
train_test_algorithms <- function(training, testing, algorithm, sample_size, data_file, current_subject, rep, output) {
  # Set parameters for each algorithm
  # CART
  if(algorithm == "rpart") {
    tuned <- train(similarity ~., data = training, method = "rpart", parms = list(split = "information"),  tuneLength = 3)  # split the list of parameters for information gain
  }
  # C5.0
  else if(algorithm == "c5.0") { 
    tuned <- train(similarity ~., data = training, method = "C5.0", verbose = FALSE)
  }
  # Naive Bayes
  else if(algorithm == "naive_bayes") {
    options(warn = -1)
    if(sum(as.numeric(as.character(training$similarity))) > 1) {  # if there is enough variability in similarity
      trctrl <- trainControl(method = "boot")  # use bootstrapping
    } else {  # if there is not enough variability
      trctrl <- trainControl(method = "cv", number = 100)  # use cross-validation
    }
    tuned <- train(similarity ~., data = training, method = "naive_bayes", trControl = trctrl)
  }
  # KNN
  else if(algorithm == "knn") {
    tuned <- train(similarity ~., data = training, method = "knn", tuneGrid = expand.grid(k = 5))
  }
  # Random forest
  else if(algorithm == "rf") {
    tuned <- train(similarity ~., data = training, method = "rf", tuneLength = 3, importance = TRUE)
  }
  # Support vector machine
  else if(algorithm == "svm"){
    tuned <- train(similarity ~., data = training, method = "svmLinear", tuneLength = 10)
  }
  # Neural network
  else if(algorithm == "nnet"){
    tuned <- train(similarity ~., data = training, method = "nnet", trace = FALSE)
  }
  # Logistic regression
  else if(algorithm == "glm"){
    tuned <-train(similarity ~., data = training, method = "glm", family = binomial())
  } else {
    stop("Unfamiliar algorithm type.")
  }
  
  # Calculate accuracy, precision, and recall
  if(output == "accuracy") {  # if the output type is accuracy
    # Extract training measures
    train_pred <- predict(tuned, newdata = training)  # use parameters to predict training set
    train_cm <- confusionMatrix(train_pred, training$similarity)  # calculate confusion matrix
    train_true_pos <- train_cm$table[2, 2]  # extract true positive
    train_true_neg <- train_cm$table[1, 1]  # extract true negative
    train_false_pos <- train_cm$table[2, 1]  # extract false positive
    train_false_neg <- train_cm$table[1, 2]  # extract false negative
    train_accuracy <- train_cm$overall['Accuracy']  # extract accuracy
    train_result_df <- c(data_file, algorithm, sample_size, "training", rep, current_subject, train_true_pos, train_true_neg, train_false_pos, train_false_neg, train_accuracy)  # create vector of results
    names(train_result_df) <- acc_column_names  # assign column names
    
    # Extract testing measures
    test_pred <- predict(tuned, newdata = testing)  # use parameters to predict testing set
    test_cm <- confusionMatrix(test_pred, testing$similarity)  # calculate confusion matrix
    test_true_pos <- test_cm$table[2, 2]  # extract true positive
    test_true_neg <- test_cm$table[1, 1]  # extract true negative
    test_false_pos <- test_cm$table[2, 1]  # extract false positive
    test_false_neg <- test_cm$table[1, 2]  # extract false negative
    test_accuracy <- test_cm$overall['Accuracy']  # extract accuracy
    test_result_df <- c(data_file, algorithm, sample_size, "testing", rep, current_subject, test_true_pos, test_true_neg, test_false_pos, test_false_neg, test_accuracy)  # create vector of results
    names(test_result_df) <- acc_column_names  # assign column names
    acc_result_df <- bind_rows(train_result_df, test_result_df)  # bind training and testing results
    return(acc_result_df)
  } 
  # Calculate predictor importance
  else if(output == "importance") {  # if output type is importance
    att_imp <- varImp(tuned)   # calculate variable importance
    data.frame("data_file" = data_file, "algorithm" = algorithm, "subject" = current_subject, "predictor" = row.names(att_imp$importance), "importance" = att_imp$importance[, 1])
  } else {
    stop("Unfamiliar output type: Must be 'accuracy' or 'importance'.")
  }
} 

###
## Complete the training and testing of each algorithm over each data split for each subject
###
algorithm_accuracy <- function(data_file, data_split, order, reps){
  data_df <- eval(parse(text = data_file))  # import data file
  subjects <- unique(data_df$subject)  # find all subject numbers
  algorithms <- c("rpart", "c5.0", "naive_bayes", "knn", "rf", "svm", "nnet", "glm")  # create vector of 
  foreach(current_subject = subjects)  %do% {  # for each subject
    subject_data <- filter(data_df, subject == current_subject) %>%   # extract subject's data
      select(-subject) %>%  # remove subject column
      mutate(similarity = as.factor(similarity))  # convert to factor
    print(current_subject)  # print subject's number
    foreach(sample_size = data_split) %do% {  # for each sample size
      # print(sample_size)
      foreach(rep = 1:reps) %do% {  # for each replicate
        set.seed(as.numeric(Sys.time()))  # set the random seed
        if(order == "sequential") {  # if analysis should be on sequential data
          training <- subject_data[1:sample_size, ]  # assign selected rows to training data frame
          shuffle_numbers <- sample((sample_size + 1):nrow(subject_data))  # shuffle the row numbers
          shuffled_rows <- subject_data[shuffle_numbers, ]  # shuffle the rows
          testing <- shuffled_rows[1:10, ]  # assign next 10 rows as testing data frame
          if(length(unique(training$similarity)) > 1) {  # if the training sample has at least one similar and one dissimilar judgment
            foreach(algorithm = algorithms) %do% {  # for each algorithm
              # print(algorithm)
              alg_acc <- train_test_algorithms(training, testing, algorithm, sample_size, data_file, current_subject, rep, "accuracy")  # build training results dataframe
              alg_acc_all <- bind_rows(alg_acc_all, alg_acc)  # bind current results to previous results
            }
          }
        } else if(order == "random") {  # if analysis should be on randomized data
          training_rows <- as.vector(createDataPartition(subject_data$similarity, p = sample_size / nrow(subject_data), list = FALSE))
          training <- subject_data[training_rows, ]  # assign selected rows to training data frame
          testing_all <- subject_data[-training_rows, ]  # assign unselected rows to testing data frame
          testing_all <- testing_all[sample(1:nrow(testing_all)), ]  # shuffle row order
          testing <- testing_all[1:10, ]  # assign next 10 rows as testing data frame
          if(length(unique(training$similarity)) > 1) {  # if the training sample has at least one similar and one dissimilar judgment
            foreach(algorithm = algorithms) %do% {  # for each algorithm
              alg_acc <- train_test_algorithms(training, testing, algorithm, sample_size, data_file, current_subject, rep, "accuracy")  # build training results dataframe
              alg_acc_all <- bind_rows(alg_acc_all, alg_acc)  # bind current results to previous results
            }
          }
        }
      }
      alg_acc_summarized <- alg_acc_all[apply(alg_acc_all, 1, function(x)any(!is.na(x))), ] %>%  # remove rows of all NAs
        group_by(data_file, algorithm, sample_size, train_test, subject) %>%   # group by everything but replicate
        summarise(mean_true_positive = mean(as.numeric(true_pos), na.rm = TRUE),  # calculate mean true positives over replicate
                  mean_true_negative = mean(as.numeric(true_neg), na.rm = TRUE),  # calculate mean true negatives over replicate
                  mean_false_positive = mean(as.numeric(false_pos), na.rm = TRUE),  # calculate mean false positives over replicate
                  mean_false_negative = mean(as.numeric(false_neg), na.rm = TRUE),  # calculate mean false negatives over replicate
                  mean_accuracy = mean(as.numeric(accuracy), na.rm = TRUE)  # calculate mean accuracy over replicate
        )
    }
  }
  return(alg_acc_summarized)  # return results
}

###
## Calculate predictor importance
###
predictor_importance <- function(data_file){
  data_df <- eval(parse(text = data_file))  # import data file
  subjects <- unique(data_df$subject)  # find all subject numbers
  algorithms = c("rpart", "c5.0", "naive_bayes", "knn", "rf", "glm", "nnet") #svm left out because it does not provide variable importance data
  foreach(current_subject = subjects)  %do% {  # for each subject
    att_imp_all <- data.frame(matrix(rep(NA, 5), nrow = 1))   # initialize accuracy data frame
    names(att_imp_all) <- imp_column_names  # rename data frame columns
    subject_data <- filter(data_df, subject == current_subject) %>%   # extract subject's data
      select(-subject) %>%  # remove subject column
      mutate(similarity = as.factor(similarity))  # convert to factor
    print(current_subject)  # print subject's number
    training <- data.frame(subject_data)
    rep <- 0  # assign number of replicates
    foreach(algorithm = algorithms) %do% {  # for each algorithm
      att_imp <- train_test_algorithms(training, testing, algorithm, sample_size, data_file, current_subject, rep, "importance")  # build training results dataframe
      att_imp_all <- bind_rows(att_imp_all, att_imp)  # bind current results to previous results
    }
    att_imp_all <- att_imp_all[apply(att_imp_all, 1, function(x)any(!is.na(x))), ]  # remove rows of all NAs
    write_csv(att_imp_all, "stevens_etal_2020_data_importance.csv", append = TRUE)
  }
}

###
## Calculate boxplot data without outliers
## From: https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
###
calc_boxplot_stat <- function(x) {
  coef <- 1.5  # use 1.5 x the interquartile range for whiskers
  n <- sum(!is.na(x))  # find numer of non-NA data points
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))  # calculate quantiles
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")  # name statistics
  iqr <- diff(stats[c(2, 4)])  # calculate difference between 25th and 75th quantile
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)  # calculte outliers
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)  # calculate 1.5 times interquartile range without outliers
  }
  return(stats)
}

###
## Create themes for plots
###
# Theme for figures without legends
theme_plots <- function () { 
  theme_bw(base_size = 25, base_family = "Arial") %+replace%   # set font size and family
    theme(
      panel.grid = element_blank(),  # remove grid lines
      legend.position = "none"  # remove legend
    )
}

# Theme for figures with legends
theme_legend <- function () { 
  theme_bw(base_size = 25, base_family = "Arial") %+replace%    # set font size and family
    theme(
      panel.grid = element_blank(),  # remove grid lines
      legend.title = element_blank(),  # remove legend title
      legend.key.width = unit(3, "line"),  # increase width of legend lines
      legend.key = element_rect(fill = "transparent", color = NA),  # make legend background transparent
      legend.background = element_rect(fill = "transparent", color = NA)  # make legend background transparent
    )
}


# Load and prepare data and calculate descriptive statistics ----------------------------------------------------

## Load and prepare data
data_all <- read_csv("stevens_etal_2020_data_raw.csv") %>% # import data
  filter(!is.na(similarity) & subject != "2-147") %>%  # remove NAs and subject 2-147 (not enough variation in choice)
  select(dataset:gender, type, small, large, difference:logistic, similarity)  # rearrange columns
data_all[data_all ==  -Inf] <- -100  # replace -Inf with -10
data_trimmed <- select(data_all, -mean_ratio, -relative_diff, -salience, -disparity)  # remove unused predictors
data1 <- filter(data_all, dataset == 1)         # ABC web panel
data2 <- filter(data_all, dataset == 2)         # UNL students--collected for this study
num_subjects1 <- length(unique(data1$subject))  # calculate number of subjects for data set 1
num_subjects2 <- length(unique(data2$subject))  # calculate number of subjects for data set 2
num_subjects_all <- length(unique(data_all$subject))    # calculate number of subjects for all data

## Summarize demographic data
data1_demo <- data1 %>% 
  filter(small == 3, large == 13, type == "Amount") %>% # get one row per subject
  select(dataset, subject, age, gender) # select demographic columns
data1_gender <- table(data1_demo$gender)  # find gender split
data1_age <- summary(data1_demo$age)   # calculate mean, median, range for age
data1_agesd <- sd(data1_demo$age)   # calculate standard deviation for age
data2_demo <- data2 %>% 
  filter(small == 9, large == 18, type == "Amount") %>% # get one row per subject
  select(dataset, subject, age, gender) # select demographic columns
data2_gender <- table(data2_demo$gender)  # find gender split
data2_age <- summary(data2_demo$age)   # calculate mean, median, range for age
data2_agesd <- sd(data2_demo$age)   # calculate standard deviation for age

## Divide data into data sets and types
data1_trimmed <- filter(data_trimmed, dataset == 1)         # ABC web panel
data2_trimmed <- filter(data_trimmed, dataset == 2)         # UNL students--collected for this study
amount_1 <- data1_trimmed %>% 
  filter(type == "Amount") %>%  # filter amount judgments
  select(-dataset, -age, -gender, -type)
delay_1 <- data1_trimmed %>% 
  filter(type == "Delay") %>%  # filter amount judgments
  select(-dataset, -age, -gender, -type)
amount_2 <- data2_trimmed %>% 
  filter(type == "Amount") %>%  # filter amount judgments
  select(-dataset, -age, -gender, -type)
delay_2 <- data2_trimmed %>% 
  filter(type == "Delay") %>%  # filter amount judgments
  select(-dataset, -age, -gender, -type)

## Remove unnecessary objects
rm(data_trimmed, data1, data2, data1_demo, data2_demo, data1_trimmed, data2_trimmed)  # remove objects
gc()  # use garbage collection to return memory to the system


# Accuracy, precision, and recall ----------------------------------------------------

## Prepare iterated variables
data_files <- c("amount_1", "delay_1", "amount_2", "delay_2")  # create vector of data files
data_splits = c(15, 20, 25, 30)  # create vector of data splits

## Algorithm accuracy
# This component of the analysis takes a very long time to run. We are including the code to run in sequence, but we actually ran this in parallel on a super computer. To speed further analysis, we include the file generated by this analysis (stevens_etal_data_accuracy.csv).
# Prepare accuracy data frame
# acc_column_names <- c("data_file", "algorithm", "sample_size", "train_test", "rep", "subject", "true_pos", "true_neg", "false_pos", "false_neg", "accuracy")  # create vector of column names
# alg_acc_all <- data.frame(matrix(rep(NA, 11), nrow = 1))   # initialize accuracy data frame
# names(alg_acc_all) <- acc_column_names  # rename data frame columns
# 
# 
# # Random
# accuracy_random <- algorithm_accuracy(data_files, data_splits, "random", 100)  # calculate accuracy, precision, and recall
# accuracy_random <- accuracy_random %>% 
#   group_by(data_file, algorithm, sample_size, train_test, subject) %>%   # for each data file, algorithm, sample size, train_test, and subject
#   summarise(true_positive = mean(mean_true_positive, na.rm = TRUE),  # calculate mean true positives
#             true_negative = mean(mean_true_negative, na.rm = TRUE),  # calculate mean true negatives
#             false_positive = mean(mean_false_positive, na.rm = TRUE),  # calculate mean false positives
#             false_negative = mean(mean_false_negative, na.rm = TRUE),  # calculate mean false negatives
#             accuracy = mean(mean_accuracy, na.rm = TRUE),  # calculate mean accuracy
#             order = "Random")  # create order column
# 
# # Sequential
# accuracy_sequential <- algorithm_accuracy(data_files, data_splits, "sequential", 100)  # calculate accuracy, precision, and recall
# accuracy_sequential <- accuracy_sequential %>% 
#   group_by(data_file, algorithm, sample_size, train_test, subject) %>%   # for each data file, algorithm, sample size, train_test, and subject
#   summarise(true_positive = mean(mean_true_positive, na.rm = TRUE),  # calculate mean true positives
#             true_negative = mean(mean_true_negative, na.rm = TRUE),  # calculate mean true negatives
#             false_positive = mean(mean_false_positive, na.rm = TRUE),  # calculate mean false positives
#             false_negative = mean(mean_false_negative, na.rm = TRUE),  # calculate mean false negatives
#             accuracy = mean(mean_accuracy, na.rm = TRUE),  # calculate mean accuracy
#             order = "Sequential")  # create order column
# 
# # Combine and write data
# accuracy_all <- bind_rows(accuracy_random, accuracy_sequential) %>%   # combine data frames
#   select(order, everything())  # place order column first
# write_csv(accuracy_all, "stevens_etal_2020_data_accuracy.csv")  # write data to file



# This component of the analysis takes a while to run. To speed further analysis, we include the file generated by this analysis (stevens_etal_data_importance.csv).
# Predictor importance ------------------------------------------

## Predictor importance
# Prepare importance data frame
# imp_column_names <- c("data_file", "algorithm", "subject", "predictor", "importance")  # create vector of column names
# write_csv(data.frame(matrix(imp_column_names, nrow = 1)), "stevens_etal_2020_data_importance.csv", col_names = FALSE)  # initiate data file
# 
# # Calculate predictor importance for random samples from all four data sets
# importance <- foreach(dataset = data_files, .combine = rbind) %do% {  # for each data set
#   predictor_importance(dataset)   # calculate predictor importance
# }



# Plots -------------------------------------------------------------------

# _Pairwise predictor correlations -------
predictor_amount_data <- data_all %>%
  filter(type == "Amount") %>%  # filter amount data
  unite("question_type", small:large, remove = FALSE) %>%  # create question_type column
  group_by(question_type) %>%  # for each question type
  summarize_at(vars(c(small, large, difference:logistic)), mean, na.rm = TRUE)  # calculate mean values

predictor_delay_data <- data_all %>%
  filter(type == "Delay" & log_ratio != -100) %>%  # filter amount data
  unite("question_type", small:large, remove = FALSE) %>%  # create question_type column
  group_by(question_type) %>%  # for each question type
  summarize_at(vars(c(small, large, difference:logistic)), mean, na.rm = TRUE)  # calculate mean values

# ggpairs(predictor_amount_data, columns = 2:12, columnLabels = c("Small", "Large", "Difference", "Ratio", "Mean\nratio", "Log ratio", "Relative\ndifference", "Salience", "Disparity", "Discrim.", "Logistic"), axisLabels = "none", upper = list(continuous = ggally_cor_v1_5))
# ggsave(filename = "figures/predictor_correlations_amount.png", width = 9, height = 9)
# 
# ggpairs(predictor_delay_data, columns = 2:12, columnLabels = c("Small", "Large", "Difference", "Ratio", "Mean\nratio", "Log ratio", "Relative\ndifference", "Salience", "Disparity", "Discrim.", "Logistic"), axisLabels = "none", upper = list(continuous = ggally_cor_v1_5))
# ggsave(filename = "figures/predictor_correlations_delay.png", width = 9, height = 9)

## _Accuracy, precision, and recall -------

# Prepare all accuracy, precision, and recall data
accuracy_all <- read_csv("stevens_etal_2020_data_accuracy.csv") %>%   # import accuracy data
  separate(data_file, into = c("type", "dataset")) %>%  # create separate columns for type and dataset
  mutate(type = fct_recode(type, "Amount" = "amount", "Delay" = "delay"),  # rename type
         dataset = fct_recode(dataset, "Data set 1" = "1", "Data set 2" = "2"),  # rename dataset
         algorithm = fct_recode(algorithm, "C5.0" = "c5.0","kNN" = "knn", "Regression" = "glm", "Naive Bayes" = "naive_bayes", "Neural Network" = "nnet", "Random Forest" = "rf", "CART" = "rpart", "SVM" = "svm"),  # rename algorithms
         algorithm = fct_relevel(algorithm, "kNN", "Regression", "CART", "Naive Bayes", "C5.0", "Neural Network", "Random Forest", "SVM"),  # rename algorithms
         precision = true_positive / (true_positive + false_positive), # calculate precision
         recall = true_positive / (true_positive + false_negative)  # calculate recall
  )

# Split into random and sequential ordering
accuracy_random_all <- filter(accuracy_all, order == "Random")  # filter random order data
accuracy_sequential_all <- filter(accuracy_all, order == "Sequential")  # filter sequential order data

# Create long data
## Training
all_training <- accuracy_random_all %>%
  filter(train_test == "training"  & sample_size == 30) %>%  # select testing data of largest sample size
  pivot_longer(true_positive:recall, names_to = "measure", values_to = "value") %>%  # pivot to long data
  mutate(measure = fct_recode(measure, "True positive" = "true_positive", "True negative" = "true_negative", "False positive" = "false_positive", "False negative" = "false_negative", "Accuracy" = "accuracy", "Precision" = "precision", "Recall" = "recall"))  # recode factor levels
## Testing
all_testing <- accuracy_random_all %>%  # testing
  filter(train_test == "testing"  & sample_size == 30) %>%  # select testing data of largest sample size
  pivot_longer(true_positive:recall, names_to = "measure", values_to = "value") %>%  # pivot to long data
  mutate(measure = fct_recode(measure, "True positive" = "true_positive", "True negative" = "true_negative", "False positive" = "false_positive", "False negative" = "false_negative", "Accuracy" = "accuracy", "Precision" = "precision", "Recall" = "recall"))  # recode factor levels

# Overall
## Summarize by subject, algorithm, and measure
acc_prec_recall <- c("Accuracy", "Precision", "Recall")
## Training
all_subject_training <- all_training %>%
  group_by(subject, algorithm, measure) %>%  # for each subject, algorithm, and measure
  summarise(mean_value = mean(value, na.rm = TRUE))  # calculate mean value
acc_prec_recall_subject_training <- filter(all_subject_training, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall
## Testing
all_subject_testing <- all_testing %>%
  group_by(subject, algorithm, measure) %>%  # for each subject, algorithm, and measure
  summarise(mean_value = mean(value, na.rm = TRUE))  # calculate mean value
acc_prec_recall_subject_testing <- filter(all_subject_testing, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Add within-subjects 95% confidence intervals
## Training
all_wsci_training <- wsci(all_subject_training, id = "subject", factors = c("algorithm", "measure"), dv = "mean_value")  # calculate 95% WSCI
all_means_training <- all_subject_training %>%
  group_by(algorithm, measure) %>%  # for each algorithm and measure
  summarise(value = mean(mean_value, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(all_wsci_training, by = c("algorithm", "measure")) %>%  # join the 95% WSCI data
  rename(ci = mean_value) %>%  # rename CI column
  mutate(uci = value + ci,  # create upper CI column
         lci = value - ci)  # create lower CI column
acc_prec_recall_means_training <- filter(all_means_training, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall
## Testing
all_wsci_testing <- wsci(all_subject_testing, id = "subject", factors = c("algorithm", "measure"), dv = "mean_value")  # calculate 95% WSCI
all_means_testing <- all_subject_testing %>%
  group_by(algorithm, measure) %>%  # for each algorithm and measure
  summarise(value = mean(mean_value, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(all_wsci_testing, by = c("algorithm", "measure")) %>%  # join the 95% WSCI data
  rename(ci = mean_value) %>%  # rename CI column
  mutate(uci = value + ci,  # create upper CI column
         lci = value - ci)  # create lower CI column
acc_prec_recall_means_testing <- filter(all_means_testing, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Plot
## Training
ggplot(acc_prec_recall_subject_training, aes(x = reorder_within(algorithm, by = mean_value, within = measure, fun = "mean"), y = mean_value, fill = algorithm)) +
  facet_wrap(~measure, scale = "free_x") +  # create panel for each measure
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots without outliers
  geom_point(data = acc_prec_recall_means_training, aes(x = reorder_within(algorithm, by = value, within = measure, fun = "mean"), y = value), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = acc_prec_recall_means_training, aes(x = reorder_within(algorithm, by = value, within = measure, fun = "mean"), ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  scale_x_reordered() +  # scale properly for reordered algorithms
  coord_cartesian(ylim = c(0.66, 1)) +  # limit scale to magnify means and CIs
  labs(x = "Algorithm", y = "Mean rate") +  # label axes
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15)) +  # adjust strip font size
  ggsave("figures/accuracy_precision_recall_training.png", width = 9, height = 6)

## Testing
ggplot(acc_prec_recall_subject_testing, aes(x = reorder_within(algorithm, by = mean_value, within = measure, fun = "mean"), y = mean_value, fill = algorithm)) +
  facet_wrap(~measure, scale = "free_x") +  # create panel for each measure
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots without outliers
  geom_point(data = acc_prec_recall_means_testing, aes(x = reorder_within(algorithm, by = value, within = measure, fun = "mean"), y = value), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = acc_prec_recall_means_testing, aes(x = reorder_within(algorithm, by = value, within = measure, fun = "mean"), ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  scale_x_reordered() +  # scale properly for reordered algorithms
  coord_cartesian(ylim = c(0.5, 1)) +  # limit scale to magnify means and CIs
  labs(x = "Algorithm", y = "Mean rate") +  # label axes
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 15),  # adjust axis fonts
        strip.text = element_text(size = 12)) +  # adjust strip font size
  ggsave("figures/accuracy_precision_recall_testing.png", width = 9, height = 6)

# Data set 1
## Training
all_training_dataset1 <- filter(all_training, dataset == "Data set 1" & subject != "1-74") #%>%   # filter subject 1-74 because they have NAs, which prevent calculation of WSCIs
acc_prec_recall_training_dataset1 <- filter(all_training_dataset1, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Training
all_testing_dataset1 <- filter(all_testing, dataset == "Data set 1" & subject != "1-74") #%>%   # filter subject 1-74 because they have NAs, which prevent calculation of WSCIs
acc_prec_recall_testing_dataset1 <- filter(all_testing_dataset1, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Add within-subjects 95% confidence intervals
## Training
all_wsci_training_dataset1 <- wsci(all_training_dataset1, id = "subject", factors = c("algorithm", "measure", "type"), dv = "value")  # calculate 95% WSCI
all_means_training_dataset1 <- all_training_dataset1 %>%
  group_by(algorithm, measure) %>%  # for each algorithm and measure
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(all_wsci_training, by = c("algorithm", "measure")) %>%  # join the 95% WSCI data
  rename(ci = mean_value.y, value = mean_value.x) %>%  # rename CI column
  mutate(uci = value + ci,  # create upper CI column
         lci = value - ci)  # create lower CI column
acc_prec_recall_means_training_dataset1 <- filter(all_means_training_dataset1, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Testing
all_wsci_testing_dataset1 <- wsci(all_testing_dataset1, id = "subject", factors = c("algorithm", "measure", "type"), dv = "value")  # calculate 95% WSCI
all_means_testing_dataset1 <- all_testing_dataset1 %>%
  group_by(algorithm, measure) %>%  # for each algorithm and measure
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(all_wsci_testing, by = c("algorithm", "measure")) %>%  # join the 95% WSCI data
  rename(ci = mean_value.y, value = mean_value.x) %>%  # rename CI column
  mutate(uci = value + ci,  # create upper CI column
         lci = value - ci)  # create lower CI column
acc_prec_recall_means_testing_dataset1 <- filter(all_means_testing_dataset1, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Plot
## Training
acc_prec_recall_plot_training_dataset1 <- ggplot(acc_prec_recall_training_dataset1, aes(x = algorithm, y = value, fill = algorithm)) +
  facet_grid(type ~ measure, scale = "free_y") +  # create panel for each type and measure
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots
  geom_point(data = acc_prec_recall_means_training_dataset1, aes(x = algorithm, y = value), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = acc_prec_recall_means_training_dataset1, aes(x = algorithm, ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  coord_cartesian(ylim = c(0.6, 1)) +  # limit scale to magnify means and CIs
  labs(x = "Algorithm", y = "Mean rate") +  # label title and x axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        title = element_text(size = 15),  # adjust title fonts
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15))   # adjust strip font size

## Testing
acc_prec_recall_plot_testing_dataset1 <- ggplot(acc_prec_recall_testing_dataset1, aes(x = algorithm, y = value, fill = algorithm)) +
  facet_grid(type ~ measure, scale = "free_y") +  # create panel for each type and measure
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots
  geom_point(data = acc_prec_recall_means_testing_dataset1, aes(x = algorithm, y = value), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = acc_prec_recall_means_testing_dataset1, aes(x = algorithm, ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  coord_cartesian(ylim = c(0.6, 1)) +  # limit scale to magnify means and CIs
  labs(x = "Algorithm", y = "Mean rate") +  # label title and x axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        title = element_text(size = 15),  # adjust title fonts
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15))   # adjust strip font size

# Data set 2
## Training
all_training_dataset2 <- filter(all_training, dataset == "Data set 2") #%>%   # filter subject 1-74 because they have NAs, which prevent calculation of WSCIs
acc_prec_recall_training_dataset2 <- filter(all_training_dataset2, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Training
all_testing_dataset2 <- filter(all_testing, dataset == "Data set 2") #%>%   # filter subject 1-74 because they have NAs, which prevent calculation of WSCIs
acc_prec_recall_testing_dataset2 <- filter(all_testing_dataset2, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Add within-subjects 95% confidence intervals
## Training
all_wsci_training_dataset2 <- wsci(all_training_dataset2, id = "subject", factors = c("algorithm", "measure", "type"), dv = "value")  # calculate 95% WSCI
all_means_training_dataset2 <- all_training_dataset2 %>%
  group_by(algorithm, measure) %>%  # for each algorithm and measure
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(all_wsci_training, by = c("algorithm", "measure")) %>%  # join the 95% WSCI data
  rename(ci = mean_value.y, value = mean_value.x) %>%  # rename CI column
  mutate(uci = value + ci,  # create upper CI column
         lci = value - ci)  # create lower CI column
acc_prec_recall_means_training_dataset2 <- filter(all_means_training_dataset2, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Testing
all_wsci_testing_dataset2 <- wsci(all_testing_dataset2, id = "subject", factors = c("algorithm", "measure", "type"), dv = "value")  # calculate 95% WSCI
all_means_testing_dataset2 <- all_testing_dataset2 %>%
  group_by(algorithm, measure) %>%  # for each algorithm and measure
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(all_wsci_testing, by = c("algorithm", "measure")) %>%  # join the 95% WSCI data
  rename(ci = mean_value.y, value = mean_value.x) %>%  # rename CI column
  mutate(uci = value + ci,  # create upper CI column
         lci = value - ci)  # create lower CI column
acc_prec_recall_means_testing_dataset2 <- filter(all_means_testing_dataset2, measure %in% acc_prec_recall)  # filter accuracy, precision, and recall

## Plot
## Training
acc_prec_recall_plot_training_dataset2 <- ggplot(acc_prec_recall_training_dataset2, aes(x = algorithm, y = value, fill = algorithm)) +
  facet_grid(type ~ measure, scale = "free_y") +  # create panel for each type and measure
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots
  geom_point(data = acc_prec_recall_means_training_dataset2, aes(x = algorithm, y = value), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = acc_prec_recall_means_training_dataset2, aes(x = algorithm, ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  coord_cartesian(ylim = c(0.6, 1)) +  # limit scale to magnify means and CIs
  labs(x = "Algorithm", y = "Mean rate") +  # label title and x axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        title = element_text(size = 15),  # adjust title fonts
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15))   # adjust strip font size

## Testing
acc_prec_recall_plot_testing_dataset2 <- ggplot(acc_prec_recall_testing_dataset2, aes(x = algorithm, y = value, fill = algorithm)) +
  facet_grid(type ~ measure, scale = "free_y") +  # create panel for each type and measure
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots
  geom_point(data = acc_prec_recall_means_testing_dataset2, aes(x = algorithm, y = value), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = acc_prec_recall_means_testing_dataset2, aes(x = algorithm, ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  coord_cartesian(ylim = c(0.6, 1)) +  # limit scale to magnify means and CIs
  labs(x = "Algorithm", y = "Mean rate") +  # label title and x axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        title = element_text(size = 15),  # adjust title fonts
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15))   # adjust strip font size

## Combine data set figures
acc_prec_recall_plot_training_dataset1 / acc_prec_recall_plot_training_dataset2 + plot_annotation(tag_levels = 'A')
ggsave("figures/accuracy_precision_recall_training_datasets.png", width = 8, height = 11)
acc_prec_recall_plot_testing_dataset1 / acc_prec_recall_plot_testing_dataset2 + plot_annotation(tag_levels = 'A')
ggsave("figures/accuracy_precision_recall_testing_datasets.png", width = 8, height = 11)


## _Predictor importance -------
# Prepare all importance data
importance <- read_csv("stevens_etal_2020_data_importance.csv") %>%  # import importance data
  separate(data_file, into = c("type", "dataset")) %>%  # create separate columns for type and dataset
  mutate(type = fct_recode(type, "Amount" = "amount", "Delay" = "delay"),  # rename type
         dataset = fct_recode(dataset, "Data set 1" = "1", "Data set 2" = "2"),  # rename dataset
         predictor = fct_recode(predictor, "Difference" = "difference", "Discriminability" = "discrim", "Log Ratio" = "log_ratio", "Logistic" = "logistic", "Ratio" = "ratio", "Small" = "small", "Large" = "large"),  # rename predictors
         predictor = fct_relevel(predictor, "Small", "Large", "Log Ratio", "Ratio", "Discriminability", "Logistic", "Difference")
  ) %>%
  group_by(dataset, type, subject, predictor, algorithm) %>%  # group by type, dataset, predictor, subject
  summarise(mean_importance = mean(importance, na.rm = TRUE))  # calculate mean importance

# Overall
importance_subject <- importance %>%
  group_by(subject, predictor) %>%  # for each subject and predictor
  summarise(importance = mean(mean_importance, na.rm = TRUE))  # calculate mean importance

## Calculate 95% within-subjects confidence intervals
importance_wsci <- wsci(importance_subject, id = "subject", factors = "predictor", dv = "importance")
importance_means <- importance_subject %>%
  group_by(predictor) %>%  # for each algorithm and measure
  summarise(mean_importance = mean(importance, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(importance_wsci, by = c("predictor")) %>%  # join the 95% WSCI data
  rename(ci = importance) %>%  # rename CI column
  mutate(uci = mean_importance + ci,  # create upper CI column
         lci = mean_importance - ci)  # create lower CI column

## Plot
ggplot(importance_subject, aes(x = reorder(predictor, importance), y = importance)) +
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots without outliers
  geom_point(data = importance_means, aes(x = reorder(predictor, mean_importance), y = mean_importance), shape = 20, size = 4) +  # plot mean values per algorithm
  geom_linerange(data = importance_means, aes(x = reorder(predictor, mean_importance), ymin = lci, ymax = uci), inherit.aes = FALSE, size = 1) +  # plot WSCIs per algorithm
  labs(x = "Predictor", y = "Importance") +  # label x-axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),   # rotate x-axis labels
        axis.text.y = element_text(size = 15),   # adjust axis text fonts
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15)) +  # adjust strip font size
  ggsave("figures/predictor_importance.png", width = 9, height = 6)

# Data set 1
predictor_importance_dataset1 <- importance %>%
  filter(dataset == "Data set 1") %>%
  group_by(subject, predictor, type) %>%
  summarise(importance = mean(mean_importance, na.rm = TRUE))

## Calculate 95% within-subjects confidence intervals
importance_wsci_dataset1 <- wsci(predictor_importance_dataset1, id = "subject", factors = c("predictor", "type"), dv = "importance")
importance_means_dataset1 <- predictor_importance_dataset1 %>%
  group_by(predictor, type) %>%  # for each algorithm and measure
  summarise(mean_importance = mean(importance, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(importance_wsci_dataset1, by = c("predictor", "type")) %>%  # join the 95% WSCI data
  rename(ci = importance) %>%  # rename CI column
  mutate(uci = mean_importance + ci,  # create upper CI column
         lci = mean_importance - ci)  # create lower CI column

## Plot
predictor_importance_plot_dataset1 <- ggplot(predictor_importance_dataset1, aes(x = predictor, y = importance)) +  # reorder algorithms by mean importance within type
  facet_grid(~ type, scale = "free") +  # facet by type
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots without outliers
  geom_point(data = importance_means_dataset1, aes(x = predictor, y = mean_importance), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = importance_means_dataset1, aes(x = predictor, ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  # scale_x_reordered() +  # scale properly for reordered algorithms
  labs(x = "Predictor", y = "Importance") +  # label x-axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        title = element_text(size = 15),  # adjust title fonts
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # adjust axis text fonts
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15))   # adjust strip font size

# Data set 2
predictor_importance_dataset2 <- importance %>%
  filter(dataset == "Data set 2") %>%
  group_by(subject, predictor, type) %>%
  summarise(importance = mean(mean_importance, na.rm = TRUE))

## Calculate 95% within-subjects confidence intervals
importance_wsci_dataset2 <- wsci(predictor_importance_dataset2, id = "subject", factors = c("predictor", "type"), dv = "importance")
importance_means_dataset2 <- predictor_importance_dataset2 %>%
  group_by(predictor, type) %>%  # for each algorithm and measure
  summarise(mean_importance = mean(importance, na.rm = TRUE)) %>%  # calculate the mean value for each measure
  left_join(importance_wsci_dataset2, by = c("predictor", "type")) %>%  # join the 95% WSCI data
  rename(ci = importance) %>%  # rename CI column
  mutate(uci = mean_importance + ci,  # create upper CI column
         lci = mean_importance - ci)  # create lower CI column

## Plot
predictor_importance_plot_dataset2 <- ggplot(predictor_importance_dataset2, aes(x = predictor, y = importance)) +  # reorder algorithms by mean importance within type
  facet_grid(~ type, scale = "free") +  # facet by type
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", color = "gray50") +  # plot boxplots without outliers
  geom_point(data = importance_means_dataset2, aes(x = predictor, y = mean_importance), shape = 20) +  # plot mean values per algorithm
  geom_linerange(data = importance_means_dataset2, aes(x = predictor, ymin = lci, ymax = uci), inherit.aes = FALSE) +  # plot WSCIs per algorithm
  # scale_x_reordered() +  # scale properly for reordered algorithms
  labs(x = "Predictor", y = "Importance") +  # label x-axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        title = element_text(size = 15),  # adjust title fonts
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # adjust axis text fonts
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15))   # adjust strip font size

## Combine data set figures
predictor_importance_plot_dataset1 / predictor_importance_plot_dataset2 + plot_annotation(tag_levels = 'A')
ggsave("figures/predictor_importance_datasets.png", width = 8, height = 10)

# Algorithm
importance_algorithm <- importance %>%
  group_by(subject, predictor, algorithm) %>%  # for each subject and predictor
  summarise(importance = mean(mean_importance, na.rm = TRUE)) %>%  # calculate mean importance
  mutate(algorithm = fct_recode(algorithm, "C5.0" = "c5.0","kNN" = "knn", "Regression" = "glm", "Naive Bayes" = "naive_bayes", "Neural Network" = "nnet", "Random Forest" = "rf", "CART" = "rpart", "SVM" = "svm"),  # rename algorithms
         algorithm = fct_relevel(algorithm, "CART", after = 1),  # rename algorithms
         algorithm = fct_relevel(algorithm, "Regression", after = 6)  # rename algorithms
  )

## Plot
ann_text <- data.frame(predictor = "Difference", importance = 50, lab = "*", algorithm = factor("Regression", levels = levels(importance_algorithm$algorithm)))
importance_algorithm_plot <- ggplot(importance_algorithm, aes(x = predictor, y = importance)) +
  facet_wrap(~algorithm) +
  stat_summary(fun.data = calc_boxplot_stat, geom = "boxplot", color = "gray50") +  # plot boxplots without outliers
  stat_summary(fun.data = "mean_cl_normal", geom = "point", shape = 20) +  # plot means and 95% CIs
  stat_summary(fun.data = "mean_cl_normal", geom = "linerange") +  # plot means and 95% CIs
  geom_text(data = ann_text, label = "*", size = 6) +
  labs(x = "Predictor", y = "Importance") +  # label x-axis
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        legend.position = "none",  # remove legend
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),   # rotate x axis labels
        axis.text.y = element_text(size = 12),   # adjust axis text fonts
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15)) +   # adjust strip font size
  ggsave("figures/predictor_importance_algorithm.png", width = 7, height = 6)


## _Sample size effects on accuracy -------

# Prepare all sample size data
accuracy_testing <- bind_rows(accuracy_random_all, accuracy_sequential_all) %>%  # combine random and sequential data
  filter(train_test == "testing") %>%  # filter testing data
  group_by(dataset, order, subject, sample_size, algorithm, type) %>%  # for each data set, order, subject, sample_size, and algorithm
  summarise(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>%   # calculate mean accuracy
  ungroup(order) %>%  # ungroup to work with order
  mutate(order = factor(order, levels = c("Random", "Sequential"))) %>%  # reorder levels
  drop_na()  # remove NaNs

# Overall
accuracy_overall <- accuracy_testing %>%
  group_by(order, subject, sample_size, algorithm) %>%  # for each order, subject, sample size, and algorithm
  summarise(accuracy = mean(mean_accuracy, na.rm = TRUE))   # calculate mean accuracy
accuracy_means <- accuracy_overall %>%
  group_by(order, algorithm, sample_size) %>%  # for each order, algorithm, and sample size
  summarise(mean_pred = mean(accuracy))  # calculate mean accuracy

## Plot
ggplot(accuracy_overall, aes(x = sample_size, y = accuracy, group = algorithm, color = algorithm, linetype = algorithm)) +
  facet_grid(~order) +  # create panel for each order
  stat_summary(fun.data = "mean_cl_normal", geom = "point", shape = 20) +  # plot means and 95% CIs
  stat_summary(fun.data = "mean_cl_normal", geom = "linerange", linetype = 1) +  # plot means and 95% CIs
  stat_summary(fun = "mean", geom = "line") +  # plot lines connecting means
  labs(x = "Sample size", y = "Mean accuracy rate") +  # label axes
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        axis.text = element_text(size = 15),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 18),  # adjust strip font size
        legend.title = element_blank(),  # remove legend title
        legend.key.size = unit(0.9, 'lines'),  # change line spacing
        legend.text = element_text(size = 15),  # change legend text size
        legend.position = c(0.85, 0.25)  # position legend
  ) +
  ggsave("figures/accuracy_samplesize.png", width = 9, height = 6)

# Random
accuracy_samplesize_random <- accuracy_testing %>%
  filter(order == "Random") %>%  # keep only testing data
  ggplot(aes(x = sample_size, y = mean_accuracy, group = algorithm, color = algorithm, linetype = algorithm)) +
  facet_grid(dataset ~ type) +  # create panel for each dataset and type
  stat_summary(fun.data = "mean_cl_normal", geom = "point", shape = 20) +  # plot means and 95% CIs
  stat_summary(fun.data = "mean_cl_normal", geom = "linerange", linetype = 1) +  # plot means and 95% CIs
  stat_summary(fun = "mean", geom = "line") +  # plot lines connecting means
  coord_cartesian(ylim = c(0.4, 1)) +  # limit scale
  labs(x = "Sample size", y = "Mean accuracy rate") +  # label axes
  guides(fill = guide_legend(nrow =  4)) +
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        title = element_text(size = 15),  # adjust title fonts
        axis.text = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15),  # adjust strip font size
        legend.title = element_blank(),  # remove legend title
        legend.key.size = unit(0.9, 'lines'),  # change line spacing
        legend.text = element_text(size = 10),  # change legend text size
  )   # position legend

# Sequential
accuracy_samplesize_sequential <- accuracy_testing %>%
  filter(order == "Sequential") %>%  # keep only testing data
  ggplot(aes(x = sample_size, y = mean_accuracy, group = algorithm, color = algorithm, linetype = algorithm)) +
  facet_grid(dataset ~ type) +  # create panel for each dataset and type
  stat_summary(fun.data = "mean_cl_normal", geom = "point", shape = 20) +  # plot means and 95% CIs
  stat_summary(fun.data = "mean_cl_normal", geom = "linerange", linetype = 1) +  # plot means and 95% CIs
  stat_summary(fun = "mean", geom = "line") +  # plot lines connecting means
  coord_cartesian(ylim = c(0.4, 1)) +  # limit scale
  labs(x = "Sample size", y = "Mean accuracy rate") +  # label axes
  guides(fill = guide_legend(nrow =  4)) +
  theme_bw() +  # set theme
  theme(panel.grid = element_blank(),  # remove grid lines
        title = element_text(size = 15),  # adjust title fonts
        axis.text = element_text(size = 12),   # rotate x axis labels
        axis.title = element_text(size = 25),  # adjust axis fonts
        strip.text = element_text(size = 15),  # adjust strip font size
        legend.title = element_blank(),  # remove legend title
        legend.key.size = unit(0.9, 'lines'),  # change line spacing
        legend.text = element_text(size = 10),  # change legend text size
  )   # position legend

# Combine order figures
accuracy_samplesize_order <- accuracy_samplesize_random / accuracy_samplesize_sequential + plot_annotation(tag_levels = 'A')
ggsave("figures/accuracy_samplesize_order.png", width = 8, height = 11)



# Tables -------------------------------------------------------------------
accuracy_random_all <- accuracy_random_all %>% 
  mutate(prop_tp = true_positive / (true_positive + true_negative + false_positive + false_negative) * 100,  # calculate true positive percentage
         prop_tn = true_negative / (true_positive + true_negative + false_positive + false_negative) * 100,  # calculate true negative percentage
         prop_fp = false_positive / (true_positive + true_negative + false_positive + false_negative) * 100,  # calculate false positive percentage
         prop_fn = false_negative / (true_positive + true_negative + false_positive + false_negative) * 100,  # calculate false negative percentage
         specificity = true_negative / (true_negative + false_positive),  # calculate specificity
         negative_pred_val = true_negative / (true_negative + false_negative)  # calculate negative predictive value
         )

confusion_matrix_means <- accuracy_random_all %>% 
  filter(sample_size == 30) %>% 
  group_by(train_test, algorithm) %>% 
  summarize(true_positive = mean(prop_tp, na.rm = TRUE),
            true_negative = mean(prop_tn, na.rm = TRUE),
            false_positive = mean(prop_fp, na.rm = TRUE),
            false_negative = mean(prop_fn, na.rm = TRUE),
            acc = mean(accuracy, na.rm = TRUE) * 100,
            rec = mean(recall, na.rm = TRUE) * 100,
            spec = mean(specificity, na.rm = TRUE) *100,
            prec = mean(precision, na.rm = TRUE) * 100,
            neg_pred_val = mean(negative_pred_val, na.rm = TRUE) * 100
  ) %>% 
  mutate(train_test = str_to_sentence(train_test),  # capitalize training and testing
         algorithm = factor(algorithm, levels = c("C5.0", "CART", "kNN", "Naive Bayes", "Neural Network", "Random Forest", "Regression", "SVM"))
  ) %>% 
  arrange(desc(train_test), algorithm)
