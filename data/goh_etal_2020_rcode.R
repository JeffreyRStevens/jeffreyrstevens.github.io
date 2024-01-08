## ---
##
## Script name: goh_etal_2020_rcode.R
##
## Purpose of script: This script analyzes tip sizes and feelings towards establishments 
##  that use tip screens vs. tip jars
##
## Authors: Francine Goh (francinegoh@gmail.com), Alexandria Jungck, 
##  Jeffrey Stevens (jeffrey.r.stevens@gmail.com)
##
## Date Finalized: 26 Apr 2020 
##
## License: This script is released under the Creative Commons 
##   Attribution-NonCommercial-ShareAlike 4.0 International license (CC BY-NC-SA 4.0). 
##   You may share and adapt this content with attribution, for non-commercial purposes 
##   if you cite and ShareAlike (distribute any contributions under the same license).
##
## ---
##
## Notes:
##   
## Instructions: Place this file and the data files in the main directory.
## 	Create folders called "data" and "figures". 
##  Set the R working directory to the main directory.  At the R command prompt, type
## 	> source("goh_etal_2020_rcode.R")
## 	This will run the script, adding all of the calculated variables to the workspace and 
##  saving figures in the figures directory. If packages do not load properly, install them 
##  with install.packages("package_name").
## Data files:
##  goh_etal_2020_data.csv
##   subject_nr = participant number
##   gender = participant gender
##   age = participant age
##   race = participant ethnicity
##   bp_ts = tip amount for tip screen, barista present condition
##   ba_ts = tip amount for tip screen, barista absent condition
##   bp_rec = tip amount for receipt, barista present condition
##   ba_rec = tip amount for receipt, barista absent condition
##   bp_tj = tip amount for tip jar, barista present condition
##   ba_tj = tip amount for tip jar, barista absent condition
##   feel_ts = participant's degree of negative feelings towards tip screens
##   feel_tj = participant's degree of negative feelings towards tip jars
##   avoid_ts = participant's frequency of avoidance of tip screens
##   avoid_tj = participant's frequency of avoidance of tip jars
##   EQ_mean = participant's mean empathy score
##
## ---


# Load libraries and functions --------------------------------------------
library(BayesFactor)
library(broom)
library(car)
library(lme4)
library(lsr)
library(moments)
library(rptR)
library(tidyverse)
library(patchwork)


# Convert BIC values to Bayes factor
bic_bf <- function(null, alternative) {
  new_bf <- exp((null - alternative) / 2) # convert BICs to Bayes factor
  names(new_bf) <- NULL  # remove BIC label
  return(new_bf)  # return Bayes factor of alternative over null hypothesis
}

# Create qqplots using ggplot2 with a linear model
gg_qqplot <- function(LM) # argument: a linear model
{
  y <- quantile(resid(LM), c(0.25, 0.75))  # extract residuals and calculate quantiles of data
  x <- qnorm(c(0.25, 0.75))  # calculate quantiles of normal distribution
  slope <- diff(y) / diff(x)  # calculate slope of expected normal distribution line
  int <- y[1L] - slope * x[1L]  # calculate intercept of expected normal distribution line
  p <- ggplot(LM, aes(sample = .resid)) +
    stat_qq(alpha = 0.5) +  # plot the residuals against standard normal distribution
    geom_abline(slope = slope, intercept = int, color="blue")  # plot expected normal distribution line
  return(p)
}

set.seed(0) # for consistency in model algorithm computations


# Import data -------------------------------------------------------------
all_data <- read_csv("goh_etal_2020_data.csv", col_names = TRUE)  # import data


# Feelings towards tip screens and tip jars --------------
# Create dataset for feelings towards tipscreen and tipjar
feelings_data <- all_data %>% 
  filter(!is.na(feel_ts), !is.na(feel_tj)) %>%  # remove participants who have NA responses for feeling questions
  select(subject_nr, feel_ts, feel_tj) %>%  # select these columns from all_data
  gather(feel_ts, feel_tj, key = "screen_jar", value = "feelings_rating") %>%  # sort participant responses by payment type (tipscreen/tipjar)
  mutate(payment_type = ifelse(grepl("ts", screen_jar), "tipscreen",  # insert "tipscreen" label in "payment_type" column if "screen_jar" column contains "ts"
                               ifelse(grepl("tj", screen_jar), "tipjar", NA)))  # insert "tipjar" label in "payment_type" column if "screen_jar" column contains "tj"

# Create boxplot to compare participants' feelings towards tip screens & tip jars (higher ratings = more negative feelings toward payment type)
feeling_plot <- ggplot(feelings_data, aes(x = payment_type, y = feelings_rating, fill = payment_type)) +  # create payment type & feelings rating graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Degree of negative feelings") +  # rename x & y axes labels
  scale_x_discrete(labels = c("tipjar" = "Tip jar", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  scale_y_continuous(limits = c(1, 6)) +  # specify y-axis limits
  scale_fill_manual(values=c("#0C7BDC", "#FFC20A")) +  # change boxplot fill colors
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 30, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# t-test & Bayes analysis: compare participant feelings towards tip jars vs. tip screens
feelings_data_analyses <-  feelings_data %>%  # create dataset for analyses
  select(subject_nr, screen_jar, feelings_rating) %>%  # select columns
  pivot_wider(names_from = screen_jar, values_from = feelings_rating)  # pivot data to wide format

t.test(feelings_data_analyses$feel_tj, feelings_data_analyses$feel_ts, paired = TRUE)   # paired t-test: compare feelings towards tip jars vs. tip screens
ttestBF(feelings_data_analyses$feel_tj, feelings_data_analyses$feel_ts, paired = TRUE)  # Bayes factor analysis
cohensD(feelings_data_analyses$feel_tj, feelings_data_analyses$feel_ts)  # calculate effect size (Cohen's d) for feelings data


# Avoidance of tip screens and tip jars --------------
# Create dataset for avoidance of tipscreen and tipjar
avoidance_data <- all_data %>% 
  filter(!is.na(avoid_ts), !is.na(avoid_tj)) %>%  # remove participants who have NA responses for avoidance questions
  select(subject_nr, avoid_ts, avoid_tj) %>%  # select these columns from all_data
  gather(avoid_ts, avoid_tj, key = "screen_jar", value = "avoidance_frequency") %>%  # sort participant responses by payment type (tip screen/tip jar)
  mutate(payment_type = ifelse(grepl("ts", screen_jar), "tipscreen",  # insert "tipscreen" label in "payment_type" column if "screen_jar" column contains "ts"
                               ifelse(grepl("tj", screen_jar), "tipjar", NA)))  # insert "tipjar" label in "payment_type" column if "screen_jar" column contains "tj"

# Create boxplot to compare participants' frequency of avoidance of tipscreens & tipjars (Qualtrics response options key: 1 = Never, 2 = Once, 3 = 2-5 times, 4 = 6-10 times, 5 = more than 10 times)
avoidance_plot <- ggplot(avoidance_data, aes(x = payment_type, y = avoidance_frequency, fill = payment_type)) +  # create payment type & frequency of avoidance of tipscreen/tipjar graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Frequency of avoidance") +  # rename x & y axes labels
  scale_x_discrete(labels = c("tipjar" = "Tip jar", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#0C7BDC", "#FFC20A")) +  # change boxplot fill colors
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 30, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# t-test & Bayes analysis: compare participant avoidance of tip jars vs. tip screens
avoidance_data_analyses <-  avoidance_data %>%  # create dataset for analyses
  select(subject_nr, screen_jar, avoidance_frequency) %>%  # select columns
  pivot_wider(names_from = screen_jar, values_from = avoidance_frequency)  # pivot data to wide format

t.test(avoidance_data_analyses$avoid_tj, avoidance_data_analyses$avoid_ts, paired = TRUE)   # paired t-test: compare avoidance of tip jars vs. tip screens
ttestBF(avoidance_data_analyses$avoid_tj, avoidance_data_analyses$avoid_ts, paired = TRUE)  # Bayes factor analysis
cohensD(avoidance_data_analyses$avoid_tj, avoidance_data_analyses$avoid_ts)  # calculate effect size (Cohen's d) for avoidance data

# Combine subfigures to one figure
(feeling_plot + avoidance_plot) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 25))
ggsave("figures/feeling_avoidance_payment_method.png", width = 13, height = 6.5)

# Tip amounts ---------------
# Create dataset for participant tip amounts, barista conditions & payment type
tipping_data <- all_data %>% 
  select(subject_nr, bp_ts:ba_tj) %>%   # select these columns from all_data 
  pivot_longer(bp_ts:ba_tj, names_to = "condition", values_to = "response") %>%   # sort participant responses by barista condition & payment condition
  mutate(barista = ifelse(grepl("bp", condition), "Barista present",  # insert "Barista present" label in "barista" column if tipping condition contains "bp" in "condition" column
                          ifelse(grepl("ba", condition), "Barista absent", NA)),  # insert "Barista absent" label in "barista" column if tipping condition contains "ba" in "condition" column
         payment_type = ifelse(grepl("ts", condition), "Tip screen",  # insert "Tip screen" in "payment_type" column if tipping condition contains "ts" in "condition" column
                               ifelse(grepl("rec", condition), "Receipt",  # insert "Receipt" in "payment_type" column if tipping condition contains "rec" in "condition" column
                                      ifelse(grepl("tj", condition), "Cash", NA))),  # insert "Cash" in "payment_type" column if tipping condition contains "tj" in "condition" column
         condition = as.factor(condition),   # convert variables to factor type
         barista = as.factor(barista),
         payment_type = as.factor(payment_type),  
  )

tipping_data_noNA <- tipping_data[complete.cases(tipping_data), ]  # remove rows with missing participant responses

# Calculate summary statistics
tip_size_mean <- mean(tipping_data_noNA$response)  # calculate mean tip size
tip_size_sd <- sd(tipping_data_noNA$response)  # calculate standard deviation for tip size

# Identify outliers (i.e., data points 3 SD above/below mean tip size)
cut_off <- tip_size_sd * 3
outliers <- tipping_data_noNA %>%
  mutate(lower = tip_size_mean - cut_off,  # calculate lower range cutoff
         upper = tip_size_mean + cut_off,  # calculate upper range cutoff
         outlier_point = ifelse(response <= lower, 1, 0) | ifelse(response >= upper, 1, 0)) %>% 
  filter(outlier_point == TRUE)  # select outlier data rows

tipping_data_clean <- tipping_data_noNA %>% 
  anti_join(outliers)  # remove outliers


# _Plots: barista condition and payment type effects on tip amounts --------
# Plot barista condition effect on tip amounts
## Calculate participants' mean tip amounts for each barista (absent/present) condition
barista_condition <- tipping_data_clean %>% 
  group_by(subject_nr, barista) %>% 
  summarise(tip_mean = mean(response, na.rm = TRUE))  # calculate participants' mean tip amounts for each barista (absent/present) condition

## Create boxplot to compare participants' mean tip amounts by barista condition
barista_plot <- ggplot(barista_condition, aes(x = barista, y = tip_mean, fill = barista)) +  # create barista absent/present & tip_mean graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Barista presence", y = "Mean tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("Barista absent" = "Absent", "Barista present" = "Present")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#f02525", "#25f0f0")) +  # change boxplot fill colors
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# Plot payment type effect on tip amounts
## Calculate participants' mean tip amounts for each payment type (credit card, receipt, cash)
payment_condition <- tipping_data_clean %>% 
  group_by(subject_nr, payment_type) %>% 
  summarise(tip_mean = mean(response, na.rm = TRUE))  # calculate participants' mean tip amounts for each payment type (credit card, receipt, cash)

## Create boxplot to compare participants' mean tip amounts by payment type
payment_plot <- ggplot(payment_condition, aes(x = payment_type, y = tip_mean, fill = payment_type)) +  # create payment type & tip_mean graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Mean tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("cash" = "Cash", "receipt" = "Receipt", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# Combine subfigures to one figure
(payment_plot + barista_plot) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 25))
ggsave("figures/tipping_payment_method_barista.png", width = 15, height = 6.5)


# _Check normality assumptions for tip amounts ----
# Calculate skewness & kurtosis of raw tip amounts
skewness(tipping_data_clean$response)  # calculate skewness of raw tip amounts
kurtosis(tipping_data_clean$response)  # calculate kurtosis of raw tip amounts

# Transform & center tip amounts to reduce model convergence issues
tipping_data_centered <- tipping_data_clean %>% 
  mutate(response_sqrt = sqrt(response))  # square root transform response

# Calculate skewness & kurtosis for square root transformed tip amounts
skewness(tipping_data_centered$response_sqrt, na.rm = TRUE)  # calculate skewness of square root transformed tip amounts
kurtosis(tipping_data_centered$response_sqrt, na.rm = TRUE)  # calculate kurtosis of square root transformed tip amounts

# Center tip amounts
tipping_data_centered <- tipping_data_centered %>%
  select(-contains("condition"), -response) %>%   # remove condition and response columns
  pivot_wider(names_from = payment_type, values_from = response_sqrt) %>%  # pivot data to wide format
  mutate(cash_cen = Cash - mean(Cash, na.rm = TRUE),  # group mean center tip amount for cash payment type
         receipt_cen = Receipt - mean(Receipt, na.rm = TRUE),  # group mean center tip amount for receipt payment type
         tipscreen_cen = `Tip screen` - mean(`Tip screen`, na.rm = TRUE)) %>%  # group mean center tip amount for tipscreen payment type
  select(-`Tip screen`, -Receipt, -Cash) %>%  # remove uncentered data
  pivot_longer(cash_cen:tipscreen_cen, names_to = "payment_type", values_to = "response_cen") %>%  # pivot data to long format
  unite("condition", barista:payment_type, remove = FALSE) %>%  # create condition column
  mutate(condition = as.factor(condition)) %>%   # convert to factor
  drop_na()


# _Intraclass correlation --------------------------------------------------
# Calculate intraclass correlation (ICC): what is the correlation between two randomly selected tip amounts made by a participant?
icc_tipamount <- lmer(response_cen ~ (1 | subject_nr), data = tipping_data_centered)  # conduct linear mixed model with participant as random effect (i.e., random intercept model)
icc_tipamount_CI <- rpt(response_cen ~ (1 | subject_nr), grname = "subject_nr", data = tipping_data_centered, datatype = "Gaussian",  nboot = 1000, npermut = 0)  # calculate ICC and 95% CIs (bootstrap)

icc_value <- as.numeric(icc_tipamount_CI$R*100)  # save ICC value
icc_lower <- icc_tipamount_CI$CI_emp$`2.5%` * 100  # save ICC lower 95% CI
icc_upper <- icc_tipamount_CI$CI_emp$`97.5%` * 100  # save ICC upper 95% CI
# Result: In an empty, random intercept model, 47.12% [40.79, 53.02] of the variation in tip amount is accounted for by individual differences among participants (i.e., participants tipped differently from one another across the six tipping conditions).


# _Model selection ---------------------------------------------------------
# First find best-fitting random effect model, then test for fixed effect predictors
# Process is backward by elimintating weakest terms sequentially, starting with full model, until only significant effects remain
# Nested model comparisons (likelihood ratio tests using anova command) are used to select best fitting models

# Random effects
## Create random effects models
barista_payment_random_full <- lmer(response_cen ~ (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered)  # full random model (i.e., include random intercepts and slopes for barista presence and payment type); fail to converge
barista_payment_random_full <- lmer(response_cen ~ (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered, control = lmerControl(optimizer = "bobyqa"))  # full random model (i.e., include random intercepts and slopes for barista presence and payment type); use bobyqa optimizer; BEST
barista_payment_random_barista_only <- lmer(response_cen ~ (1 + barista | subject_nr), data = tipping_data_centered)  # drop payment type intercept and slope
barista_payment_random_payment_only <- lmer(response_cen ~ (1 + payment_type | subject_nr), data = tipping_data_centered)  # drop barista intercept and slope; fail to converge
barista_payment_random_payment_only <- lmer(response_cen ~ (1 + payment_type | subject_nr), data = tipping_data_centered, control = lmerControl(optimizer = "bobyqa"))  # drop barista intercept and slope; use bobyqa optimizer
barista_payment_random_intercept_only <- lmer(response_cen ~ (1 | subject_nr), data = tipping_data_centered)  # drop random slopes (i.e., include only intercept of random factor subject_nr in model)

## Likelihood ratio tests for model comparison
barista_payment_random_anova <- anova(barista_payment_random_full, barista_payment_random_barista_only, barista_payment_random_payment_only, barista_payment_random_intercept_only)
barista_payment_random_anova_tidy <- tidy(barista_payment_random_anova)  # create tidy table of model parameters 
## Result: The full random effects model with random intercepts and slopes for barista presence and payment type is better than the other models

# Fixed effects
## Create fixed-effects models
barista_payment_full <- lmer(response_cen ~ barista * payment_type + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered) # full model; fail to converge
barista_payment_full <- lmer(response_cen ~ barista * payment_type + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered, control = lmerControl(optimizer = "bobyqa")) # full model; use bobyqa optimizer
barista_payment_fixed <- lmer(response_cen ~ barista + payment_type + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered) # drop interaction first, then drop weakest variables; fail to converge
barista_payment_fixed <- lmer(response_cen ~ barista + payment_type + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered) # drop interaction first, then drop weakest variables; use bobyqa optimizer; fail to converge
barista_payment_fixed2 <- lmer(response_cen ~ barista + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered)  # drop payment type; fail to converge
barista_payment_fixed2 <- lmer(response_cen ~ barista + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered, control = lmerControl(optimizer = "bobyqa")) # drop payment type; use bobyqa optimizer
barista_payment_fixed3 <- lmer(response_cen ~ payment_type + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered)  # drop barista condition; fail to converge
barista_payment_fixed3 <- lmer(response_cen ~ payment_type + (1 + barista | subject_nr) + (1 + payment_type | subject_nr), data = tipping_data_centered, control = lmerControl(optimizer = "bobyqa"))  # drop barista condition; use bobyqa optimizer

## Likelihood ratio tests for model comparison
barista_payment_fixed_anova <- anova(barista_payment_full, barista_payment_fixed2, barista_payment_fixed3, barista_payment_random_full)
barista_payment_fixed_anova_tidy <- tidy(barista_payment_fixed_anova)  # create tidy table of model parameters 
# Result: The fixed effects model with barista only is better than the other models

## Bayes factors for fixed effects
### Extract BICs
barista_payment_random_bic <- barista_payment_fixed_anova_tidy$BIC[1]  # empty random model
barista_payment_full_bic <- barista_payment_fixed_anova_tidy$BIC[4]  # full random model
barista_payment_fixed2_bic <- barista_payment_fixed_anova_tidy$BIC[2]  # dropped payment type model
barista_payment_fixed3_bic <- barista_payment_fixed_anova_tidy$BIC[3]  # dropped barista condition model

### Convert BICs to BFs compared to empty random model
barista_payment_full_bf <- bic_bf(barista_payment_random_bic, barista_payment_full_bic)
barista_payment_fixed2_bf <- bic_bf(barista_payment_random_bic, barista_payment_fixed2_bic)
barista_payment_fixed3_bf <- bic_bf(barista_payment_random_bic, barista_payment_fixed3_bic)

## Create model comparison table
### Random effects model comparison
barista_payment_random_anova_table <- barista_payment_random_anova_tidy %>%  
  mutate("Model specification" = c("Participant", "Participant + barista presence slope", "Participant + payment type slope", "Participant + barista presence slope + payment type slope"),  # create specification column
         "Fixed effects" = rep("-", 4),  # create fixed effects column
         "Random effects" = c("(1 | subject_nr)", "(1 + barista | subject_nr)", "(1 + payment_type | subject_nr)", "(1 + barista | subject_nr) + (1 + payment_type | subject_nr)"), # create random effects column
         BF = NA) %>%  # create Bayes factor column
  select("Model specification", "Random effects", "Fixed effects", AIC, BIC, logLik, df, statistic, df, "p-value" = p.value, BF)  # select and rearrange columns

### Fixed effects model comparision
barista_payment_fixed_anova_table <- barista_payment_fixed_anova_tidy %>%  
  mutate("Model specification" = c("RE only", "Barista presence", "Payment type", "Barista presence * payment type"),  # create specification column
         "Fixed effects" = c("-", "barista", "payment_type", "barista * payment_type"),  # create fixed effects column
         "Random effects" = c("(1 + barista | subject_nr) + (1 + payment_type | subject_nr)", "(1 + barista | subject_nr) + (1 + payment_type | subject_nr)", "(1 + barista | subject_nr) + (1 + payment_type | subject_nr)", "(1 + barista | subject_nr) + (1 + payment_type | subject_nr)"),  # create random effects column
         BF = c(NA, barista_payment_fixed2_bf, barista_payment_fixed3_bf, barista_payment_full_bf)) %>%   # create Bayes factor column
  select("Model specification", "Random effects", "Fixed effects", AIC, BIC, logLik, df, statistic, df, "p-value" = p.value, BF)  # select and rearrange columns

### Combine tables
rand_row <- tibble("Model specification" = "Random effect models")  # create row for random effect models
fixed_row <- tibble("Model specification" = "Fixed effect models")  # create row for fixed effect models
barista_payment_anova_table <- bind_rows(rand_row, barista_payment_random_anova_table, fixed_row, barista_payment_fixed_anova_table)  # combine tables

## Check assumptions
boxplot(residuals(barista_payment_fixed2) ~ tipping_data_centered$condition) # view distributions
leveneTest(residuals(barista_payment_fixed2) ~ tipping_data_centered$condition) # check assumption of equal variances
plot(barista_payment_fixed2)  # plot residuals vs. predicted values
plot(density(residuals(barista_payment_fixed2)))  # plot distribution of residuals
gg_qqplot(barista_payment_fixed2) # plot qqplot


# Empathy ------------------
# Prepare data
outliers_subjects <- outliers %>% 
  distinct(subject_nr)  # give only subject_nr of outliers

tipping_data_clean_eq <-  tipping_data_clean %>%  # create dataset to merge with EQ data
  select(subject_nr, condition:payment_type)

eq_data <- all_data %>% 
  filter(!is.na(EQ_mean), !is.na(bp_ts), !is.na(ba_ts), !is.na(bp_rec), !is.na(ba_rec), !is.na(bp_tj), !is.na(ba_tj)) %>%  # remove participants who do not have EQ_mean scores, tip amounts for the 6 tipping conditions
  select(subject_nr, bp_ts:ba_tj, EQ_mean) %>%  # select only these columns
  pivot_longer(bp_ts:ba_tj, names_to = "condition", values_to = "response") %>%
  left_join(tipping_data_clean_eq, by = c("subject_nr", "condition", "response")) %>%  # combine datasets
  anti_join(outliers_subjects) %>%  # remove all cases for outlier participants (i.e., not just outlier response values)
  mutate(response_sqrt = sqrt(response)) %>%  # square root transform response variable to reduce model convergence issues
  select(-condition, -response) %>%   # remove condition and response columns
  pivot_wider(names_from = payment_type, values_from = response_sqrt) %>%
  mutate(cash_cen = Cash - mean(Cash, na.rm = TRUE),  # group mean center tip amount for cash payment type
         receipt_cen = Receipt - mean(Receipt, na.rm = TRUE),  # group mean center tip amount for receipt payment type
         tipscreen_cen = `Tip screen` - mean(`Tip screen`, na.rm = TRUE)) %>%  # group mean center tip amount for tipscreen payment type
  select(-`Tip screen`, -Receipt, -Cash) %>%  # remove uncentered data
  pivot_longer(cash_cen:tipscreen_cen, names_to = "payment_type", values_to = "response_cen") %>%   # remove extra rows
  group_by(subject_nr, EQ_mean, barista) %>% 
  summarise(response_mean = mean(response_cen, na.rm = TRUE))  # calculate participants' mean tip amounts for each barista (absent/present) condition

eq_data_plot <- eq_data %>%
  left_join(barista_condition, by = c("subject_nr", "barista")) # combine eq_data with barista_condition


# _Plot EQ data ------------------------------------------------------------
# Create scatterplot to compare EQ effects on participant mean tip amounts by barista condition
ggplot(eq_data_plot, aes(x = EQ_mean, y = tip_mean)) +  # create mean EQ scores & mean tip amount axes labels
  geom_jitter(width = 0.01, height = 0.01) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~barista) +
  labs(x = "Mean empathy score", y = "Mean tip size ($)") +  # rename x & y axes labels
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        strip.text.x = element_text(size = 30),  # change facet label size
        legend.position = "none")  # remove legend from graph
ggsave("figures/tipping_empathy_barista.png", width = 13.8, height = 6.5)  # save plot


# _Model selection ---------------------------------------------------------
# First find best-fitting random effect model, then test for fixed effect predictors
# Process is backward by elimintating weakest terms sequentially, starting with full model, until only significant effects remain
# Nested model comparisons (likelihood ratio tests using anova command) are used to select best fitting models

# Random effects
## Create random effects models
# For the full random effects model (see below), the random effects variances are confounded with the residual variance because trying to estimate 2 random effects per subject_nr when there are only 2 estimations per subject_nr. Therefore, we cannot use this model
# eq_barista_random_full <- lmer(response_mean ~ (1 + barista | subject_nr) + (1 | EQ_mean), data = eq_data)  # full random model (i.e., include random intercept and slope for barista presence)
eq_barista_random2 <- lmer(response_mean ~ (1 | subject_nr) + (1 | EQ_mean), data = eq_data)  # drop random slope
eq_barista_random2 <- lmer(response_mean ~ (1 | subject_nr) + (1 | EQ_mean), data = eq_data, control = lmerControl(optimizer = "bobyqa"))  # drop random slope; use bobyqa optimizer; overfit
eq_barista_random3 <- lmer(response_mean ~ (1 | subject_nr), data = eq_data)  # drop mean EQ; BEST
eq_barista_random4 <- lm(response_mean ~ 1, data = eq_data) # empty model

## Likelihood ratio tests for model comparison
eq_barista_random_anova <- anova(eq_barista_random2, eq_barista_random3, eq_barista_random4)  
eq_barista_random_anova_tidy <- tidy(eq_barista_random_anova)  # create tidy table of model parameters 
## Result: A random intercept for each participant is better than the empty model

# Fixed effects
## Create fixed-effects models
eq_barista_full <- lmer(response_mean ~ barista * EQ_mean + (1 | subject_nr), data = eq_data) # full model
eq_barista_fixed <- lmer(response_mean ~ barista + EQ_mean + (1 | subject_nr), data = eq_data) # drop interaction first, then drop weakest variables
eq_barista_fixed2 <- lmer(response_mean ~ barista + (1 | subject_nr), data = eq_data)  # drop mean EQ; BEST
eq_barista_fixed3 <- lmer(response_mean ~ EQ_mean + (1 | subject_nr), data = eq_data)  # drop barista condition

## Likelihood ratio tests for model comparison
## Including mean EQ does not improve model fit, therefore model without mean EQ is chosen for parsimony
eq_barista_fixed_anova <- anova(eq_barista_full, eq_barista_fixed, eq_barista_fixed2, eq_barista_fixed3, eq_barista_random3)  
eq_barista_fixed_anova_tidy <- tidy(eq_barista_fixed_anova)  # create tidy table of model parameters 

## Compare models: barista condition only vs. barista condition + EQ
eq_barista_fixed_comparison <- anova(eq_barista_fixed, eq_barista_fixed2)  
eq_barista_fixed_comparison_tidy <- tidy(eq_barista_fixed_comparison)  # create tidy table of model parameters 
# Result: The fixed effects model with barista only (no EQ) is better than the other models

## Bayes factors for fixed effects
### Extract BICs
eq_barista_random_bic <- eq_barista_fixed_anova_tidy$BIC[1]  # random model
eq_barista_full_bic <- eq_barista_fixed_anova_tidy$BIC[5]  # full random model
eq_barista_fixed_bic <- eq_barista_fixed_anova_tidy$BIC[4]  # dropped interaction model
eq_barista_fixed2_bic <- eq_barista_fixed_anova_tidy$BIC[2]  # dropped mean EQ model
eq_barista_fixed3_bic <- eq_barista_fixed_anova_tidy$BIC[3]  # dropped barista condition model

### Convert BICs to BFs
eq_barista_full_bf <- bic_bf(eq_barista_random_bic, eq_barista_full_bic)
eq_barista_fixed_bf <- bic_bf(eq_barista_random_bic, eq_barista_fixed_bic)
eq_barista_fixed2_bf <- bic_bf(eq_barista_random_bic, eq_barista_fixed2_bic)
eq_barista_fixed3_bf <- bic_bf(eq_barista_random_bic, eq_barista_fixed3_bic)

eq_barista_comparison_fixed_bf <- eq_barista_fixed_bf / eq_barista_fixed2_bf  # BF for barista condition only vs. barista condition + EQ comparison

## Create model comparison table
### Random effects model comparison
eq_barista_random_anova_table <- eq_barista_random_anova_tidy %>%  
  mutate("Model specification" = c("Empty", "Participant", "Participant + EQ"),  # create specification column
         fixed = rep("-", 3),  # create fixed effects column
         random = c("1", "(1 | subject_nr)", "(1 | subject_nr) + (1 | EQ_mean)"), # create random effects column
         BF = NA) %>%  # create Bayes factor column
  select("Model specification", random, fixed, AIC, BIC, logLik, df, statistic, df, "p-value" = p.value, BF)  # select and rearrange columns

### Fixed effects model comparision
eq_barista_fixed_anova_table <- eq_barista_fixed_anova_tidy %>%  
  mutate("Model specification" = c("RE only", "Barista presence", "EQ", "Barista presence + EQ", "Barista presence * EQ"),  # create specification column
         fixed = c("-", "barista", "EQ_mean", "barista + EQ_mean", "barista * EQ_mean"),  # create fixed effects column
         random = c("(1 | subject_nr)", "(1 | subject_nr)", "(1 | subject_nr)", "(1 | subject_nr)", "(1 | subject_nr)"),  # create random effects column
         BF = c(NA, eq_barista_fixed2_bf, eq_barista_fixed3_bf, eq_barista_fixed_bf, eq_barista_full_bf)) %>%   # create Bayes factor column
  select("Model specification", random, fixed, AIC, BIC, logLik, df, statistic, df, "p-value" = p.value, BF)  # select and rearrange columns

### Combine tables
rand_row <- tibble("Model specification" = "Random effect models")  # create row for random effect models
fixed_row <- tibble("Model specification" = "Fixed effect models")  # create row for fixed effect models
eq_barista_anova_table <- bind_rows(rand_row, eq_barista_random_anova_table, fixed_row, eq_barista_fixed_anova_table)  # combine tables

## Check assumptions
boxplot(residuals(eq_barista_fixed2) ~ eq_data$barista)  # view distributions
leveneTest(residuals(eq_barista_fixed2) ~ eq_data$barista) # check assumption of equal variances
plot(eq_barista_fixed2)  # plot residuals vs. predicted values
plot(density(residuals(eq_barista_fixed2)))  # plot distribution of residuals
gg_qqplot(eq_barista_fixed2) # plot qqplot


# Descriptives ------------------------------------------------------------
# Participant gender descriptives
all_data %>% 
  count(gender)

# Participant age descriptives
mean(all_data$age, na.rm = TRUE)  # calculate participant mean age
sd(all_data$age, na.rm = TRUE) # calculate standard deviation for participant mean age

# Participant ethnicity descriptives
all_data %>% 
  count(race)

# Barista condition and payment type
## Barista condition
barista_condition_absent <- barista_condition %>% 
  filter(barista == "Barista absent")  # select barista absent condition
mean(barista_condition_absent$tip_mean)  # calculate mean tip amount for barista absent condition
sd(barista_condition_absent$tip_mean)  # calculate standard deviation for barista absent condition
length(unique(barista_condition_absent$subject_nr))  # find N

barista_condition_present <- barista_condition %>% 
  filter(barista == "Barista present")  # select barista present condition
mean(barista_condition_present$tip_mean)  # calculate mean tip amount for barista present condition
sd(barista_condition_present$tip_mean)  # calculate standard deviation for barista present condition
length(unique(barista_condition_present$subject_nr))  # find N

## Payment type
payment_condition_tipscreen <- payment_condition %>% 
  filter(payment_type == "Tip screen")  # select tip screen condition
mean(payment_condition_tipscreen$tip_mean)  # calculate mean tip amount for tip screen condition
sd(payment_condition_tipscreen$tip_mean)  # calculate standard deviation for tip screen condition
length(unique(payment_condition_tipscreen$subject_nr))  # find N

payment_condition_receipt <- payment_condition %>% 
  filter(payment_type == "Receipt")  # select receipt condition
mean(payment_condition_receipt$tip_mean)  # calculate mean tip amount for receipt condition
sd(payment_condition_receipt$tip_mean)  # calculate standard deviation for receipt condition
length(unique(payment_condition_receipt$subject_nr))  # find N

payment_condition_tj <- payment_condition %>% 
  filter(payment_type == "Cash")  # select cash condition
mean(payment_condition_tj$tip_mean)  # calculate mean tip amount for cash condition
sd(payment_condition_tj$tip_mean)  # calculate standard deviation for cash condition
length(unique(payment_condition_tj$subject_nr))  # find N

# Participant EQ descriptives
mean(eq_data$EQ_mean, na.rm = TRUE)  # calculate participant mean EQ
sd(eq_data$EQ_mean, na.rm = TRUE)  # calculate standard deviation for participant mean EQ
length(unique(eq_data$subject_nr))  # retrieve N for EQ data
