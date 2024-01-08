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
## Date Finalized: 2022-04-19
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
## Instructions: Create folders called "data" and "figures". Place this file in the main 
## 	directory. Place the data file in the data directory. Set the R working directory to 
##  the main directory.  At the R command prompt, type
## 	> source("goh_etal_2020_rcode.R")
## 	This will run the script, adding all of the calculated variables to the workspace and 
##  saving figures in the figures directory. If packages do not load properly, install them 
##  with install.packages("package_name").
## 
## Data file:
##  goh_etal_2020_data.csv
##   study = study number
##   subject_nr = participant number
##   gender = participant gender
##   age = participant age
##   race = participant ethnicity
##   condition = study 1: first condition experienced; study 2: condition experienced
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
library(afex)
library(BayesFactor)
library(car)
library(emmeans)
library(lme4)
library(lsr)
library(moments)
library(papaja)
library(patchwork)
library(tidyverse)


set.seed(0) # for consistency in model algorithm computations


# Import data -------------------------------------------------------------
all_data <- read_csv("goh_etal_2022_data.csv")  # import data
all_data_1 <- filter(all_data, study == 1) %>%   # separate study 1 data
  rename(first_condition = condition)  # rename condition column
all_data_2 <- filter(all_data, study == 2) %>%   # separate Study 2 data
  rename(condition_name = condition)  # rename condition column


# Study 1: Feelings towards tip screens and tip jars ----------------------
# Create dataset for feelings towards tipscreen and tipjar
feelings_data <- all_data_1 %>% 
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
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 30, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# t-test & Bayes analysis: compare participant feelings towards tip jars vs. tip screens
feelings_data_analyses <-  feelings_data %>%  # create dataset for analyses
  select(subject_nr, screen_jar, feelings_rating) %>%  # select columns
  pivot_wider(names_from = screen_jar, values_from = feelings_rating)  # pivot data to wide format

feelings_ttest_1 <- t.test(feelings_data_analyses$feel_tj, feelings_data_analyses$feel_ts, paired = TRUE)   # paired t-test: compare feelings towards tip jars vs. tip screens
feelings_bf_1 <- ttestBF(feelings_data_analyses$feel_tj, feelings_data_analyses$feel_ts, paired = TRUE)  # Bayes factor analysis
feelings_cohensd_1 <- cohensD(feelings_data_analyses$feel_tj, feelings_data_analyses$feel_ts)  # calculate effect size (Cohen's d) for feelings data


# Study 1: Avoidance of tip screens and tip jars --------------------------
# Create dataset for avoidance of tipscreen and tipjar
avoidance_data <- all_data_1 %>% 
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
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 30, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# t-test & Bayes analysis: compare participant avoidance of tip jars vs. tip screens
avoidance_data_analyses <-  avoidance_data %>%  # create dataset for analyses
  select(subject_nr, screen_jar, avoidance_frequency) %>%  # select columns
  pivot_wider(names_from = screen_jar, values_from = avoidance_frequency)  # pivot data to wide format

avoidance_ttest_1 <- t.test(avoidance_data_analyses$avoid_tj, avoidance_data_analyses$avoid_ts, paired = TRUE)   # paired t-test: compare avoidance of tip jars vs. tip screens
avoidance_bf_1 <- ttestBF(avoidance_data_analyses$avoid_tj, avoidance_data_analyses$avoid_ts, paired = TRUE)  # Bayes factor analysis
avoidance_cohensd_1 <- cohensD(avoidance_data_analyses$avoid_tj, avoidance_data_analyses$avoid_ts)  # calculate effect size (Cohen's d) for avoidance data

# Combine subfigures to one figure
(feeling_plot + avoidance_plot) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 25))
ggsave("figures/feeling_avoidance_payment_method.png", width = 13, height = 6.5)


# Study 1: Tip amounts ----------------------------------------------------
# Create dataset for participant tip amounts, barista conditions & payment type
tipping_data <- all_data_1 %>% 
  select(subject_nr, bp_ts:ba_tj, first_condition, EQ_mean) %>%   # select these columns from all_data_1 
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


# _Plots: Barista condition and payment type effects on tip amounts --------
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
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
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
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

barista_payment_plot_1 <- ggplot(tipping_data_clean, aes(x = barista, y = response, fill = barista)) +  # create barista absent/present & tip size graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  facet_wrap(~payment_type) +
  labs(x = "Barista presence", y = "Tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("Barista absent" = "Absent", "Barista present" = "Present")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#f02525", "#25f0f0")) +  # change boxplot fill colors
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph


# _Plots: Barista condition and payment type effects on tip amounts (first condition) -------
tipping_data_1_fc <- tipping_data_clean %>% 
  filter(first_condition == condition)  # select participants' first tipping condition 

# Create boxplot to compare participants' tip amounts by barista condition
barista_plot_first_condition <- ggplot(tipping_data_1_fc, aes(x = barista, y = response, fill = barista)) +  # create barista absent/present & tip size graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Barista presence", y = "Tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("Barista absent" = "Absent", "Barista present" = "Present")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#f02525", "#25f0f0")) +  # change boxplot fill colors
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# Create boxplot to compare participants' tip amounts by payment type
payment_plot_first_condition <- ggplot(tipping_data_1_fc, aes(x = payment_type, y = response, fill = payment_type)) +  # create payment type & tip size graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("cash" = "Cash", "receipt" = "Receipt", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# Combine within-subjects and first condition figures to one figure
payment_plot + barista_plot + payment_plot_first_condition + barista_plot_first_condition +
  plot_layout(heights = unit(c(10, 10), c("cm", "cm"))) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 20))
ggsave("figures/tipping_payment_method_barista_combined.png", width = 15, height = 12)


# _Check normality assumptions for tip amounts ----
# Calculate skewness & kurtosis of raw tip amounts
skewness(tipping_data_clean$response)  # calculate skewness of raw tip amounts
kurtosis(tipping_data_clean$response)  # calculate kurtosis of raw tip amounts

# Transform tip amounts to reduce skewness
tipping_data_sqrt <- tipping_data_clean %>%
  mutate(response_sqrt = sqrt(response))  # square root transform response

# Calculate skewness & kurtosis for square root transformed tip amounts
skewness(tipping_data_sqrt$response_sqrt, na.rm = TRUE)  # calculate skewness of square root transformed tip amounts
kurtosis(tipping_data_sqrt$response_sqrt, na.rm = TRUE)  # calculate kurtosis of square root transformed tip amounts


# _ANOVA analyses -------------------------------------------------
# Within-subjects ANOVA for barista presence and payment type effects on tip size
tipping_anova_1 <- aov_car(response_sqrt ~ barista * payment_type + Error(subject_nr/ (barista * payment_type)), data = tipping_data_sqrt)  # ANOVA
tipping_anova_1_table <- apa_print(tipping_anova_1)

# emm_barista_1 <- emmeans(tipping_anova_1, ~ barista)  # post-hoc contrasts: barista presence
# pairs(emm_barista_1, adjust = "tukey")  # pairwise comparisons
emm_payment_type_1 <- emmeans(tipping_anova_1, ~ payment_type)  # post-hoc contrasts: payment type
emm_payment_type_results_1 <- pairs(emm_payment_type_1, adjust = "tukey")  # pairwise comparisons

tipping_anova_1_bf <- anovaBF(response_sqrt ~ barista * payment_type, data = tipping_data_sqrt, whichRandom = "subject_nr")  # Bayes factor analyses

# First tipping condition ANOVA for barista presence and payment type effects on tip size
tipping_anova_1_fc <- aov_car(response ~ barista * payment_type + Error(subject_nr), data = tipping_data_1_fc)  # ANOVA
tipping_anova_1_fc_table <- apa_print(tipping_anova_1_fc)
tipping_anova_1_fc_bf <- anovaBF(response ~ barista * payment_type, data = tipping_data_1_fc)  # Bayes factor analyses


# Study 1: Empathy --------------------------------------------------------
# Prepare data
outliers_subjects <- outliers %>% 
  distinct(subject_nr)  # get only subject_nr of outliers

tipping_data_clean_eq <-  tipping_data_clean %>%  # create dataset to merge with EQ data
  select(subject_nr, condition:payment_type)

eq_data_1 <- all_data_1 %>% 
  filter(!is.na(EQ_mean), !is.na(bp_ts), !is.na(ba_ts), !is.na(bp_rec), !is.na(ba_rec), !is.na(bp_tj), !is.na(ba_tj)) %>%  # remove participants who do not have EQ_mean scores, tip amounts for the 6 tipping conditions
  select(subject_nr, bp_ts:ba_tj, EQ_mean) %>%  # select only these columns
  pivot_longer(bp_ts:ba_tj, names_to = "condition", values_to = "response") %>%
  left_join(tipping_data_clean_eq, by = c("subject_nr", "condition", "response")) %>%  # combine datasets
  anti_join(outliers_subjects) %>%  # remove all cases for outlier participants (i.e., not just outlier response values)
  mutate(response_sqrt = sqrt(response)) %>%  # square root transform response
  select(-condition, -response) %>%   # remove condition and response columns
  group_by(subject_nr, EQ_mean, barista) %>% 
  summarise(response_mean = mean(response_sqrt, na.rm = TRUE))  # calculate participants' mean tip amounts for each barista (absent/present) condition


# _Plot: Empathy effects on barista presence ------------------------------
eq_data_plot <- eq_data_1 %>%
  left_join(barista_condition, by = c("subject_nr", "barista")) # combine eq_data with barista_condition

# Create scatterplot to compare EQ effects on participant mean tip amounts by barista condition
empathy_plot <- ggplot(eq_data_plot, aes(x = EQ_mean, y = tip_mean)) +  # create mean EQ scores & mean tip amount axes labels
  geom_jitter(width = 0.01, height = 0.01) +
  geom_smooth(method = "lm") +
  facet_wrap(~barista) +
  labs(x = "Mean empathy score", y = "Mean tip size ($)") +  # rename x & y axes labels
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        strip.text.x = element_text(size = 28),  # change facet label size
        legend.position = "none")  # remove legend from graph


# _Plot: Empathy effects on barista presence (first condition) ------------------------------
# Create scatterplot to compare EQ effects on participant tip amounts by barista condition
empathy_plot_first_condition <- ggplot(tipping_data_1_fc, aes(x = EQ_mean, y = response)) +  # create mean EQ scores & tip amount axes labels
  geom_jitter(width = 0.01, height = 0.01) +
  geom_smooth(method = "lm") +
  facet_wrap(~barista) +
  labs(x = "Mean empathy score", y = "Tip size ($)") +  # rename x & y axes labels
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        strip.text.x = element_text(size = 28),  # change facet label size
        legend.position = "none")  # remove legend from graph

# Combine within-subjects and first condition figures to one figure
(empathy_plot + empathy_plot_first_condition) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 25))
ggsave("figures/tipping_empathy_barista_combined.png", width = 15, height = 6.5)


# _Regression analyses ----------------------------------------------------
eq_lmer_1 <- lmer(response_mean ~ barista * EQ_mean + (1 | subject_nr), data = eq_data_1)
eq_lmer_1_table <- apa_print(eq_lmer_1)
eq_lmer_1_gentestbf <- generalTestBF(response_mean ~ barista * EQ_mean, data = eq_data_1, whichRandom = "subject_nr")  # Bayes factor analyses
eq_lmer_1_bf <- eq_lmer_1_gentestbf[4]/eq_lmer_1_gentestbf[3]

# First tipping condition ANOVA for effect of empathy on barista presence on tip size
eq_lm_1_fc <- lm(response ~ barista * EQ_mean, data = tipping_data_1_fc)
eq_lm_1_fc_table <- apa_print(eq_lm_1_fc)
eq_lm_1_fc_bf <- lmBF(response ~ barista * EQ_mean, data = tipping_data_1_fc) / lmBF(response ~ barista + EQ_mean, data = tipping_data_1_fc)  # Bayes factor analyses


# Study 1: Descriptives ---------------------------------------------------
# Participant gender descriptives
gender_1 <- all_data_1 %>% 
  count(gender)

# Participant age descriptives
age_1 <- mean(all_data_1$age, na.rm = TRUE)  # calculate participant mean age
age_sd_1 <- sd(all_data_1$age, na.rm = TRUE) # calculate standard deviation for participant mean age

# Participant ethnicity descriptives
ethnicity_1 <- all_data_1 %>% 
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

payment_condition_cash <- payment_condition %>% 
  filter(payment_type == "Cash")  # select cash condition
mean(payment_condition_cash$tip_mean)  # calculate mean tip amount for cash condition
sd(payment_condition_cash$tip_mean)  # calculate standard deviation for cash condition
length(unique(payment_condition_cash$subject_nr))  # find N

# Participant EQ descriptives
mean(eq_data_1$EQ_mean, na.rm = TRUE)  # calculate participant mean EQ
sd(eq_data_1$EQ_mean, na.rm = TRUE)  # calculate standard deviation for participant mean EQ
length(unique(eq_data_1$subject_nr))  # retrieve N for EQ data


# Study 2: Tip amounts ----------------------------------------------------
# Create dataset for participant tip amounts, barista conditions & payment type
tipping_data_2 <- all_data_2 %>% 
  pivot_longer(bp_ts:ba_tj, names_to = "condition", values_to = "response") %>%   # sort participant responses by barista condition & payment condition
  filter(condition_name == condition) %>%   # drop duplicate study condition rows
  mutate(barista = ifelse(grepl("bp", condition), "Barista present",  # insert "Barista present" label in "barista" column if tipping condition contains "bp" in "condition" column
                          ifelse(grepl("ba", condition), "Barista absent", NA)),  # insert "Barista absent" label in "barista" column if tipping condition contains "ba" in "condition" column
         payment_type = ifelse(grepl("ts", condition), "Tip screen",  # insert "Tip screen" in "payment_type" column if tipping condition contains "ts" in "condition" column
                               ifelse(grepl("rec", condition), "Receipt",  # insert "Receipt" in "payment_type" column if tipping condition contains "rec" in "condition" column
                                      ifelse(grepl("tj", condition), "Cash", NA))),  # insert "Cash" in "payment_type" column if tipping condition contains "tj" in "condition" column
         response = Recode(response, "NA = 0"),  # recode NA responses in "response" column to 0
         condition = as.factor(condition),   # convert variables to factor type
         barista = as.factor(barista),
         payment_type = as.factor(payment_type))

# Calculate summary statistics
tip_size_mean_2 <- mean(tipping_data_2$response)  # calculate mean tip size
tip_size_sd_2 <- sd(tipping_data_2$response)  # calculate standard deviation for tip size

# Identify outliers (cut-off value based on exclusion criteria from Study 1 that excluded tip amounts >= $2)
outliers_2 <- tipping_data_2 %>%
  filter(response >= 2)

tipping_data_clean_2 <- tipping_data_2 %>% 
  anti_join(outliers_2)  # remove outliers


# _Plots: Barista condition and payment type effects on tip amounts --------
# Plot barista condition effect on tip amounts
barista_plot_2 <- ggplot(tipping_data_clean_2, aes(x = barista, y = response, fill = barista)) +  # create barista absent/present & tip size graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Barista presence", y = "Tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("Barista absent" = "Absent", "Barista present" = "Present")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#f02525", "#25f0f0")) +  # change boxplot fill colors
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# Plot payment type effect on tip amounts
payment_plot_2 <- ggplot(tipping_data_clean_2, aes(x = payment_type, y = response, fill = payment_type)) +  # create payment type & tip size graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("cash" = "Cash", "receipt" = "Receipt", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# Plot barista condition and payment type interaction effect on tip amounts
barista_payment_plot_2 <- ggplot(tipping_data_clean_2, aes(x = barista, y = response, fill = barista)) +  # create barista absent/present & tip size graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  facet_wrap(~payment_type) +
  labs(x = "Barista presence", y = "Tip size ($)") +  # rename x & y axes labels
  scale_x_discrete(labels = c("Barista absent" = "Absent", "Barista present" = "Present")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#f02525", "#25f0f0")) +  # change boxplot fill colors
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 24, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 20),  # change x & y axes tick mark font & size
        strip.text.x = element_text(size = 22),  # change facet label size
        legend.position = "none")  # remove legend from graph

# Combine subfigures to one figure
arrange_plots <- c(
  area(1, 1, 1, 2),
  area(1, 3, 1, 5),
  area(2, 2, 2, 4)
)
payment_plot_2 + barista_plot_2 + barista_payment_plot_2 + plot_layout(design = arrange_plots) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 20))
ggsave("figures/tipping_payment_method_barista_2.png", width = 15, height = 12)
# (payment_plot_2 + barista_plot_2) / barista_payment_plot_2 +
#   plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
#   theme(plot.tag = element_text(size = 25))


# _Check normality assumptions for tip amounts ----
# Calculate skewness & kurtosis of raw tip amounts
skewness(tipping_data_clean_2$response)  # calculate skewness of raw tip amounts
kurtosis(tipping_data_clean_2$response)  # calculate kurtosis of raw tip amounts


# _ANOVA analyses -------------------------------------------------
# Between subjects ANOVA for barista presence and payment type effects on tip size
tipping_anova_2 <- aov_car(response ~ barista * payment_type + Error(subject_nr), data = tipping_data_clean_2)  # ANOVA
tipping_anova_2_table <- apa_print(tipping_anova_2)

# emm_barista_2 <- emmeans(tipping_anova_2, ~ barista)  # post-hoc contrasts: barista presence
# pairs(emm_barista_2, adjust = "tukey")  # pairwise comparisons
# emm_payment_type_2<- emmeans(tipping_anova_2, ~ payment_type)  # post-hoc contrasts: payment type
# pairs(emm_payment_type_2, adjust = "tukey")  # pairwise comparisons
emm_interaction_2<- emmeans(tipping_anova_2, ~ barista | payment_type)  # post-hoc contrasts: interaction
emm_interaction_results_2 <- pairs(emm_interaction_2, adjust = "tukey")  # pairwise comparisons

tipping_anova_2_bf <- anovaBF(response ~ barista * payment_type, data = tipping_data_clean_2)  # Bayes factor analyses

res <- residuals(tipping_anova_2)
hist(res, main = "Histogram of residuals", xlab = "Residuals")
lines(density(res),col=2)
plot(emm_interaction_2, comparisons = TRUE)
plot(density(residuals(tipping_anova_2)), ask = FALSE)  # plot distribution of residuals


# Study 2: Empathy --------------------------------------------------------
# _Plot: Empathy effects on barista presence ------------------------------
ggplot(tipping_data_clean_2, aes(x = EQ_mean, y = response)) +  # create mean EQ scores & tip amount axes labels
  geom_jitter(width = 0.01, height = 0.01) +
  geom_smooth(method = "lm") +
  facet_wrap(~barista) +
  labs(x = "Mean empathy score", y = "Tip size ($)") +  # rename x & y axes labels
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 32, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 32, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        strip.text.x = element_text(size = 30),  # change facet label size
        legend.position = "none")  # remove legend from graph
ggsave("figures/tipping_empathy_barista_2.png", width = 7.5, height = 6.5)  # save plot
# ggsave("figures/tipping_empathy_barista_2.png", width = 15, height = 6.5)  # save plot for Markdown


# _Regression analyses ----------------------------------------------------
eq_lm_2 <- lm(response ~ barista * EQ_mean, data = tipping_data_clean_2)
eq_lm_2_table <- apa_print(eq_lm_2)
eq_lm_2_bf <- lmBF(response ~ barista * EQ_mean, data = tipping_data_clean_2) / lmBF(response ~ barista + EQ_mean, data = tipping_data_clean_2)  # Bayes factor analyses


# Study 2: Feelings towards tip screens and tip jars ----------------------
# Create dataset for feelings towards tipscreen and tipjar
feelings_data_2 <- tipping_data_clean_2 %>% 
  filter(!is.na(feel_ts), !is.na(feel_tj)) %>%  # remove participants who have NA responses for feeling questions
  select(subject_nr, feel_ts, feel_tj) %>%  # select these columns from all_data
  gather(feel_ts, feel_tj, key = "screen_jar", value = "feelings_rating") %>%  # sort participant responses by payment type (tipscreen/tipjar)
  mutate(payment_type = ifelse(grepl("ts", screen_jar), "tipscreen",  # insert "tipscreen" label in "payment_type" column if "screen_jar" column contains "ts"
                               ifelse(grepl("tj", screen_jar), "tipjar", NA)))  # insert "tipjar" label in "payment_type" column if "screen_jar" column contains "tj"

# Create boxplot to compare participants' feelings towards tip screens & tip jars (higher ratings = more negative feelings toward payment type)
feeling_plot_2 <- ggplot(feelings_data_2, aes(x = payment_type, y = feelings_rating, fill = payment_type)) +  # create payment type & feelings rating graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Degree of negative feelings") +  # rename x & y axes labels
  scale_x_discrete(labels = c("tipjar" = "Tip jar", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  scale_y_continuous(limits = c(1, 6)) +  # specify y-axis limits
  scale_fill_manual(values=c("#0C7BDC", "#FFC20A")) +  # change boxplot fill colors
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 30, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# t-test & Bayes analysis: compare participant feelings towards tip jars vs. tip screens
feelings_data_analyses_2 <-  feelings_data_2 %>%  # create dataset for analyses
  select(subject_nr, screen_jar, feelings_rating) %>%  # select columns
  pivot_wider(names_from = screen_jar, values_from = feelings_rating)  # pivot data to wide format

feelings_ttest_2 <- t.test(feelings_data_analyses_2$feel_tj, feelings_data_analyses_2$feel_ts, paired = TRUE)   # paired t-test: compare feelings towards tip jars vs. tip screens
feelings_bf_2 <- ttestBF(feelings_data_analyses_2$feel_tj, feelings_data_analyses_2$feel_ts, paired = TRUE)  # Bayes factor analysis
feelings_cohensd_2 <- cohensD(feelings_data_analyses_2$feel_tj, feelings_data_analyses_2$feel_ts)  # calculate effect size (Cohen's d) for feelings data


# Study 2: Avoidance of tip screens and tip jars --------------------------
# Create dataset for avoidance of tipscreen and tipjar
avoidance_data_2 <- tipping_data_clean_2 %>% 
  filter(!is.na(avoid_ts), !is.na(avoid_tj)) %>%  # remove participants who have NA responses for avoidance questions
  select(subject_nr, avoid_ts, avoid_tj) %>%  # select these columns from all_data
  gather(avoid_ts, avoid_tj, key = "screen_jar", value = "avoidance_frequency") %>%  # sort participant responses by payment type (tip screen/tip jar)
  mutate(payment_type = ifelse(grepl("ts", screen_jar), "tipscreen",  # insert "tipscreen" label in "payment_type" column if "screen_jar" column contains "ts"
                               ifelse(grepl("tj", screen_jar), "tipjar", NA)))  # insert "tipjar" label in "payment_type" column if "screen_jar" column contains "tj"

# Create boxplot to compare participants' frequency of avoidance of tipscreens & tipjars (Qualtrics response options key: 1 = Never, 2 = Once, 3 = 2-5 times, 4 = 6-10 times, 5 = more than 10 times)
avoidance_plot_2 <- ggplot(avoidance_data_2, aes(x = payment_type, y = avoidance_frequency, fill = payment_type)) +  # create payment type & frequency of avoidance of tipscreen/tipjar graph axes labels
  geom_boxplot(outlier.shape = NA) +  # create boxplot & remove outlier points
  stat_summary(fun.data = "mean_cl_boot", size = 1.25) +  # generate mean and confidence level bars
  labs(x = "Payment method", y = "Frequency of avoidance") +  # rename x & y axes labels
  scale_x_discrete(labels = c("tipjar" = "Tip jar", "tipscreen" = "Tip screen")) +  # rename x-axis tick mark labels
  scale_fill_manual(values=c("#0C7BDC", "#FFC20A")) +  # change boxplot fill colors
  theme_classic(base_family = "Arial") +
  theme(axis.title.y = element_text(face = "bold", size = 30, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # change y axis title font, size & margin
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t = 10, r = 0, b = 0, l = 0)),  # change x axis title font, size & margin
        axis.text = element_text(face = "bold", size = 26),  # change x & y axes tick mark font & size
        legend.position = "none")  # remove legend from graph

# t-test & Bayes analysis: compare participant avoidance of tip jars vs. tip screens
avoidance_data_analyses_2 <-  avoidance_data_2 %>%  # create dataset for analyses
  select(subject_nr, screen_jar, avoidance_frequency) %>%  # select columns
  pivot_wider(names_from = screen_jar, values_from = avoidance_frequency)  # pivot data to wide format

avoidance_ttest_2 <- t.test(avoidance_data_analyses_2$avoid_tj, avoidance_data_analyses_2$avoid_ts, paired = TRUE)   # paired t-test: compare avoidance of tip jars vs. tip screens
avoidance_bf_2 <- ttestBF(avoidance_data_analyses_2$avoid_tj, avoidance_data_analyses_2$avoid_ts, paired = TRUE)  # Bayes factor analysis
avoidance_cohensd_2 <- cohensD(avoidance_data_analyses_2$avoid_tj, avoidance_data_analyses_2$avoid_ts)  # calculate effect size (Cohen's d) for avoidance data

# Combine subfigures to one figure
(feeling_plot_2 + avoidance_plot_2) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ')') &
  theme(plot.tag = element_text(size = 25))
ggsave("figures/feeling_avoidance_payment_method_2.png", width = 13, height = 6.5)


# Study 2: Descriptives ---------------------------------------------------
# Participant gender descriptives
gender_2 <- tipping_data_clean_2 %>% 
  count(gender)

# Participant age descriptives
age_2 <- mean(tipping_data_clean_2$age, na.rm = TRUE)  # calculate participant mean age
age_sd_2 <- sd(tipping_data_clean_2$age, na.rm = TRUE) # calculate standard deviation for participant mean age

# Participant ethnicity descriptives
ethnicity_2 <- tipping_data_clean_2 %>% 
  count(race)

# Barista condition and payment type
## Barista condition
barista_condition_absent_2 <- tipping_data_clean_2 %>% 
  filter(barista == "Barista absent")  # select barista absent condition
mean(barista_condition_absent_2$response)  # calculate mean tip amount for barista absent condition
sd(barista_condition_absent_2$response)  # calculate standard deviation for barista absent condition
length(unique(barista_condition_absent_2$subject_nr))  # find N

barista_condition_present_2 <- tipping_data_clean_2 %>% 
  filter(barista == "Barista present")  # select barista present condition
mean(barista_condition_present_2$response)  # calculate mean tip amount for barista present condition
sd(barista_condition_present_2$response)  # calculate standard deviation for barista present condition
length(unique(barista_condition_present_2$subject_nr))  # find N

## Payment type
payment_condition_tipscreen_2 <- tipping_data_clean_2 %>% 
  filter(payment_type == "Tip screen")  # select tip screen condition
mean(payment_condition_tipscreen_2$response)  # calculate mean tip amount for tip screen condition
sd(payment_condition_tipscreen_2$response)  # calculate standard deviation for tip screen condition
length(unique(payment_condition_tipscreen_2$subject_nr))  # find N

payment_condition_receipt_2 <- tipping_data_clean_2 %>% 
  filter(payment_type == "Receipt")  # select receipt condition
mean(payment_condition_receipt_2$response)  # calculate mean tip amount for receipt condition
sd(payment_condition_receipt_2$response)  # calculate standard deviation for receipt condition
length(unique(payment_condition_receipt_2$subject_nr))  # find N

payment_condition_cash_2 <- tipping_data_clean_2 %>% 
  filter(payment_type == "Cash")  # select cash condition
mean(payment_condition_cash_2$response)  # calculate mean tip amount for cash condition
sd(payment_condition_cash_2$response)  # calculate standard deviation for cash condition
length(unique(payment_condition_cash_2$subject_nr))  # find N

# Participant EQ descriptives
mean(tipping_data_clean_2$EQ_mean, na.rm = TRUE)  # calculate participant mean EQ
sd(tipping_data_clean_2$EQ_mean, na.rm = TRUE)  # calculate standard deviation for participant mean EQ
length(unique(tipping_data_clean_2$subject_nr))  # retrieve N for EQ data

save.image("tipping_workspace.RData")
