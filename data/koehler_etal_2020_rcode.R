## ---
##
## Script name: koehler_etal_2020_rcode.R
##
## Purpose of script: This script computes the publishable analyses for our study examining the effects of exercise on food intertemporal choice.
##
## Authors: Elise Thayer (eliserthayer@gmail.com), Jeffrey R. Stevens  (jeffrey.r.stevens@gmail.com)
##
## Date Created: 2020-05-31
##
## Date Finalized: 2020-12-11
##
## License: All materials presented here are released under the Creative Commons Attribution 4.0 International Public License (CC BY 4.0).
##  You are free to:
##   Share — copy and redistribute the material in any medium or format
##   Adapt — remix, transform, and build upon the material for any purpose, even commercially.
##  Under the following terms:
##   Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
##   No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
##
## ---
##
## Notes:
## Instructions: Place this file and the data files in the main directory.
## 	Set the R working directory to the main directory. At the R command prompt, type
## 	> source("koehler_etal_2020_rcode.R")
## 	This will run the script, adding all of the calculated variables to the workspace and
##  calculating the inferential statistics. If packages do not load properly, install them
##  with install.packages("package_name").
## Data files:
##  koehler_etal_2020_data.csv
##   datatype - Flag for type of data (Amount or Choice)
##   participant - participant ID#
##   condition - Experimental condition (exercise or rest)
##   time - Timeframe for decision (now or later)
##   diff - Difference score type (prepost = post - pre, prepost30 = post+30 - pre)
##   mean - Mean difference score
##
## ---


# Load libraries ---------------------------------------------------------------

library(afex) # for interaction plotslibrary(lsr) # for Cohen's d
library(ggpubr)  # for ggqqplot
library(lsr) # for cohens d
library(papaja) # for within-subjects confidence intervals (available at https://github.com/crsh/papaja)
library(tidyverse) # tidyverse


# Define functions ----------------------------------------------------------

## Analyzes amount ("food amount preference") data after collapsing over food
analyze_amount <- function(df, dv) {
  # Filter dataframe by measure of interest
  dv_df <- df %>%
    filter(diff == dv)

  # Calculate ANOVA using afex::aov_ez
  dv_aov_ez <-  aov_ez("participant", "mean", data = dv_df, within = c("condition", "time"), between = NULL, anova_table = list(es = "ges"))

  # Text model assumptions
  ## Plot qqplot for assumption check
  print(residuals_qqplot(dv_aov_ez))

  ## Plot density plot over condition and time for assumption check
  print(ggplot(dv_df, aes(mean, fill = time)) +
          geom_density(alpha = 0.2) +  # plot density plot
          facet_grid(.~ condition) +  # facet by condition
          ggtitle(dv) +  # print title
          theme_classic())  # use classic theme

  # Print ANOVA output to console
  print("#### ANOVA summary ####")  # print to screen
  print(summary(dv_aov_ez))  # print ANOVA results
  print(knitr::kable(nice(dv_aov_ez)))  # print ANOVA results

  # Plot afex interaction plot
  print(afex_plot(dv_aov_ez, x = "time", trace = "condition",  # interaction plot
                  error = "within",
                  mapping = c("color", "fill"),
                  data_geom = geom_boxplot, data_arg = list(width = 0.4),
                  point_arg = list(size = 1.5), line_arg = list(size = 1)) +
          theme_classic())

  # Generate means and within-subject confidence intervals
  print("#### Means and CIs ####")   # print to screen
  ## Condition
  ### Calculate means per participant and condition (over time)
  dv_condition <- dv_df %>%
    group_by(participant, condition) %>%
    summarise(mean = mean(mean), .groups = "keep")
  ### Calculate means per condition (over participant and time)
  dv_condition_cis <- dv_condition %>%
    group_by(condition) %>%
    summarise(mean = mean(mean), .groups = "keep")
  ### Calculate within-subjects 95% confidence intervals
  dv_condition_cis$wsci <- wsci(data = dv_condition, id = "participant", factors = "condition", dv = "mean")$mean  # calculate within-subject CIs
  dv_condition_cis <- mutate(dv_condition_cis, uci = mean + wsci, lci = mean - wsci)  # create columns of upper and lower bounds
  print(dv_condition_cis)

  ## Time
  ### Calculate means per participant and time (over condition)
  dv_time <- dv_df %>% # prep
    group_by(participant, time) %>%
    summarise(mean = mean(mean), .groups = "keep")
  ### Calculate means per time (over participant and condition)
  dv_time_cis <- dv_time %>% # prep means
    group_by(time) %>%
    summarise(mean = mean(mean), .groups = "keep")
  ### Calculate within-subjects 95% confidence intervals
  dv_time_cis$wsci <- wsci(data = dv_time, id = "participant", factors = "time", dv = "mean")$mean  # calculate within-subject CIs
  dv_time_cis <- mutate(dv_time_cis, uci = mean + wsci, lci = mean - wsci)  # create columns of upper and lower bounds
  print(dv_time_cis)

}

## Analyzes choice ("intertemporal food preference") data after filtering by only same food
analyze_choice <- function(df, dv) {
  # Filter data by measure of interest
  dv_df <- df %>%
    filter(diff == dv)
  dv_df_wide <- pivot_wider(dv_df, names_from = "condition", values_from = "mean") %>%  # convert to wide form
    mutate(diff = exercise - rest)  # calculate column of difference scores

  # Print t-test output
  print("#### T-test summary ####")  # print to screen
  dv_ttest <- t.test(mean ~ condition, data = dv_df, paired = TRUE)  # calculate t-test
  print(dv_ttest)

  # Check assumptions
  print("#### Check asumptions ####")  # print to screen
  ## Calculate Shapiro test of normality of differences
  dv_shapiro <- with(dv_df, mean[condition == "exercise"] - mean[condition == "rest"])  # compute differences
  print(shapiro.test(dv_shapiro))  # test normality of differences

  ## Plot qqplot of difference scores
  print(ggqqplot(dv_df_wide, x = "diff"))

  # Effect size (Cohen's d)
  print("#### Effect size (Cohen's d) ####")  # print to screen
  print(cohensD(mean ~ condition, data = dv_df, method = "paired"))

  # Means and confidence intervals
  ## Calculate means per condition (over participants)
  dv_cis <- dv_df %>%  # prep means
    group_by(condition) %>%
    summarise(mean = mean(mean), .groups = "keep")

  ## Calculate within-subjects 95% confidence intervals
  dv_cis$wsci <- wsci(data = dv_df, id = "participant", factors = "condition", dv = "mean")$mean  # calculate within-subject CIs
  dv_cis <- mutate(dv_cis, uci = mean + wsci, lci = mean - wsci)  # create columns of upper and lower bounds
  print(dv_cis)
}

# Calculates residuals from afex_aov model (from https://github.com/singmann/afex/issues/64)
residuals.afex_aov <- function(object, model = "multivariate", ...) {
  if (length(attr(object, "within")) == 0 || model == "multivariate") {
    return(residuals(object$lm))
  } else {
    data <- object$data$long
    dv <- attr(object, "dv")
    id <- attr(object, "id")
    between <- names(attr(object, "between"))
    within <- names(attr(object, "within"))

    # within
    combs <- expand.grid(lapply(within, function(x) c(x, NA)))
    combs$id <- id
    combs <- head(combs, -1)
    within_res <- list()
    for (i in seq_len(nrow(combs))) {
      tem_fs <- as.vector(na.omit(t(combs[i, ])))
      ag_data <- aggregate(data[, dv], data [, tem_fs], mean)
      temp_name <- paste0(head(tem_fs, -1), collapse = "*")
      form <- formula(paste0("x~", temp_name))
      within_res[[temp_name]] <- residuals(lm(form, ag_data))
    }

    all_residuals <- within_res

    # between
    if (!is.null(between)) {
      ag_data <- aggregate(data[, dv], data[, c(between, id)], mean)
      form <- formula(paste0("x~", paste0(c(between), collapse = "*")))
      all_residuals[[id]] <- residuals(lm(form, ag_data))
    }
    return(all_residuals)
  }
}

# Plots qqplot from afex_aov model  (adapted from https://github.com/singmann/afex/issues/64)
residuals_qqplot <- function(object) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package ggplot2 is required.", call. = FALSE)
  }
  all_residuals <- residuals(object, model = "univariate")

  if (is.list(all_residuals)) {
    all_residuals <- lapply(names(all_residuals), function(x) data.frame(residuals = all_residuals[[x]], proj = x))

    plot_data <- do.call("rbind", all_residuals)
  } else {
    plot_data <- data.frame(residuals = all_residuals,
                            proj = "Error")
  }

  ggpubr::ggqqplot(plot_data, x = "residuals", facet.by = "proj")

}

# Import and prepare data ----------------------------------------------------
all_data <- read_csv("koehler_etal_2020_data.csv")  # import data
amount_data <- filter(all_data, datatype == "Amount")  # filter amount data
choice_data <- filter(all_data, datatype == "Choice")  # filter choice data


# Analyze data ----------------------------------------------------
# _Amount "Amount food preference" ----
# Post - pre
print("############ Amount post - pre ############")  # print to screen
analyze_amount(amount_data, "prepost")

# Post + 30 - pre
print("############ Amount post+30 - pre ############")  # print to screen
analyze_amount(amount_data, "prepost30")

# _Choice "Intertemporal food preferences" ----
# Post - pre
print("############ Choice post - pre ############")  # print to screen
analyze_choice(choice_data, "prepost")

# Post + 30 - Pre
print("############ Choice post+30 - pre ############")  # print to screen
analyze_choice(choice_data, "prepost30")
