###################################################
### stevens_etal_2018_rcode.R
### Created by Jeffrey R. Stevens on 07 Aug 2015 (jeffrey.r.stevens@gmail.com),
###  Current version: 30 May 2018
### Citation: Stevens, Stevens JR, Woike JK, Schooler LJ, Lindner S, Pachur T. 2018 
###   Social contact patterns can buffer costs of forgetting in the evolution of 
###   cooperation. Proceedings of the Royal Society of London: Series B, 20180407. 
###   http://dx.doi.org/10.1098/rspb.2018.0407
### Summary: This script generates plots of memory functions and effects of
###   skew and forgetting on cooperation rates and strategy evolution
### Instructions: Place this file and the data file (stevens_etal_2018_data.csv)
### 	in the same directory. Create a directory called 'figures'.
###   Set the R working directory to the directory with this R script.  
###   At the R command prompt, type
### 	> source("stevens_etal_2018_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PNG versions of the figures in the figures directory.
### Uses: This script can be reproduced and modified for personal and scientific use.
### Data files:
###  stevens_etal_2018_data--data for agent-based, evolutionary simulation of repeated
###		prisoner's dilemma games
###		measure - describes the type of measure: Cooperation (proportion of interactions
###			with cooperation in last generation), strategies (ALLC, ALLD, etc.: proportion
###			of agents in the population in last generation that play that strategy),
###			Robustness (proportion of interactions with cooperation in the last generation
###			for the additional simultions checking the robustness to interaction number)
###		lambda - starting point parameter for Wickelgren's (1974) forgetting function
###		psi - decay parameter for Wickelgren's (1974) forgetting function
###		memory_error_type - type of memory error (guessing error, skipping error,
###			or perfect memory with no error)
###		skew - degree of skew in the distribution of contact frequency (high, low, no skew)
###		interactions - number of interactions experienced by agents (100, 300, 600, or 1100)
###		mean - mean proportion of cooperation or mean proportion of agents in population
###		sd - standard deviation of proportion of cooperation or proportion of agents in population
###################################################

################
### Clear and define variables, load libraries
################
options(warn = -1)  # disable warnings

suppressPackageStartupMessages(library(car))      # needed for Recode
library(ggplot2)  # needed for plots

col.blind <- c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#F0E442", "#CC79A7", "#56B4E9", "black")   # create vector of color-blind friendly colors for figures

total_interactions <- round(39 * 4500 * 1000 * 250 / 1000000000, 1)	# calculate total number of interactions in simulations (excluding robustness simulations)

################
### Prepare simulation data
################
sim_data <- read.csv("stevens_etal_2018_data.csv") # read in data file
sim_data$percent <- sim_data$mean * 100       # convert proportions to percentages
sim_data$quant_skew <- as.numeric(as.character(sim_data$skew)) # copy skew to numeric quantitative skew column
sim_data$skew <- factor(sim_data$skew, levels = c("No skew", "Low skew", "High skew"))  # reorder skew levels
sim_data$forget_rate <- as.factor(paste(paste("λ=", round(sim_data$lambda, 2), sep = ""), paste("Ψ=", round(sim_data$psi, 2), sep = ""), sep = "\n"))  # create forget_rate column
sim_data$memory_condition <- paste(sim_data$memory_error_type, sim_data$forget_rate, sep = "\n") # create memory_condition column
sim_data$forget_skew <- as.factor(paste(sim_data$memory_error_type, sim_data$skew, sep = ", ")) # create forget_skew column
sim_data$forget_rate <- factor(sim_data$forget_rate, levels = c("λ=1\nΨ=0", "λ=1\nΨ=0.1", "λ=1\nΨ=0.25", "λ=0.87\nΨ=0.22", "λ=0.6\nΨ=0.25", "λ=0.4\nΨ=0.25", "λ=0.25\nΨ=0.5"))  # reorder forget_rate levels
sim_data <- sim_data[order(sim_data$forget_rate), ] # reorder data frame by forget_rate

################
### Create cooperation rate plots
################

########
### Prepare cooperation data
########
coop_data <- subset(sim_data, measure == "Cooperation") # subset cooperation data
coop_data$forget_skew <- factor(coop_data$forget_skew)  # remove unused factor levels
coop_forget <- aggregate(mean ~ memory_error_type, data = coop_data, FUN = "mean")  # calculate mean cooperation rate for each memory eror type (Skipping or guessing)

# Calculate differences in cooperation rates between memory conditions
coop_guessing_min <- 100 - round(coop_data$mean[coop_data$memory_condition == "Guessing\nλ=1\nΨ=0.25" & coop_data$skew == "No skew"] * 100, 0)
coop_guessing_max <- 100 - round(coop_data$mean[coop_data$memory_condition == "Guessing\nλ=1\nΨ=0.25" & coop_data$skew == "High skew"] * 100, 0)
coop_guessing_diff_max <- round((coop_data$mean[coop_data$memory_condition == "Guessing\nλ=1\nΨ=0.25" & coop_data$skew == "High skew"] - coop_data$mean[coop_data$memory_condition == "Guessing\nλ=1\nΨ=0.25" & coop_data$skew == "No skew"]) * 100, 0)
coop_guessing_diff_emp <- round((coop_data$mean[coop_data$memory_condition == "Guessing\nλ=0.87\nΨ=0.22" & coop_data$skew == "High skew"] - coop_data$mean[coop_data$memory_condition == "Guessing\nλ=0.87\nΨ=0.22" & coop_data$skew == "No skew"]) * 100, 0)
coop_skipping_diff_max <- round((coop_data$mean[coop_data$memory_condition == "Skipping\nλ=0.4\nΨ=0.25" & coop_data$skew == "High skew"] - coop_data$mean[coop_data$memory_condition == "Skipping\nλ=0.4\nΨ=0.25" & coop_data$skew == "No skew"]) * 100, 0)
coop_skipping_diff_emp <- round((coop_data$mean[coop_data$memory_condition == "Skipping\nλ=0.87\nΨ=0.22" & coop_data$skew == "High skew"] - coop_data$mean[coop_data$memory_condition == "Skipping\nλ=0.87\nΨ=0.22" & coop_data$skew == "No skew"]) * 100, 0)

## Robustness data
coop_robust_means <- subset(sim_data, (measure == "Cooperation" & lambda == 0.87 & psi == 0.22 & memory_error_type == "Skipping" & interactions == 100) | measure == "Robustness")	# subset data from robustness simulations
coop_robust_means <- coop_robust_means[order(coop_robust_means$interactions), ] # reorder data frame by interactions

## Matrix 2 data
coop2_data <- subset(coop_data, psi != 0.1 & psi != 0.5 & skew != "Low skew") # subset conditions used for matrix 2
matrix2_data <- subset(sim_data, measure == "Matrix2")  # subset matrix 2 data
coop_matrix2_data <- rbind(coop2_data, matrix2_data)  # combine matrix 1 and 2 data
coop_matrix2_data <- subset(coop_matrix2_data, memory_error_type != "Perfect")  # subset conditions with memory errors
coop_matrix2_data$matrix <- Recode(coop_matrix2_data$measure, "'Cooperation'='Original matrix';'Matrix2'='New matrix'") # recode matrix names
coop_matrix2_data$matrix <- factor(coop_matrix2_data$matrix, levels = c("Original matrix", "New matrix")) # reorder matrix labels
coop_matrix2_data$memory_error_type <- factor(coop_matrix2_data$memory_error_type)	# remove unused levels
coop_matrix2_data$skew <- factor(coop_matrix2_data$skew)	# remove unused levels

# Operator data
coop3_data <- subset(coop_data, (lambda == 0.87 | lambda == 0.6) & skew != "Low skew")  # subset conditions used for different selection operators
coop3_data$measure <- "Roulette-wheel"  # rename measure
operator_data <- subset(sim_data, measure == "Truncation" | measure == "SUS") # subset different selection operator data
coop_operator_data <- rbind(coop3_data, operator_data)  # combine selection operator data
coop_operator_data$measure <- factor(coop_operator_data$measure)	# remove unused levels

# GTFT data
coop_gtft_data <- subset(sim_data, (measure == "Cooperation" | measure == "GTFTp1") & (skew == "High skew" | skew == "No skew") & psi %in% c(0, 0.25, 0.22))

# Change in action data
abs_change_data <- subset(sim_data, measure == "AbsChange")
net_change_data <- subset(sim_data, measure == "NetChange")

########
### Plot cooperation plots
########
skew_levels <- c("High skew", "Low skew", "No skew")  # assign name of skew conditions

# Plot cooperation rates by forgetting rate, type of memory error, and contact pattern skew
coop_data$skew <- factor(coop_data$skew, levels = c("High skew", "Low skew", "No skew"))
coop_data$memory_error_type <- factor(coop_data$memory_error_type, levels = c("Skipping", "Guessing", "Perfect"))

robust_coop_plot <- ggplot(coop_data, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  geom_point(size = 8) +  # plot points
  geom_rect(aes(xmin = 3.75, xmax = 4.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean cooperation percentage") +  # label axes
  ylim(-10, 101) +  # limit y axis
  scale_color_manual(values = col.blind[1:3], name = "") +  # define colors
  scale_shape_manual(values = c(16, 17, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2, 0), name = "", labels = c("Skipping", "Guessing", "")) +  # define linetypes
  geom_segment(aes(x = 1, y = -10, xend = 7, yend = -10), color = "grey30", arrow = arrow(length = unit(1, "cm")), inherit.aes = FALSE) +  # add arrow
  annotate("text", x = 4, y = -7, label = "More forgetting", size = 12) + # add text "More forgetting"
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.15, 0.3), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/cooperation_rate.png", width = 1000, height = 900)  # open plot device
plot(robust_coop_plot)	# print plot
dev.off()		# close plot device

# Plot cooperation rates by forgetting rate, type of memory error, and contact pattern skew for robustness simulations
coop_robust_means$skew <- factor(coop_robust_means$skew, levels = c("High skew", "Low skew", "No skew"))  # reorder levels
robust_coop_plot <- ggplot(coop_robust_means, aes(x = interactions, y = percent, color = skew, shape = skew)) +  # create plot
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = skew), size = 1.5) +  # plot lines
  labs(x = "Number of interactions", y = "Mean cooperation percentage") +  # label axes
  lims(x = c(0, 1600), y = c(0, 100)) +  # limit y axis
  scale_color_manual(values = col.blind[1:3], name = "") +  # define colors
  scale_shape_manual(values = c(16, 17, 15), name = "") +  # define shapes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.15, 0.15), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/robust_cooperation.png", width = 1000, height = 900)  # open plot device
plot(robust_coop_plot)	# print plot
dev.off()		# close plot device

# Plot cooperation rates by payoff matrix, forgetting rate, type of memory error, and contact pattern skew
coop_matrix2_data$skew <- factor(coop_matrix2_data$skew, levels = c("High skew", "Low skew", "No skew"))   # reorder skew factors
coop_matrix2_data$memory_error_type <- factor(coop_matrix2_data$memory_error_type, levels = c("Skipping", "Guessing", "Perfect"))   # reorder memory error types

coop_matrix_plot <- ggplot(coop_matrix2_data, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  facet_wrap(~matrix) +  # create separate facets for matrices
  geom_point(size = 8) +  # plot points
  geom_rect(aes(xmin = 1.75, xmax = 2.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean cooperation percentage") +  # label axes
  ylim(0, 100) +  # limit y axis
  scale_color_manual(values = col.blind[c(1, 3)], name = "") +  # define colors
  scale_shape_manual(values = c(16, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2), name = "") +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.9, 0.5), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, strip.text.x = element_text(size = 50, margin = margin(3,0,4,0, "mm")), strip.background = element_rect(fill="grey80"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/cooperation_matrix.png", width = 1500, height = 900)  # open plot device
plot(coop_matrix_plot)	# print plot
dev.off()		# close plot device

# Plot cooperation rates by selection operator, forgetting rate, type of memory error, and contact pattern skew
coop_operator_data$skew <- factor(coop_operator_data$skew, levels = c("High skew", "No skew"))   # reorder skew factors
coop_operator_data$memory_error_type <- factor(coop_operator_data$memory_error_type, levels = c("Skipping", "Guessing"))   # reorder memory error types

coop_operator_plot <- ggplot(coop_operator_data, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  facet_wrap(~measure) +  # create separate facets for operators
  geom_point(size = 8) +  # plot points
  geom_rect(aes(xmin = 0.75, xmax = 1.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean cooperation percentage") +  # label axes
  ylim(0, 100) +  # limit y axis
  scale_color_manual(values = col.blind[c(1, 3)], name = "") +  # define colors
  scale_shape_manual(values = c(16, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2), name = "") +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.22, 0.3), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, strip.text.x = element_text(size = 50), strip.background = element_rect(fill="grey80"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/cooperation_operator.png", width = 1500, height = 900)  # open plot device
plot(coop_operator_plot)	# print plot
dev.off()		# close plot device

# Plot cooperation rates by GTFT type, forgetting rate, type of memory error, and contact pattern skew
coop_gtft_data$skew <- factor(coop_gtft_data$skew, levels = c("High skew", "No skew"))   # reorder skew factors
coop_gtft_data$memory_error_type <- factor(coop_gtft_data$memory_error_type, levels = c("Skipping", "Guessing", "Perfect"))   # reorder memory error types
coop_gtft_data$measure <- factor(coop_gtft_data$measure, labels = c("GTFT: p = 0.99", "GTFT: p = 1.0"))

coop_gtft_plot <- ggplot(coop_gtft_data, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  facet_wrap(~measure) +  # create separate facets for GTFT types
  geom_point(size = 8) +  # plot points
  geom_rect(aes(xmin = 2.75, xmax = 3.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean cooperation percentage") +  # label axes
  ylim(0, 100) +  # limit y axis
  scale_color_manual(values = col.blind[c(1, 3)], name = "") +  # define colors
  scale_shape_manual(values = c(16, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2, 0), name = "", labels = c("Skipping", "Guessing", "")) +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.1, 0.6), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, strip.text.x = element_text(size = 50, margin = margin(3,0,5,0, "mm")), strip.background = element_rect(fill="grey80"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/cooperation_gtft.png", width = 1500, height = 900)  # open plot device
plot(coop_gtft_plot)	# print plot
dev.off()		# close plot device

################
### Create strategy frequency plots
################

########
### Three strategies
########
# Prepare data
three_strategies <- subset(sim_data, measure %in% c("ALLD", "GRIM", "TFT family"))  # subset data from ALLD, GRIM, and TFT family
three_strategies$measure <- factor(three_strategies$measure, levels = c("ALLD", "TFT family", "GRIM"))  # reorder strategy levels
three_strategies$forget_skew <- factor(three_strategies$forget_skew)  # remove unused factor levels
three_strategies$skew <- factor(three_strategies$skew, levels = c("High skew", "Low skew", "No skew"))  # reorder strategy levels
three_strategies$memory_error_type <- factor(three_strategies$memory_error_type, levels = c("Skipping", "Guessing", "Perfect"))  # reorder strategy levels

# Plot strategy frequencies by forgetting rate, type of memory error, and contact pattern skew
three_strategies_plot <- ggplot(three_strategies, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  facet_wrap(~measure, as.table = F) +  # create separate facets for strategies
  geom_point(size = 5) +  # plot points
  geom_rect(aes(xmin = 3.75, xmax = 4.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean strategy percentage") +  # label axes
  # ylim(0, 100) +  # limit y axis
  scale_color_manual(values = col.blind[1:3], name = "") +  # define colors
  scale_shape_manual(values = c(16, 17, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2, 0), name = "", labels = c("Skipping", "Guessing", "")) +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.9, 0.75), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, strip.text.x = element_text(size = 50, margin = margin(3,0,4,0, "mm")), strip.background = element_rect(fill="grey80"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/three_strategies.png", width = 2400, height = 900)  # open plot device
plot(three_strategies_plot)	# print plot
dev.off()		# close plot device

########
### All strategies
########
# Prepare data
all_strategies <- subset(sim_data, measure %in% c("ALLC", "ALLD", "CTFT", "GRIM", "GTFT", "RAND", "TF2T", "TFT", "WSLS"))  # subset data from nine strategies
all_strategies$forget_skew <- factor(all_strategies$forget_skew)  # remove unused factor levels
all_strategies$skew <- factor(all_strategies$skew, levels = c("High skew", "Low skew", "No skew"))  # reorder strategy levels
all_strategies$memory_error_type <- factor(all_strategies$memory_error_type, levels = c("Skipping", "Guessing", "Perfect"))  # reorder strategy levels

# Plot strategy frequencies by forgetting rate, type of memory error, and contact pattern skew
all_strategies_plot <- ggplot(all_strategies, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  facet_wrap(~measure) +  # create separate facets for strategies
  geom_point(size = 5) +  # plot points
  geom_rect(aes(xmin = 3.75, xmax = 4.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean strategy percentage") +  # label axes
  # ylim(0, 100) +  # limit y axis
  scale_color_manual(values = col.blind[1:3], name = "") +  # define colors
  scale_shape_manual(values = c(16, 17, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2, 0), name = "", labels = c("Skipping", "Guessing", "")) +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=90), axis.text=element_text(size=31), legend.text=element_text(size=50), legend.position = c(0.075, 0.9), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, strip.text.x = element_text(size = 70, margin = margin(3,0,4,0, "mm")), strip.background = element_rect(fill="grey80"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 1.5, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/all_strategies.png", width = 2400, height = 2700)  # open plot device
plot(all_strategies_plot)	# print plot
dev.off()		# close plot device

########
### Payoff matrix comparison
########
# Prepare data
payoff_strategies <- subset(sim_data, measure %in% c("WSLS", "TFT family", "WSLS2", "TFT family2"))  # subset data from ALLD, GRIM, and TFT family
payoff_strategies$measure <- factor(payoff_strategies$measure, levels = c("WSLS", "TFT family", "WSLS2", "TFT family2"))  # reorder strategy levels
payoff_strategies$forget_skew <- factor(payoff_strategies$forget_skew)  # remove unused factor levels
payoff_mem_cond <- unique(subset(payoff_strategies, measure %in% c("WSLS2", "TFT family2"))$memory_condition)  # find memory conditions for matrix 2 data
payoff_forget_skew <- unique(subset(payoff_strategies, measure %in% c("WSLS2", "TFT family2"))$forget_skew)  # find forgetting type and skew conditions for matrix 2 data
payoff_strategies_matrix2 <- subset(payoff_strategies, memory_condition %in% payoff_mem_cond & forget_skew %in% payoff_forget_skew)  # remove conditions not tested with matrix 2
payoff_strategies_matrix2$matrix <- as.factor(ifelse(payoff_strategies_matrix2$measure %in% c("WSLS2", "TFT family2"), "Matrix 2", "Matrix 1"))  # subset relevant data

payoff_strategies_matrix2$skew <- factor(payoff_strategies_matrix2$skew, levels = c("High skew", "No skew"))   # reorder skew factors
payoff_strategies_matrix2 <- subset(payoff_strategies_matrix2, memory_error_type != "Perfect")  # remove perfect memory conditions
payoff_strategies_matrix2$memory_error_type <- factor(payoff_strategies_matrix2$memory_error_type, levels = c("Skipping", "Guessing"))   # reorder memory error types
payoff_strategies_matrix2$measure <- factor(payoff_strategies_matrix2$measure, levels = c("WSLS", "WSLS2", "TFT family", "TFT family2"))  # reorder strategy order
payoff_strategies_matrix2$measure <- factor(payoff_strategies_matrix2$measure, labels = c("WSLS\nOriginal matrix", "WSLS\nNew matrix", "TFT family\nOriginal matrix", "TFT family\nNew matrix")) # relabel measures

# Plot strategy frequencies by forgetting rate, type of memory error, and contact pattern skew
payoff_strategies_plot <- ggplot(payoff_strategies_matrix2, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  facet_wrap(~measure, as.table = F) +  # create separate facets for strategies
  geom_point(size = 5) +  # plot points
  geom_rect(aes(xmin = 1.75, xmax = 2.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean strategy percentage") +  # label axes
  # ylim(0, 100) +  # limit y axis
  scale_color_manual(values = col.blind[c(1, 3)], name = "") +  # define colors
  scale_shape_manual(values = c(16, 17), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2), name = "") +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.88, 0.3), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, strip.text.x = element_text(size = 50, margin = margin(3,0,4,0, "mm")), strip.background = element_rect(fill="grey80"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/matrix_strategies.png", width = 1200, height = 1200)  # open plot device
plot(payoff_strategies_plot)	# print plot
dev.off()		# close plot device

################
### Create change in action plot
################
## Absolute change
abs_change_data$skew <- factor(abs_change_data$skew, levels = c("High skew", "No skew"))  # reorder levels
abs_change_data$memory_error_type <- factor(abs_change_data$memory_error_type, levels = c("Skipping", "Guessing"))  # reorder levels

abs_change_plot <- ggplot(abs_change_data, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  geom_point(size = 8) +  # plot points
  geom_rect(aes(xmin = 1.75, xmax = 2.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  labs(x = "Forgetting parameters", y = "Mean total percent\nchange in actions") +  # label axes
  ylim(0, 25) +  # limit y axis
  scale_color_manual(values = col.blind[c(1, 3)], name = "") +  # define colors
  scale_shape_manual(values = c(16, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2), name = "", labels = c("Skipping", "Guessing", "")) +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.8, 0.5), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/abs_change.png", width = 1000, height = 900)  # open plot device
plot(abs_change_plot)	# print plot
dev.off()		# close plot device

## Net change
net_change_data$skew <- factor(net_change_data$skew, levels = c("High skew", "No skew"))  # reorder levels
net_change_data$memory_error_type <- factor(net_change_data$memory_error_type, levels = c("Skipping", "Guessing"))  # reorder levels

net_change_plot <- ggplot(net_change_data, aes(x = forget_rate, y = percent, color = skew, shape = skew, linetype = memory_error_type)) +  # create plot
  geom_point(size = 8) +  # plot points
  geom_rect(aes(xmin = 1.75, xmax = 2.25), ymin = -Inf, ymax = Inf, fill = "grey90", inherit.aes = FALSE) +  # add shading for empirical parameters
  geom_point(size = 8) +  # plot points
  geom_line(aes(group = forget_skew), size = 1.5) +  # plot lines
  geom_hline(yintercept = 0, col = "grey50") + # add line at 0
  labs(x = "Forgetting parameters", y = "Mean percent shift\ntoward cooperation") +  # label axes
  scale_color_manual(values = col.blind[c(1, 3)], name = "") +  # define colors
  scale_shape_manual(values = c(16, 15), name = "") +  # define shapes
  scale_linetype_manual(values = c(1, 2), name = "", labels = c("Skipping", "Guessing", "")) +  # define linetypes
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = c(0.8, 0.5), legend.key.size = unit(2.75, 'lines'), legend.title = NULL, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/net_change.png", width = 1000, height = 900)  # open plot device
plot(net_change_plot)	# print plot
dev.off()		# close plot device

################
### Create memory functions plot
################
# Prepare data
memory_data <- data.frame(lambda = rep(c(1, 1, 1, 0.87, 0.6, 0.4, 0.25), each = 50), psi = rep(c(0, 0.1, 0.25, 0.25, 0.22, 0.25, 0.5), each = 50), interactions = rep(0:49, times = 7))  # create memory data frame
memory_data$accuracy <- memory_data$lambda * (memory_data$interactions + 1) ^ - (memory_data$psi)   # calculate memory accuracy
memory_data$forget_rate <- as.factor(paste(paste("λ=", round(memory_data$lambda, 2), sep = ""), paste("Ψ=", round(memory_data$psi, 2), sep = ""), sep = "\n"))      # create forget_rate column
memory_data$forget_rate <- factor(memory_data$forget_rate)	# remove unused levels

# Plot memory functions for different forgetting parameter values
memory_function_plot <- ggplot(memory_data, aes(x = interactions, y = accuracy, color = forget_rate)) + # create plot
  geom_line(aes(group = forget_rate), size = 1.5) + # plot lines
  labs(x = "Intervening interactions (k)", y = "Probability of remembering") +  # label axes
  ylim(0, 1) +  # set y-axis limits
  scale_color_manual(values = col.blind, name = "") +  # define colors
  annotate("text", x = 40, y = 0.96, label = "λ=1, Ψ=0", size = 12) + # add line label
  annotate("text", x = 40, y = 0.74, label = "λ=1, Ψ=0.1", size = 12) + # add line label
  annotate("text", x = 40, y = 0.44, label = "λ=1, Ψ=0.25", size = 12) + # add line label
  annotate("text", x = 40, y = 0.31, label = "λ=0.87, Ψ=0.22", size = 12) + # add line label
  annotate("text", x = 40, y = 0.22, label = "λ=0.60, Ψ=0.25", size = 12) + # add line label
  annotate("text", x = 40, y = 0.12, label = "λ=0.40, Ψ=0.25", size = 12) + # add line label
  annotate("text", x = 40, y = 0.01, label = "λ=0.25, Ψ=0.5", size = 12) + # add line label
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = "none", legend.title = NULL, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/memory_functions.png", width = 1000, height = 900)  # open plot device
plot(memory_function_plot)	# print plot
dev.off()			# close plot device

################
### Create skew plot
################
# Prepare data
no_skew <- rep(10, 10)														# create vector of no skew contact frequencies
low_skew <- c(19, 16, 14, 11, 10, 8, 7, 6, 5, 4)	# create vector of low skew contact frequencies
high_skew <- c(33, 23, 15, 10, 7, 5, 3, 2, 1, 1)	# create vector of high skew contact frequencies
members <- 1:10																		# create sequence of network members
p <- seq(0.01, 1, 0.1)														# create sequence of probabilities
empirical <- 100 * (log(p / 1.09)) / -13.5				# create empirical skew pattern
skews <- c("No skew", "Low skew", "High skew", "Empirical skew")		# create vector of skew labels
skew_functions <- data.frame(skew = rep(skews, each = 10), member = rep(members, times = 4), contacts = c(no_skew, low_skew, high_skew, empirical))	# create data frame of skew data
skew_functions$skew <- factor(skew_functions$skew, levels = skews)	# remove unused skew levels

# Plot number of interactions for each network member for each skew type
skew_plot <- ggplot(skew_functions, aes(x = member, y = contacts, color = skew, linetype = skew)) + # create plot
  geom_line(aes(group = skew), size = 3) + # plot lines
  labs(x = "Network member", y = "Number of interactions\nwith network member") +  # label axes
  scale_color_manual(values = c(col.blind[3:1], "black"), name = "") +  # define colors
  scale_linetype_manual(values = c(1, 1, 1, 2), name = "") +  # define colors
  scale_x_continuous(breaks = 1:10) +   # control tick mark placement
  annotate("text", x = 3.25, y = 25, label = "High skew", size = 15, col = col.blind[1]) + # add line label
  annotate("text", x = 9, y = 7.5, label = "Low skew", size = 15, col = col.blind[2]) + # add line label
  annotate("text", x = 8, y = 11.5, label = "No skew", size = 15, col = col.blind[3]) + # add line label
  annotate("text", x = 4.7, y = 20, label = "Empirical skew", size = 15, col = "black") + # add line label
  geom_segment(aes(x = 2.05, xend = 3, y = 18, yend = 20), col = "black") + # add line pointing to empirical line
  theme_bw() +  # use bw theme
  theme(axis.title=element_text(size=55), axis.text=element_text(size=30), legend.text=element_text(size=35), legend.position = "none", legend.title = NULL, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.75, 0.5), "lines"))  # adjust font sizes, legends
png(filename = "figures/skew.png", width = 1000, height = 900)  # open plot device
plot(skew_plot)	# print plot
dev.off()		# close plot device
