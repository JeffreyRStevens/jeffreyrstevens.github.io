###
### stevens_duque_2018_rcode.R
### Created by Jeffrey R. Stevens on 17 Jun 2016 (jeffrey.r.stevens@gmail.com),
###	 finalized on 12 Jul 2018
### Summary: This script analyzes and produces plots for citation bias analysis.
### Instructions: Place this file and the data files in the main directory.
### 	Create a folder called "figures". Set the R working directory to the 
### 	main directory.  At the R command prompt, type
### 	> source("stevens_duque_2018_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving figures in the figures directory. If packages do not
###		load properly, install them with install.packages("package_name").
### Data files:
###  stevens_duque_2018_data.csv
###   data_set - signals first or second data set
### 	journal - name of the journal that published the article
###   field - field for journal
###   citation_style - citation style for journal at that time
###   year - publication year for the article
### 	author - article's first author's surname and initials
### 	times_cited - number of times the article was cited according to Web of Science
###			as of Jun 2016 (data set 1) or Apr 2017 (data set 2)
### License: This script is released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
###   International license (CC BY-NC-SA 4.0). You may share and adapt this content with attribution, for 
###   non-commercial purposes if you ShareAlike (distribute any contributions under the same license).
###


# Load libraries ----------------------------------------------------------

library(BayesFactor)  # needed for Bayesian analysis
library(car)			    # needed for Recode
library(foreach)	    # needed for foreach
library(grid)			    # needed for grid.text
library(lme4)         # needed for linear mixed models
library(papaja)			  # needed for APA formatting
library(scales)       # needed for comma
library(tidyverse)    # needed for tidyverse functions

# Define functions ----------------------------------------------------------

###
## Create function that is the inverse of %in%, that is it finds items that do not match the vector
###
"%notin%" <- function(x, table) {  # create function analogous to %in% that searches for items not in a vector
  match(x, table, nomatch = 0) == 0
}

###
## Convert BIC values to Bayes factor (from Wagenmakers, E.-J. (2007). A practical solution to the pervasive problems of p values. Psychonomic Bulletin & Review, 14(5), 779–804. https://dx.doi.org/10.3758/BF03194105)
###
bic_bf10 <- function(null, alternative) {
  new_bf <- exp((null - alternative) / 2) # convert BICs to Bayes factor
  names(new_bf) <- NULL   # remove BIC label
  return(new_bf)  # return Bayes factor of alternative over null hypothesis
}

###
## Extract the top quantile of data for each letter
###
extractTopQuantile <- function(citation_df) {
  top_quantile <- citation_df[1, ]  # initiate data frame with first row of data
  foreach(letter = LETTERS) %do% {  # for each letter
    current_letter_data <- subset(citation_df, author == letter)  # subset the data for this letter
    current_quantile <- quantile(current_letter_data$times_cited, prob = seq(0, 0.9, length = num_quantiles))[quant]  # find cutoff for specified quantile
    current_top_data <- subset(current_letter_data, times_cited >= current_quantile)  # subset the data equal to or above the quantile cutoff
    top_quantile <- rbind(top_quantile, current_top_data) # append data
  }
  top_quantile <- top_quantile[-1, ] # remove initial row
}

###
## Create function that calculates correlation coefficent, slope, intercept, and BF for each data group
###
calculateCorrelationPerGroup <- function(citation_rate_group, group_names) {
  # Prepare new data frames
  group_citation_rate_df <- data.frame(matrix(rep(NA, 6), ncol = 6))		# initiate data frame for author data
  names(group_citation_rate_df) <- c("group", "author", "author_num", "mean_prob_cited", "mean_citation_rate", "percent_citation_rate")	# label data frame columns
  names(citation_rate_group) <- names(group_citation_rate_df)[-length(group_citation_rate_df[1, ])]   # remove last label
  group_citation_rate_stats <- group_prob_cited_stats <- data.frame(matrix(rep(NA, 7 * length(group_names)), ncol = 7)) # initiate data frame for statistics
  names(group_citation_rate_stats) <- names(group_prob_cited_stats) <- c("group", "tau", "tau_lowerCI", "tau_upperCI", "beta", "intercept", "bf")  # label data frame columns
  group_counter <- 1 	# initiate group counter
  
  # Calculate stats for each group
  foreach(current_group = group_names) %do% {	# for each group
    # Filter data for this group and calculate percent citation rate
    if(length(group_names) > 1) { # if there is more than one group
      citation_letters_group <- subset(citation_rate_group, group == current_group)	# filter data based on group
    } else {  # if it is the overall analysis
      citation_letters_group <- citation_rate_group  # use all data
    }
    citation_letters_group$percent_citation_rate <- citation_letters_group$mean_citation_rate / sum(citation_letters_group$mean_citation_rate, na.rm = TRUE) * 100 	# calculate percentage for each letter
    citation_letters_group <- data.frame(citation_letters_group)	# convert to data frame
    group_citation_rate_df <- rbind(group_citation_rate_df, citation_letters_group)	# append to data frame of all groups
    # Create data frame of citation rate statistics per group
    group_citation_rate_stats$group[group_counter] <- group_prob_cited_stats$group[group_counter] <- current_group   # assign group
    group_citation_rate_stats$tau[group_counter] <- cor.test(citation_letters_group$percent_citation_rate, citation_letters_group$author_num, method = "kendall")$estimate		# extract Kendall correlation coefficient
    group_citation_rate_stats$tau_lowerCI[group_counter] <- as.numeric(credibleIntervalKendallTau(group_citation_rate_stats$tau[group_counter], 26)[1])  # calculate lower Bayesian credible interval
    group_citation_rate_stats$tau_upperCI[group_counter] <- as.numeric(credibleIntervalKendallTau(group_citation_rate_stats$tau[group_counter], 26)[3])  # calculate upper Bayesian credible interval
    group_citation_rate_stats$bf[group_counter] <- bfCorrieKernelKendallTau(group_citation_rate_stats$tau[group_counter], 26)$bf10  # calculate Bayes factor for Kendall's tau of percent times cited ~ author
    group_citation_rate_stats$beta[group_counter] <- lm(percent_citation_rate ~ author_num, citation_letters_group)$coefficients[2] # calculate linear regression slope
    group_citation_rate_stats$intercept[group_counter] <- lm(percent_citation_rate ~ author_num, citation_letters_group)$coefficients[1]  # calculate linear regression intercept
    # Create data frame of probability of citation statistics per group
    group_prob_cited_stats$tau[group_counter] <- cor.test(citation_letters_group$mean_prob_cited, citation_letters_group$author_num, method = "kendall")$estimate		# extract Kendall correlation coefficient
    group_prob_cited_stats$tau_lowerCI[group_counter] <- as.numeric(credibleIntervalKendallTau(group_prob_cited_stats$tau[group_counter], 26)[1])  # calculate lower Bayesian credible interval
    group_prob_cited_stats$tau_upperCI[group_counter] <- as.numeric(credibleIntervalKendallTau(group_prob_cited_stats$tau[group_counter], 26)[3])  # calculate upper Bayesian credible interval
    group_prob_cited_stats$beta[group_counter] <- lm(mean_prob_cited ~ author_num, citation_letters_group)$coefficients[2] # calculate linear regression slope
    group_prob_cited_stats$intercept[group_counter] <- lm(mean_prob_cited ~ author_num, citation_letters_group)$coefficients[1]  # calculate linear regression intercept
    group_prob_cited_stats$bf[group_counter] <- bfCorrieKernelKendallTau(group_prob_cited_stats$tau[group_counter], 26)$bf10  # calculate Bayes factor for Bayesian linear regression of percent times cited ~ author
    group_counter <- group_counter + 1	# increment counter
  }
  group_citation_rate_df <<- group_citation_rate_df[-1, ]  # remove initiation row
  group_citation_rate_stats <<- group_citation_rate_stats   # save citation rate statistics
  group_prob_cited_stats <<- group_prob_cited_stats   # save probability of citation statistics
}

###
## Calculate Bayes factor for Kendall's Tau (from van Doorn, J., Ly, A., Marsman, M., & Wagenmakers, E.-J. (in press). Bayesian inference for Kendall’s rank correlation coefficient. The American Statistician. https://doi.org/10.1080/00031305.2016.1264998). Script from https://osf.io/b9qhj/.
###
# Prior specification Kendall's Tau
scaledBetaTau <- function(tau, alpha=1, beta=1){
  result <-   ((pi*2^(-2*alpha))/beta(alpha,alpha))  * cos((pi*tau)/2)^(2*alpha-1)
  return(result)
}

priorTau <- function(tau, kappa){
  scaledBetaTau(tau, alpha = (1/kappa), beta = (1/kappa))
}

priorTauPlus <- function(tau, kappa=1) {
  non.negative.index <- tau >=0
  less.than.one.index <- tau <=1
  value.index <- as.logical(non.negative.index*less.than.one.index)
  result <- tau*0
  result[value.index] <- 2*priorTau(tau[value.index], kappa)
  return(result)
}

priorTauMin <- function(tau, kappa=1) {
  negative.index <- tau <=0
  greater.than.min.one.index <- tau >= -1
  value.index <- as.logical(negative.index*greater.than.min.one.index)
  result <- tau*0
  result[value.index] <- 2*priorTau(tau[value.index], kappa)
  return(result)
}

# Posterior specification Kendall's Tau
postDensKendallTau <- function(delta,Tstar,n,kappa=1,var=var,test="two-sided"){ 
  if(test == "two-sided"){priorDens <- priorTau(delta,kappa)
  } else if(test == "positive"){priorDens <- priorTauPlus(delta,kappa)
  } else if(test == "negative"){priorDens <- priorTauMin(delta,kappa)}
  priorDens <- priorTau(delta,kappa)
  dens <- dnorm(Tstar,(1.5*delta*sqrt(n)),sd=sqrt(var))* priorDens
  return(dens)
}
posteriorTau <- function(delta,kentau,n,kappa=1,var=1,test="two-sided"){
  Tstar <- (kentau * ((n*(n-1))/2))/sqrt(n*(n-1)*(2*n+5)/18)
  var <- min(1,var)
  if(test == "two-sided"){lims <- c(-1,1)
  } else if(test == "positive"){lims <- c(0,1)
  } else if(test == "negative"){lims <- c(-1,0)}
  logicalCensor <- (delta >= lims[1] & delta <= lims[2])
  dens <- logicalCensor*postDensKendallTau(delta,Tstar,n,kappa,var,test=test)/
    integrate(function(delta){postDensKendallTau(delta,Tstar,n,kappa,var,test=test)},lims[1],lims[2])$value
} 

# Bayes factor computation Kendall's Tau
bfCorrieKernelKendallTau <- function(tau, n, kappa=1, var=1, ciValue=0.95){ 
  tempList <- list(vector())
  output <- list(n=n, r=tau, bf10=NA, bfPlus0=NA, bfMin0=NA)
  output$bf10 <- priorTau(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="two-sided")
  output$bfPlus0 <- priorTauPlus(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="positive")
  output$bfMin0 <- priorTauMin(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="negative")
  return(output)
}

# Compute credible intervals kendalls tau
credibleIntervalKendallTau <- function(kentau,n,kappa=1,var=1, test="two-sided", ciValue = 0.95){
  nSeqs <- 1000
  lowCI <- (1-ciValue)/2
  upCI <- (1+ciValue)/2
  taus <- seq(-1,1,length.out = (nSeqs-1))
  densVals <- posteriorTau(taus, kentau, n, kappa = kappa, var = var, test = test)
  densVals <- cumsum((densVals[1:(nSeqs-1)]+densVals[2:nSeqs])*0.5*(taus[2]-taus[1]))
  lowerCI <- taus[which(densVals>=lowCI)[1]]
  upperCI <- taus[which(densVals>=upCI)[1]]
  median <- taus[which(densVals>=0.5)[1]]
  return(list(lowerCI = lowerCI, median = median, upperCI = upperCI))
}

sampleTausA <- function(myTau,myN,nSamples = 3e3, var = 1){
  nSeqs <- 1000
  tauSamples <- NULL
  taus <- seq(-1,1,length.out = nSeqs)
  densVals <- posteriorTau(taus, myTau, myN, var = var)
  ceiling <- max(densVals)
  lowerB <- taus[which(round(densVals,digits=6) != 0 )][1]
  upperB <- rev(taus[which(round(densVals,digits=6) != 0 )])[1]
  
  while(length(tauSamples) < nSamples){
    prop <- runif(1,lowerB,upperB)
    propDens <- posteriorTau(prop, myTau, myN, var = var)
    if(propDens > runif(1,0,ceiling)){tauSamples <- c(tauSamples,prop)}
  }
  return(tauSamples)
}

# Data set 1 ----------------------------------------------------------------

###
## Input and process data file
###
citation_data <- read_csv2("stevens_duque_2018_data.csv", progress = FALSE)     # input data
citation_data <- citation_data %>% filter(!is.na(citation_style))               # filter data_set 1 data
citation_data1 <- citation_data %>% filter(data_set == 1)               # filter data_set 1 data
citation_data1 <- citation_data1 %>% mutate(author_num = as.numeric(as.factor(author)) - 1, year_num = year - 2000)  # convert letters to numbers 0-25 and years to numbers 0-15
citation_data1$prob_cited <- ifelse(citation_data1$times_cited > 0, 1, 0)   # convert times cited to probability of citation

journal_names1 <- unique(citation_data1$journal)  					# create vector of journals
citation_style_names1 <- unique(citation_data1$citation_style) 	# create vector of citation styles
field_names1 <- unique(citation_data1$field)     				# create vector of fields

###
## Aggregate data for sample size calculations
###
cited_df1 <- citation_data1 %>% group_by(journal, field, citation_style) %>% summarize(mean_citation_rate = round(mean(times_cited), 2), num_articles = length(times_cited))			# summarize citation count by journal, field, and citation style
times_cited_overall_n1 <- format(sum(cited_df1$num_articles), big.mark = ",", trim=TRUE)	# calculate overall number of articles
times_cited_field_n1 <- format((cited_df1 %>% group_by(field) %>% summarize(num_articles = sum(num_articles)))$num_articles, big.mark = ",", trim=TRUE)	# calculate number of articles for each field
times_cited_citation_style_n1 <- format((cited_df1 %>% group_by(citation_style) %>% summarize(num_articles = sum(num_articles)))$num_articles, big.mark = ",", trim=TRUE)[1:3]	# calculate number of articles for each citation style
times_cited_journal1 <- cited_df1 %>% group_by(journal, field) %>% summarize(num_articles = sum(num_articles))
times_cited_journal1$n <- format(times_cited_journal1$num_articles, big.mark = ",", trim=TRUE)	# calculate number of articles per journal

###
## Calculate effect of letter on citation rate per citation style
###
# Calculate citation rate per citation style
citation_rate_citation_style1 <- citation_data1 %>% group_by(citation_style, author, author_num) %>% summarize(mean_prob_cited = mean(prob_cited), mean_citation_rate = mean(times_cited))	# sum times cited by letter

# Calculate correlations, slopes, intercepts, and Bayes factors for citation style
calculateCorrelationPerGroup(citation_rate_citation_style1, citation_style_names1)
citation_style_df1 <- group_citation_rate_df    # save data
citation_style_citation_rate_stats1 <- group_citation_rate_stats  # save stats
citation_style_prob_cited_stats1 <- group_prob_cited_stats  # save stats

# Calculate Bayesian linear regression comparing citation style * letter interaction to without interaction
citation_style_null_bf1 <- lmBF(percent_citation_rate ~ author_num + group, data =  citation_style_df1, progress = FALSE) # calculate Bayes factor for null model with just main effects
citation_style_bf1 <- lmBF(percent_citation_rate ~ author_num * group, data =  citation_style_df1, progress = FALSE) # calculate Bayes factor for null model with main effects and interaction
citation_style_bf1 / citation_style_null_bf1  # calculate Bayes factor comparing interaction to null model

# Plot mean citation rate per letter and citation style
citation_style_tau_label1 <- paste("tau ==", round(citation_style_citation_rate_stats1$tau, 2))
times_cited_citation_style_ggplot1 <- ggplot(citation_style_df1, aes(x = author, y = percent_citation_rate)) +
  facet_grid(~group) +  # for each group
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), citation_style_citation_rate_stats1, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  annotate("text", 13, 5.8, label = citation_style_tau_label1, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 5.5, label = paste("BF =", round(citation_style_citation_rate_stats1$bf, 2)), size = 18) +  # add Bayes factor
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=40), strip.text.x = element_text(size = 60, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/citation_styles1.png"), width = 2000, height = 800)	# open PNG device
plot(times_cited_citation_style_ggplot1)		# plot figure
dev.off()		# close device

###
## Calculate effect of letter on citation rate per field
###
# Calculate citation rate per field
citation_rate_field1 <- citation_data1 %>% group_by(field, author, author_num) %>% summarize(mean_prob_cited = mean(prob_cited), mean_citation_rate = mean(times_cited))	# sum times cited by letter

# Calculate correlations, slopes, intercepts, and Bayes factors for citation style
calculateCorrelationPerGroup(citation_rate_field1, field_names1)
field_df1 <- group_citation_rate_df    # save data
field_citation_rate_stats1 <- group_citation_rate_stats  # save stats
field_prob_cited_stats1 <- group_prob_cited_stats  # save stats
field_df1$group <- factor(field_df1$group, levels = c("Psychology", "Biology")) # reorder fields
field_citation_rate_stats1$group <- factor(field_citation_rate_stats1$group, levels = c("Psychology", "Biology")) # reorder fields
field_prob_cited_stats1$group <- factor(field_prob_cited_stats1$group, levels = c("Psychology", "Biology")) # reorder fields

# Calculate Bayesian linear regression comparing field * letter interaction to without interaction
field_null_bf1 <- lmBF(percent_citation_rate ~ author_num + group, data =  field_df1, progress = FALSE) # calculate Bayes factor for null model with just main effects
field_bf1 <- lmBF(percent_citation_rate ~ author_num * group, data =  field_df1, progress = FALSE) # calculate Bayes factor for null model with main effects and interaction
field_bf1 / field_null_bf1  # calculate Bayes factor comparing interaction to null model

# Plot mean citation rate per letter and field
field_tau_label1 <- paste("tau ==", round(field_citation_rate_stats1$tau, 2))
times_cited_field_ggplot1 <- ggplot(field_df1, aes(x = author, y = percent_citation_rate)) +
  facet_grid(~group) +  # for each group
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), field_citation_rate_stats1, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  annotate("text", 13, 5.8, label = field_tau_label1, size = 18, parse = TRUE) +  # add tau value  
  annotate("text", 13, 5.5, label = paste("BF =", round(field_citation_rate_stats1$bf, 2)), size = 18) +  # add Bayes factor
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=40), strip.text.x = element_text(size = 60, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/fields1.png"), width = 1400, height = 800)	# open PNG device
plot(times_cited_field_ggplot1)		# plot figure
dev.off()		# close device

###
## Calculate effect of letter on citation rate per journal
###
# Calculate citation rate per journal
citation_rate_journal1 <- citation_data1 %>% group_by(journal, author, author_num) %>% summarize(mean_prob_cited = mean(prob_cited), mean_citation_rate = mean(times_cited))		# sum times cited by letter and journal

# Calculate correlations, slopes, intercepts, and Bayes factors for citation style
calculateCorrelationPerGroup(citation_rate_journal1, journal_names1)
journal_df1 <- group_citation_rate_df    # save data
journal_citation_rate_stats1 <- group_citation_rate_stats  # save stats
journal_df1 <- merge(journal_df1, times_cited_journal1[, c("journal", "field")], by.x = "group", by.y = "journal")  # add field column
journal_citation_rate_stats1 <- merge(journal_citation_rate_stats1, times_cited_journal1[, c("journal", "field")], by.x = "group", by.y = "journal")  # add field column

## Psychology data
# Prepare data
psychology_journal_labels1 <- c("American Psychologist", "Animal Cognition", "Annual Review of Psychology", "Current Directions in\nPsychological Science", "Journal of Comparative Psychology", "Journal of Experimental Psychology\nAnimal Learning and Cognition", "Journal of the\nExperimental Analysis of Behavior", "Learning and Behavior", "Perspectives on\nPsychological Science", "Psychological Bulletin", "Psychological Review", "Psychological Science", "Trends in Cognitive Sciences")		# create vector of psychology journal names

# Filter psychology journals
psychology_journal_df1 <- filter(journal_df1, field == "Psychology")  # filter out psychology journals
psychology_journal_df1$journal <- factor(psychology_journal_df1$group, labels = psychology_journal_labels1)  # relabel journals
psychology_citation_rate_stats1 <- filter(journal_citation_rate_stats1, field == "Psychology")  # filter out psychology journals for statistics
psychology_citation_rate_stats1$journal <- psychology_journal_labels1  # relabel journals
psychology_journal_stats1 <- filter(times_cited_journal1, field == "Psychology")  # filter out psychology journals

# Plot mean citation rate per psychology journal
psychology_tau_label1 <- paste("tau ==", round(psychology_citation_rate_stats1$tau, 2))
psychology_ggplot1 <- ggplot(psychology_journal_df1, aes(x = author, y = percent_citation_rate)) +
  facet_wrap(~journal) +  # for each journal
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), psychology_citation_rate_stats1, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  annotate("text", 13, 11, label = psychology_tau_label1, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 9.5, label = paste("BF =", round(psychology_citation_rate_stats1$bf, 1)), size = 18) +  # add Bayes factor
  annotate("text", 13, 8, label = paste("N =", psychology_journal_stats1$n), size = 18) +  # add tau value
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=30), strip.text.x = element_text(size = 25, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/psychology_journals1.png"), width = 1800, height = 2200)	# open PNG device
plot(psychology_ggplot1)		# plot figure
dev.off()		# close device

## Biology data
# Prepare data
biology_journal_labels1 <- c("Animal Behaviour", "Annual Review of\nEcology, Evolution,and Systematics", "Behavioral Ecology and Sociobiology", "Behavioral Ecology", "Behaviour", "Biological Reviews", "Ethology", "Journal of Ethology", "Philosophical Transactions\nof the Royal Society", "PLOS Biology", "Proceedings of the Royal Society", "Quarterly Review of Biology", "Trends in Ecology and Evolution")		# create vector of biology journal names

# Filter biology journals
biology_journal_df1 <- filter(journal_df1, field == "Biology")  # filter out biology journals
biology_journal_df1$journal <- factor(biology_journal_df1$group, labels = biology_journal_labels1)  # relabel journals
biology_citation_rate_stats1 <- filter(journal_citation_rate_stats1, field == "Biology")  # filter out biology journals for statistics
biology_citation_rate_stats1$journal <- biology_journal_labels1  # relabel journals
biology_journal_stats1 <- filter(times_cited_journal1, field == "Biology")  # filter out psychology journals

# Plot mean citation rate per biology journal
biology_tau_label1 <- paste("tau ==", round(biology_citation_rate_stats1$tau, 2))
biology_ggplot1 <- ggplot(biology_journal_df1, aes(x = author, y = percent_citation_rate)) +
  facet_wrap(~journal) +  # for each journal
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), biology_citation_rate_stats1, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  annotate("text", 13, 15, label = biology_tau_label1, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 13, label = paste("BF =", round(biology_citation_rate_stats1$bf, 1)), size = 18) +  # add Bayes factor
  annotate("text", 13, 11, label = paste("N =", biology_journal_stats1$n), size = 18) +  # add N
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=30), strip.text.x = element_text(size = 25, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/biology_journals1.png"), width = 1800, height = 2200)	# open PNG device
plot(biology_ggplot1)		# plot figure
dev.off()		# close device

# Data set 2 -----------------------------------------------

###
## Input and process data file
###
citation_data2 <- citation_data %>% filter(data_set == 2)               # filter data_set 1 data
citation_data2 <- citation_data2 %>% mutate(author_num = as.numeric(as.factor(author)) - 1, year_num = year - 2000)  # convert letters to numbers 0-25 and years to numbers 0-15
citation_data2$prob_cited <- ifelse(citation_data2$times_cited > 0, 1, 0)   # convert times cited to probability of citation

field_names2 <- unique(citation_data2$field)     				# create vector of fields

###
## Aggregate data for sample size calculations
###
cited_df2 <- citation_data2 %>% group_by(journal, field, citation_style) %>% summarize(mean_citation_rate = round(mean(times_cited), 2), num_articles = length(times_cited))			# summarize citation count by journal, field, and citation style
times_cited_overall_n2 <- format(sum(cited_df2$num_articles), big.mark = ",", trim=TRUE)	# calculate overall number of articles
times_cited_field_n2 <- format((cited_df2 %>% group_by(field) %>% summarize(num_articles = sum(num_articles)))$num_articles, big.mark = ",", trim=TRUE)	# calculate number of articles for each field
times_cited_journal2 <- cited_df2 %>% group_by(journal, field) %>% summarize(num_articles = sum(num_articles))	# calculate number of articles per journal
times_cited_journal2$n <- format(times_cited_journal2$num_articles, big.mark = ",", trim=TRUE)  # format number of articles
journal_names2 <- unique(citation_data2$journal)  					# create vector of journals

###
## Calculate effect of letter on citation rate per field
###
# Calculate citation rate per field
citation_rate_field2 <- citation_data2 %>% group_by(field, author, author_num) %>% summarize(mean_prob_cited = mean(prob_cited), mean_citation_rate = mean(times_cited))	# sum times cited by letter
field_df2 <- field_df_nooutliers2 <- data.frame(matrix(rep(NA, 6), ncol = 6))		# initiate data frame for author data

# Calculate correlations, slopes, intercepts, and Bayes factors for citation style
calculateCorrelationPerGroup(citation_rate_field2, field_names2)
field_df2 <- group_citation_rate_df    # save data
field_citation_rate_stats2 <- group_citation_rate_stats  # save stats
field_prob_cited_stats2 <- group_prob_cited_stats  # save stats
field_df2$group <- factor(field_df2$group, levels = c("Psychology", "Geoscience")) # reorder fields
field_citation_rate_stats2$group <- factor(field_citation_rate_stats2$group, levels = c("Psychology", "Geoscience")) # reorder fields
field_prob_cited_stats2$group <- factor(field_prob_cited_stats2$group, levels = c("Psychology", "Geoscience")) # reorder fields

# Calculate Bayesian linear regression comparing field * letter interaction to without interaction
field_null_bf2 <- lmBF(percent_citation_rate ~ author_num + group, data =  field_df2, progress = FALSE) # calculate Bayes factor for null model with just main effects
field_bf2 <- lmBF(percent_citation_rate ~ author_num * group, data =  field_df2, progress = FALSE) # calculate Bayes factor for null model with main effects and interaction
field_bf2 / field_null_bf2  # calculate Bayes factor comparing interaction to null model

# Plot mean citation rate per letter and field
field_tau_label2 <- paste("tau ==", round(field_citation_rate_stats2$tau, 2))
times_cited_field_ggplot2 <- ggplot(field_df2, aes(x = author, y = percent_citation_rate)) +
  facet_grid(~group) +  # for each group
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), field_citation_rate_stats2, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  annotate("text", 13, 6, label = field_tau_label2, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 5.7, label = paste("BF =", round(field_citation_rate_stats2$bf, 2)), size = 18) +  # add Bayes factor
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=40), strip.text.x = element_text(size = 60, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/fields2.png"), width = 1400, height = 800)	# open PNG device
plot(times_cited_field_ggplot2)		# plot figure
dev.off()		# close device

###
## Calculate effect of letter on citation rate per journal
###
# Calculate citation rate per journal
citation_rate_journal2 <- citation_data2 %>% group_by(journal, author, author_num) %>% summarize(mean_prob_cited = mean(prob_cited), mean_citation_rate = mean(times_cited))		# sum times cited by letter and journal

# Calculate correlations, slopes, intercepts, and Bayes factors for citation style
calculateCorrelationPerGroup(citation_rate_journal2, journal_names2)
journal_df2 <- group_citation_rate_df    # save data
journal_citation_rate_stats2 <- group_citation_rate_stats  # save stats
journal_df2 <- merge(journal_df2, times_cited_journal2[, c("journal", "field")], by.x = "group", by.y = "journal")  # add field column
journal_citation_rate_stats2 <- merge(journal_citation_rate_stats2, times_cited_journal2[, c("journal", "field")], by.x = "group", by.y = "journal")  # add field column

## Psychology data
# Prepare data
psychology_journal_labels2 <- c("American Journal \nof Orthopsychiatry", "Asian American Journal \nof Psychology", "Attention Perception \n& Psychophysics", "Behavioral Neuroscience", "Behavior Research Methods", "Canadian Journal of \nExperimental Psychology", "Canadian Psychology-\nPsychologie Canadienne", "Cognitive Affective & \nBehavioral Neuroscience", "Contemporary Psychology-\nAPA Review of Books", "Cultural Diversity & \nEthnic Minority Psychology", "Developmental Psychology", "Dreaming", "Emotion", "Experimental and Clinical \nPsychopharmacology", "Group Dynamics-Theory \nResearch and Practice", "Health Psychology", "History of Psychology", "Journal of \nAbnormal Psychology", "Journal of Applied Psychology", "Journal of Consulting and \nClinical Psychology", "Journal of \nCounseling Psychology", "Journal of Diversity \nin Higher Education", "Journal of \nEducational Psychology", "Journal of Experimental \nPsychology-Applied", "Journal of Experimental \nPsychology-General", "Journal of Experimental \nPsychology-HPP", "Journal of Experimental \nPsychology-LMC", "Journal of Family Psychology", "Journal of Occupational \nHealth Psychology", "Journal of Personality \nand Social Psychology", "Law and Human Behavior", "Memory & Cognition", "Military Psychology", "Neuropsychology", "Personality Disorders-Theory \nResearch and Treatment", "Professional Psychology-\nResearch and Practice", "Psychiatric Rehabilitation \nJournal", "Psychoanalytic Psychology", "Psychological Assessment", "Psychological Methods", "Psychological Services", "Psychological Trauma-Theory \nResearch Practice and Policy", "Psychology and Aging", "Psychology of \nAddictive Behaviors", "Psychology of Aesthetics \nCreativity and the Arts", "Psychology of Men \n& Masculinity", "Psychology of Religion \nand Spirituality", "Psychology of Violence", "Psychology Public Policy \nand Law", "Psychotherapy", "Rehabilitation Psychology", "Review of General Psychology", "School Psychology Quarterly", "Training and Education in \nProfessional Psychology")

# Filter psychology journals
psychology_journal_df2 <- filter(journal_df2, field == "Psychology")  # filter out psychology journals
psychology_journal_df2$journal <- factor(psychology_journal_df2$group, labels = psychology_journal_labels2)  # relabel journals
psychology_citation_rate_stats2 <- filter(journal_citation_rate_stats2, field == "Psychology")  # filter out psychology journals for statistics
psychology_citation_rate_stats2$journal <- psychology_journal_labels2  # relabel journals
psychology_journal_stats2 <- filter(times_cited_journal2, field == "Psychology")  # filter out psychology journals

# Plot mean citation rate per psychology journal
psychology_tau_label2 <- paste("tau ==", round(psychology_citation_rate_stats2$tau, 2))
psychology_ggplot2 <- ggplot(psychology_journal_df2, aes(x = author, y = percent_citation_rate)) +
  facet_wrap(~journal) +  # for each journal
  geom_point(size = 5, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), psychology_citation_rate_stats2, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  ylim(0, 15) +
  annotate("text", 13, 14, label = psychology_tau_label2, size = 12, parse = TRUE) +  # add tau value
  annotate("text", 13, 12, label = paste("BF =", round(psychology_citation_rate_stats2$bf, 1)), size = 12) +  # add Bayes factor
  annotate("text", 13, 10, label = paste("N =", psychology_journal_stats2$n), size = 12) +  # add N
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text.x=element_text(size=10), axis.text.y=element_text(size=20), strip.text.x = element_text(size = 15, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/psychology_journals2.png"), width = 1800, height = 2200)	# open PNG device
plot(psychology_ggplot2)		# plot figure
dev.off()		# close device

## Geoscience data
# Prepare data
geoscience_journal_labels2 <- c("Atmospheric Science Letters", "Bulletin of the Seismological \nSociety of America", "Earth Interactions", "Environmental & \nEngineering Geoscience", "Geochemistry-Exploration \nEnvironment Analysis", "Geological Society of \nAmerica Bulletin", "Geology", "Geosphere", "International Journal \nof Climatology", "Italian Journal of Geosciences", "Journal of Applied \nMeteorology and Climatology", "Journal of Atmospheric \nand Oceanic Technology", "Journal of Climate", "Journal of Hydrometeorology", "Journal of \nMicropalaeontology", "Journal of \nPhysical Oceanography", "Journal of the \nAtmospheric Sciences", "Journal of the \nGeological Society", "Lithosphere", "Meteorological Applications", "Monthly Weather Review", "Petroleum Geoscience", "Proceedings of the \nYorkshire Geological Society", "Quarterly Journal of Engineering \nGeology and Hydrogeology", "Quarterly Journal of the \nRoyal Meteorological Society", "Scottish Journal of Geology", "Seismological Research Letters", "South African Journal of Geology", "Vadose Zone Journal", "Weather and Forecasting", "Weather Climate and Society")
# Filter geoscience journals
geoscience_journal_df2 <- filter(journal_df2, field == "Geoscience")  # filter out geoscience journals
geoscience_journal_df2$journal <- factor(geoscience_journal_df2$group, labels = geoscience_journal_labels2)  # relabel journals
geoscience_citation_rate_stats2 <- filter(journal_citation_rate_stats2, field == "Geoscience")  # filter out geoscience journals for statistics
geoscience_citation_rate_stats2$journal <- geoscience_journal_labels2  # relabel journals
geoscience_journal_stats2 <- filter(times_cited_journal2, field == "Geoscience")  # filter out psychology journals

# Plot mean citation rate per geoscience journal
geoscience_tau_label2 <- paste("tau ==", round(geoscience_citation_rate_stats2$tau, 2))
geoscience_ggplot2 <- ggplot(geoscience_journal_df2, aes(x = author, y = percent_citation_rate)) +
  facet_wrap(~journal) +  # for each journal
  geom_point(size = 5, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), geoscience_citation_rate_stats2, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean citation rate") +  # create axis labels
  ylim(0, 15) +
  annotate("text", 13, 14, label = geoscience_tau_label2, size = 12, parse = TRUE) +  # add tau value
  annotate("text", 13, 12, label = paste("BF =", round(geoscience_citation_rate_stats2$bf, 1)), size = 12) +  # add Bayes factor
  annotate("text", 13, 10, label = paste("N =", geoscience_journal_stats2$n), size = 12) +  # add N
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=18), strip.text.x = element_text(size = 18, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/geoscience_journals2.png"), width = 1800, height = 2200)	# open PNG device
plot(geoscience_ggplot2)		# plot figure
dev.off()		# close device

# Follow-up analysis -----------------------------------------------

###
# Probability of citation
###
# Plot mean probability of citation per letter and citation style for data set 1
prob_style_tau_label1 <- paste("tau ==", round(citation_style_prob_cited_stats1$tau, 2))
prob_cited_citation_style_ggplot1 <- ggplot(citation_style_df1, aes(x = author, y = mean_prob_cited)) +
  facet_grid(~group) +  # for each group
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), citation_style_prob_cited_stats1, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean probability of citation") +  # create axis labels
  annotate("text", 13, 1, label = prob_style_tau_label1, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 0.985, label = paste("BF =", round(citation_style_prob_cited_stats1$bf, 2)), size = 18) +  # add Bayes factor
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=40), strip.text.x = element_text(size = 60, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/citation_styles_prob1.png"), width = 2000, height = 1000)	# open PNG device
plot(prob_cited_citation_style_ggplot1)		# plot figure
dev.off()		# close device

# Plot mean probability of citation per letter and field for data set 1
prob_field_tau_label1 <- paste("tau ==", round(field_prob_cited_stats1$tau, 2))
prob_cited_field_ggplot1 <- ggplot(field_df1, aes(x = author, y = mean_prob_cited)) +
  facet_grid(~group) +  # for each group
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), field_prob_cited_stats1, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean probability of citation") +  # create axis labels
  annotate("text", 13, 0.99, label = prob_field_tau_label1, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 0.976, label = paste("BF =", round(field_prob_cited_stats1$bf, 2)), size = 18) +  # add Bayes factor
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=40), strip.text.x = element_text(size = 60, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/fields_prob1.png"), width = 2000, height = 1000)	# open PNG device
plot(prob_cited_field_ggplot1)		# plot figure
dev.off()		# close device

# Plot mean probability of citation per letter and field for data set 2
prob_field_tau_label2 <- paste("tau ==", round(field_prob_cited_stats2$tau, 2))
prob_cited_field_ggplot2 <- ggplot(field_df2, aes(x = author, y = mean_prob_cited)) +
  facet_grid(~group) +  # for each group
  geom_point(size = 7, color = "grey50") +  # set point size
  geom_abline(aes(slope = beta, intercept = intercept), field_prob_cited_stats2, color = "black", size = 2) +  # create linear regression line
  labs(x = "First letter of author's last name",y = "Mean probability of citation") +  # create axis labels
  annotate("text", 13, 0.99, label = prob_field_tau_label2, size = 18, parse = TRUE) +  # add tau value
  annotate("text", 13, 0.971, label = paste("BF =", round(field_prob_cited_stats2$bf, 2)), size = 18) +  # add Bayes factor
  theme_bw() +  # use theme
  theme(axis.title=element_text(size=70), axis.text=element_text(size=40), strip.text.x = element_text(size = 60, margin = margin(0.5, 0, 0.5, 0, "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/fields_prob2.png"), width = 2000, height = 1000)	# open PNG device
plot(prob_cited_field_ggplot2)		# plot figure
dev.off()		# close device

# Article-level analyses -----------------------------------------------

#### 
## Prepare data
#### 
citation_data_all <- rbind(citation_data1, citation_data2)  # combine data sets
citation_data_all_nozero <- subset(citation_data_all, times_cited > 0)  # remove data with no citations
citation_data_all_nozero$field <- factor(citation_data_all_nozero$field, levels = c("Psychology", "Biology", "Geoscience"))

#### 
# Conduct LMM with log-transformed citation counts
#### 
## Full model includes interactions between letter and citation style and letter and field (but not citation style and field)
## Investigate citation style and field interactions with letter
citation_data_all_nozero$log_times_cited <- log(citation_data_all_nozero$times_cited) # log transform citation counts
letter_field_year_journal_lmer <- lmer(log_times_cited ~ author_num * field + (1 | year_num) + (1 | journal), data = citation_data_all_nozero, REML = FALSE) # fit letter and field interaction
letter_field_nointer_year_journal_lmer <- lmer(log_times_cited ~ author_num + field + (1 | year_num) + (1 | journal), data = citation_data_all_nozero, REML = FALSE) # fit letter and field main effects (no interaction)
bic_bf10(BIC(letter_field_nointer_year_journal_lmer), BIC(letter_field_year_journal_lmer))  # Calculate model comparison Bayes factors for field by letter interaction

## Assess whether random slopes should be included
letter_field_year_journal_rs.journal_lmer <- lmer(log_times_cited ~ author_num * field + (author_num | journal) + (1 | year_num), data = citation_data_all_nozero, REML = FALSE) # fit random slopes for journal
anova(letter_field_year_journal_lmer, letter_field_year_journal_rs.journal_lmer) # conduct likelihood ratio test comparing this model with model without random slopes
bic_bf10(BIC(letter_field_year_journal_lmer), BIC(letter_field_year_journal_rs.journal_lmer)) # compare Bayes factors for this model and model without random slopes
# Likelihood ratio tests suggests that including random slopes is not warranted at p = 0.112 (BF = 0.000067)

letter_field_year_journal_rs.year_lmer <- lmer(log_times_cited ~ author_num * field + (1 | journal) + (author_num | year_num), data = citation_data_all_nozero, REML = FALSE) # fit random slopes for year
anova(letter_field_year_journal_lmer, letter_field_year_journal_rs.year_lmer) # conduct likelihood ratio test comparing this model with model without random slopes
bic_bf10(BIC(letter_field_year_journal_lmer), BIC(letter_field_year_journal_rs.year_lmer)) # compare Bayes factors for this model and model without random slopes
# Likelihood ratio tests suggests that including random slopes is not warranted at p = 0.972 (BF = 0.0000077)

#### 
# Iterate LMMs across quantiles
#### 
## Prepare data
num_quantiles <- 10 # assign number of quantiles
quantile_index <- 1:num_quantiles  # create index of quantiles
field_quantiles_bf <- quant_min <- quant_max <- NA  # initiate vectors of Bayes factors and beta coefficients for citation style and field for each quantile

## Iterate LMMs for each quantile
foreach(quant = quantile_index) %do% {  # for each quantile
  # Extract top quantile data for each field
  psych_top <- extractTopQuantile(subset(citation_data_all_nozero, field == "Psychology"))
  bio_top <- extractTopQuantile(subset(citation_data_all_nozero, field == "Biology"))
  geo_top <- extractTopQuantile(subset(citation_data_all_nozero, field == "Geoscience"))
  all_top <- rbind(psych_top, bio_top)
  all_top <- rbind(all_top, geo_top)
  all_top$log_times_cited <- log(all_top$times_cited)
  # Conduct LMM for full model and model without letter by citation style interaction
  letter_field_journal_year_iter_lmer <- lmer(log_times_cited ~ author_num * field + (1 | journal) + (1 | year_num), data = all_top, REML = FALSE)
  letter_field_nointer_journal_year_iter_lmer <- lmer(log_times_cited ~ author_num + field + (1 | journal) + (1 | year_num), data = all_top, REML = FALSE)
  # Calculate Bayes factor for full/no interaction comparison
  field_quantiles_bf[quant] <- bic_bf10(BIC(letter_field_nointer_journal_year_iter_lmer), BIC(letter_field_journal_year_iter_lmer))
  # Extract minimum citation count per letter
  letter_mins <- all_top %>% group_by(author) %>% summarise(min_times_cited = min(times_cited))
  quant_min[quant] <- min(letter_mins$min_times_cited)  # find smallest minimum citation count
  quant_max[quant] <- max(letter_mins$min_times_cited)  # find largest minimum citation count
}

# Create dataframe of quantile data
quantile_data <- tibble(quantile_values = seq(from = 0, to = 90, by = 10), # create column of quantiles
  field_bf = field_quantiles_bf, # create column of Bayes factors
  min_times_cited = quant_min, # create column of minimum times cited
  max_times_cited = quant_max) %>%  # create column of maximum times cited
  mutate(ranges = paste("[", min_times_cited, ",", max_times_cited, "]", sep = ""),
    xaxis = c(0.335, 0.405, 0.475, 0.544, 0.615, 0.683, 0.753, 0.823, 0.892, 0.962)
)

## Plot BFs as a function of percentile
y_axis_breaks <- c(1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6) # create vector of y-axis breaks
# Plot figure
field_bf_plot <- ggplot(quantile_data, aes(x = quantile_values, y = field_bf)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = log10(1/3), ymax = log10(3), fill = "grey80") + # plot shaded rectangle for area with no evidence
  geom_point(size = 9) +  # plot points
  scale_y_log10(breaks = y_axis_breaks, labels = comma(y_axis_breaks)) +  # create y-axis breaks
  scale_x_continuous(breaks = seq(0, 90, 10)) + # create x-axis breaks
  geom_hline(yintercept = 3) +  # plot line for BF = 3
  geom_hline(yintercept = 10) +  # plot line for BF = 10
  geom_hline(yintercept = 30) +  # plot line for BF = 30
  geom_hline(yintercept = 100) +  # plot line for BF = 100
  geom_hline(yintercept = 1/3, linetype = 2) +  # plot line for BF = 1/3
  geom_hline(yintercept = 1/10, linetype = 2) +  # plot line for BF = 1/10
  geom_hline(yintercept = 1/30, linetype = 2) +  # plot line for BF = 1/30
  geom_hline(yintercept = 1/100, linetype = 2) +  # plot line for BF = 1/100
  # geom_text(aes(label = ranges, x = quantile_values, y = 0.000003), size = 8.5) +
  labs(x = "Percentile", y = "Bayes factor") +  # label axes
  theme_classic() +  # use theme 
  theme(axis.title.x = element_text(size = 70, margin = unit(c(1, 0, 0, 0), "lines")), axis.title.y = element_text(size = 70), axis.text=element_text(size = 50), plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "lines"))    # set font sizes and buffers
png(filename = paste("figures/field_bf.png"), width = 1200, height = 1000)	# open PNG device
plot(field_bf_plot)		# plot figure
for(i in 1:nrow(quantile_data)) {
  grid.text(quantile_data$ranges[i], x = quantile_data$xaxis[i], y = unit(3.2, "lines"), gp = gpar(fontsize = 25, col = "grey30"))
}
dev.off()		# close device

