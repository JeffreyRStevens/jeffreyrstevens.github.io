## ---
##
## Script name: miller_etal_2021_rcode.R
## 
## Manuscript: Miller et al. Individual repeatability, species differences and the influence of socio-ecological factors on neophobia in 10 corvid species
##
## Purpose of script: This script analyzes the data investigating the effect of novel objects and food on corvid neophobia.
##
## Authors: Stephan Reber (mail@stephanreber.com), Rachael Miller (neophobiaproject@gmail.com), and Jeffrey R. Stevens (jeffrey.r.stevens@gmail.com)
##
## Date Finalized: 2021-10-13
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
## Data file:
## ---
## miller_etal_2021_data.csv
##
## individual - Individual subject name
## species - Species common name
## lab - Experimental lab
## sex - Subject sex (m = male, f = female)
## age - Subject age (Adult or Juvenile)
## condition - Experimental condition
## round - Round number (1-3)
## range - Species range (0 = mainland, 1 = island)
## urban_habitat - Usage of urban habitats (0 = no, 1 = yes)
## adult_sociality - Level of sociality for adults (0 = territorial pairs, 1 = family groups)
## caching - Level of caching (0 = moderate, 1 = specialized)
## source - Source of subject (0 = captivity, 1 = wild)
## live_hunting - Exhibits live hunting (0 = no, 1 = yes)
## flock_size - Average flock size (0 = less than 100, 1 = greater than 100)
## touch - Whether subject touched familar food in 10 minutes (0 = no, 1 = yes)
## latency - Latency to touch familiar food (s)
##
## ---

# Load packages ----------------------------------------------------------------

library(multcomp)
library(lme4)
library(car)
library(ggsignif)
library(ggbeeswarm)
library(patchwork)
library(tidyverse)



# Import and prepare data -------------------------------------------------------------

all_data <- read.csv("miller_etal_2021_data.csv")
all_data <- all_data %>% 
  mutate(across(individual:condition, factor),
         across(range:flock_size, factor),
         individual = str_trim(individual),
         condition = fct_recode(condition, "Novel food" = "Novel_food", "Novel object" = "Novel_object"),
         round = factor(round, ordered = TRUE),
         range = fct_recode(range, "Mainland" = "0", "Island" = "1"),
         range = fct_relevel(range, c("Mainland", "Island")),
         urban_habitat = fct_recode(urban_habitat, "No" = "0", "Yes" = "1"),
         urban_habitat = fct_relevel(urban_habitat, c("No", "Yes")),
         adult_sociality = fct_recode(adult_sociality, "Territorial \npairs" = "0", "Family \ngroups" = "1"),
         adult_sociality = fct_relevel(adult_sociality, c("Territorial \npairs", "Family \ngroups")),
         caching = fct_recode(caching, "Moderate \ncacher" = "0", "Specialised \ncacher" = "1"),
         caching = fct_relevel(caching, c("Moderate \ncacher", "Specialised \ncacher")),
         live_hunting = fct_recode(live_hunting, "No" = "0", "Yes" = "1"),
         live_hunting = fct_relevel(live_hunting, c("No", "Yes")),
         flock_size = fct_recode(flock_size, "Small" = "0", "Large" = "1"),
         flock_size = fct_relevel(flock_size, c("Small", "Large")),
         source = fct_recode(source, "Captive" = "0", "Wild" = "1"),
         source = fct_relevel(source, c("Captive", "Wild")),
         species = fct_recode(species, "Large-billed crow" = "large_billed_crow", "Clark's nutcracker" = "clarks_nutcracker", "Blue jay" = "blue_jay", "New Caledonian crow" = "New_caledonian_crow", "Azure-winged magpie" = "azure_winged_magpie", "Pinyon jay" = "pinyon_jay", "Carrion crow" = "carrion_crow", "Eurasian jay" = "eurasian_jay", "Common raven" = "common_raven", "'Alalā" = "'alalā"),
         genus = factor(ifelse(species %in% c("Large-billed crow", "New Caledonian crow", "Carrion crow", "Common raven", "'Alalā"), "Corvus", "Other")),
         loglatency = log10(latency + 1)
  ) %>% 
  relocate(genus, .after = species)

adult_data <- all_data %>% 
  filter(age == "Adult")

individual_demo_data <- all_data %>% 
  select(-c(latency, loglatency, condition)) %>% 
  group_by(individual, round) %>% 
  slice_head(n = 1)

diff_scores_log <- all_data %>% 
  unite(cond_round, condition, round) %>% 
  pivot_wider(id_cols = individual, names_from = cond_round, values_from = loglatency) %>% 
  mutate(novelobjectdiff_1 = `Novel object_1` - Control_1,
         novelobjectdiff_2 = `Novel object_2` - Control_2,
         novelobjectdiff_3 = `Novel object_3` - Control_3,
         novelfooddiff_1 = `Novel food_1` - Control_1,
         novelfooddiff_2 = `Novel food_2` - Control_2,
         novelfooddiff_3 = `Novel food_3` - Control_3) %>% 
  select(individual, contains("diff")) %>% 
  pivot_longer(contains("diff"), names_to = "cond_round", values_to = "diff_score") %>% 
  separate(cond_round, c("condition", "round")) %>% 
  filter(!is.na(diff_score)) %>% 
  rename(diff_score_log = diff_score)

diff_scores <- all_data %>% 
  unite(cond_round, condition, round) %>% 
  pivot_wider(id_cols = individual, names_from = cond_round, values_from = latency) %>% 
  mutate(novelobjectdiff_1 = `Novel object_1` - Control_1,
         novelobjectdiff_2 = `Novel object_2` - Control_2,
         novelobjectdiff_3 = `Novel object_3` - Control_3,
         novelfooddiff_1 = `Novel food_1` - Control_1,
         novelfooddiff_2 = `Novel food_2` - Control_2,
         novelfooddiff_3 = `Novel food_3` - Control_3) %>% 
  select(individual, contains("diff")) %>% 
  pivot_longer(contains("diff"), names_to = "cond_round", values_to = "diff_score") %>% 
  separate(cond_round, c("condition", "round")) %>% 
  filter(!is.na(diff_score)) %>% 
  left_join(individual_demo_data, by = c("individual", "round")) %>% 
  left_join(diff_scores_log, by = c("individual", "condition", "round")) %>% 
  select(individual, species, genus, lab, sex, age, round, condition, diff_score, diff_score_log, range:flock_size)

novelobject_control_diff <- diff_scores %>% 
  filter(condition == "novelobjectdiff")

novelobject_control_diff_adult <- novelobject_control_diff %>% 
  filter(age == "Adult")

novelfood_control_diff <- diff_scores %>% 
  filter(condition == "novelfooddiff")

novelfood_control_diff_adult <- novelfood_control_diff %>% 
  filter(age == "Adult")

individual_condition_data <- all_data %>% 
  group_by(individual, species, genus, lab, condition, range, urban_habitat, adult_sociality, caching, live_hunting, source, flock_size) %>% 
  summarise(mean_latency = mean(latency), log_latency = mean(loglatency), .groups = "drop")

individual_round_data <- all_data %>% 
  group_by(individual, species, genus, lab, round, range, urban_habitat, adult_sociality, caching, live_hunting, source, flock_size) %>% 
  summarise(mean_latency = mean(latency), log_latency = mean(loglatency), .groups = "drop")

individual_data <- all_data %>% 
  group_by(individual, species, genus, lab, range, urban_habitat, adult_sociality, caching, live_hunting, source, flock_size) %>% 
  summarise(mean_latency = mean(latency), log_latency = mean(loglatency), .groups = "drop")

individual_diff_data <- diff_scores %>% 
  group_by(individual, species, genus, range, urban_habitat, adult_sociality, caching, live_hunting, source, flock_size) %>% 
  summarise(mean_latency = mean(diff_score), log_latency = mean(diff_score_log), .groups = "drop") 

individual_novelobjectdiff_data <- novelobject_control_diff %>% 
  group_by(individual, species, genus, range, urban_habitat, adult_sociality, caching, live_hunting, source, flock_size) %>% 
  summarise(mean_latency = mean(diff_score), log_latency = mean(diff_score_log), .groups = "drop") 

individual_novelfooddiff_data <- novelfood_control_diff %>% 
  group_by(individual, species, genus, range, urban_habitat, adult_sociality, caching, live_hunting, source, flock_size) %>% 
  summarise(mean_latency = mean(diff_score), log_latency = mean(diff_score_log), .groups = "drop") 


# Q1 Condition, round, species ----------------------------------------------------------------------

## Raw scores --------------------------------------------------------------

# LMM 1 without condition x species interaction
latency_cond_species_round_lm0 <- lmer(loglatency ~ condition + species + round + (1 | individual), data = all_data)
latency_cond_species_round_lm1 <- lmer(loglatency ~ condition + species + round + (1 | lab:individual), data = all_data)
latency_cond_species_round_lm2 <- lmer(loglatency ~ condition + species + round + (1 | lab/individual), data = all_data)
latency_cond_species_round_lm3 <- lmer(loglatency ~ condition + species + round + (1 | lab) + (1 | individual), data = all_data)
anova(latency_cond_species_round_lm0, latency_cond_species_round_lm1, latency_cond_species_round_lm2, latency_cond_species_round_lm3)

summary(latency_cond_species_round_lm0)
drop1(latency_cond_species_round_lm0, test = "Chisq")

# Tukey contrasts
latency_cond_glht <- summary(glht(latency_cond_species_round_lm0, linfct = mcp(condition = "Tukey")), test = adjusted("none"))

latency_species_glht <- summary(glht(latency_cond_species_round_lm0, linfct = mcp(species = "Tukey")), test = adjusted("none"))

latency_round_glht <- summary(glht(latency_cond_species_round_lm0, linfct = mcp(round = "Tukey")), test = adjusted("none"))

# Adult data
latency_cond_species_round_adult_lm0 <- lmer(loglatency ~ condition + species + round + (1 | individual), data = adult_data)

summary(latency_cond_species_round_adult_lm0)
drop1(latency_cond_species_round_adult_lm0, test = "Chisq")

latency_cond_adult_glht <- summary(glht(latency_cond_species_round_adult_lm0, linfct = mcp(condition = "Tukey")), test = adjusted("none"))

latency_species_adult_glht <- summary(glht(latency_cond_species_round_adult_lm0, linfct = mcp(species = "Tukey")), test = adjusted("none"))

latency_round_adult_glht <- summary(glht(latency_cond_species_round_adult_lm0, linfct = mcp(round = "Tukey")), test = adjusted("none"))

## Difference score novel object --------------------------------------------------------

# LMM 2, difference score novel object, all variables
diffscore_species_round_novelobject_lm0 <- lmer(diff_score_log ~ species + round + (1 | individual), data = novelobject_control_diff)
diffscore_species_round_novelobject_lm1 <- lmer(diff_score_log ~ species + round + (1 | lab:individual), data = novelobject_control_diff)
diffscore_species_round_novelobject_lm2 <- lmer(diff_score_log ~ species + round + (1 | lab/individual), data = novelobject_control_diff)
diffscore_species_round_novelobject_lm3 <- lmer(diff_score_log ~ species + round + (1 | lab) + (1 | individual), data = novelobject_control_diff)
anova(diffscore_species_round_novelobject_lm0, diffscore_species_round_novelobject_lm1, diffscore_species_round_novelobject_lm2, diffscore_species_round_novelobject_lm3)

summary(diffscore_species_round_novelobject_lm0)
drop1(diffscore_species_round_novelobject_lm0, test = "Chisq")

# Tukey contrasts
diffscore_species_novelobject_glht <- summary(glht(diffscore_species_round_novelobject_lm0, linfct = mcp(species = "Tukey")), test = adjusted("none"))

# Adult data
diffscore_species_round_novelobject_adult_lm0 <- lmer(diff_score_log ~ species + round + (1 | individual), data = novelobject_control_diff_adult)

summary(diffscore_species_round_novelobject_adult_lm0)
drop1(diffscore_species_round_novelobject_adult_lm0, test = "Chisq")

diffscore_species_novelobject_adult_glht <- summary(glht(diffscore_species_round_novelobject_adult_lm0, linfct = mcp(species = "Tukey")), test = adjusted("none"))

## Difference score novel food --------------------------------------------------------

# LMM 3, difference score novel food, all variables
diffscore_species_round_novelfood_lm0 <- lmer(diff_score_log ~ species + round + (1 | individual), data = novelfood_control_diff)
diffscore_species_round_novelfood_lm1 <- lmer(diff_score_log ~ species + round + (1 | lab:individual), data = novelfood_control_diff)
diffscore_species_round_novelfood_lm2 <- lmer(diff_score_log ~ species + round + (1 | lab/individual), data = novelfood_control_diff)
anova(diffscore_species_round_novelfood_lm0, diffscore_species_round_novelfood_lm1, diffscore_species_round_novelfood_lm2)

summary(diffscore_species_round_novelfood_lm0)
drop1(diffscore_species_round_novelfood_lm0, test = "Chisq")

# Tukey contrasts
diffscore_species_novelfood_glht <- summary(glht(diffscore_species_round_novelfood_lm0, linfct = mcp(species = "Tukey")), test = adjusted("none"))

# Adult data
diffscore_species_round_novelfood_adult_lm0 <- lmer(diff_score_log ~ species + round + (1 | individual), data = novelfood_control_diff_adult)

summary(diffscore_species_round_novelfood_adult_lm0)
drop1(diffscore_species_round_novelfood_adult_lm0, test = "Chisq")

diffscore_species_novelfood_adult_glht <- summary(glht(diffscore_species_round_novelfood_adult_lm0, linfct = mcp(species = "Tukey")), test = adjusted("none"))


## Lab effects ------------------------------------------------------------

### LMM for lab effects i.e. 3 species tested at more than 1 lab

# LMM a, novel object, pinyon jay, lab effect
latency_round_side_pinyonjay_lm <- lmer(loglatency ~ lab + (1 | individual), data = filter(all_data, species == "Pinyon jay"))

(lab_pinyonjay <- drop1(latency_round_side_pinyonjay_lm, test = "Chisq"))
names(lab_pinyonjay)[4] <- "p"
Anova(latency_round_side_pinyonjay_lm)

# LMM b, novel object, carrion crow, lab effect
latency_round_side_carrioncrow_lm <- lmer(loglatency ~ lab + (1 | individual), data = filter(all_data, species == "Carrion crow"))

(lab_carrioncrow <- drop1(latency_round_side_carrioncrow_lm, test = "Chisq"))
names(lab_carrioncrow)[4] <- "p"
Anova(latency_round_side_carrioncrow_lm)

# LMM c, novel object, azure-winged magpie, lab effect
latency_round_side_azurewingedmagpie_lm <- lmer(loglatency ~ lab + (1 | individual), data = filter(all_data, species == "Azure-winged magpie"))

(lab_azurewingedmagpie <- drop1(latency_round_side_azurewingedmagpie_lm, test = "Chisq"))
names(lab_azurewingedmagpie)[4] <- "p"
Anova(latency_round_side_azurewingedmagpie_lm)


# Q2 Socio-ecology ----------------------------------------------------------------------

## Novel object ----------------------------------------------------------------------

# LMM 4, difference score data, socio-ecological factors
diffscore_socioecology_novelobject_lm0 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | individual), data = novelobject_control_diff)
diffscore_socioecology_novelobject_lm1 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | lab:individual), data = novelobject_control_diff)
diffscore_socioecology_novelobject_lm2 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | species:individual), data = novelobject_control_diff)
diffscore_socioecology_novelobject_lm3 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | species:lab:individual), data = novelobject_control_diff)
anova(diffscore_socioecology_novelobject_lm0, diffscore_socioecology_novelobject_lm1, diffscore_socioecology_novelobject_lm2, diffscore_socioecology_novelobject_lm3)

summary(diffscore_socioecology_novelobject_lm0)
drop1(diffscore_socioecology_novelobject_lm0, test = "Chisq")

# Adult data

diffscore_socioecology_novelobject_adult_lm0 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | individual), data = novelobject_control_diff_adult)

summary(diffscore_socioecology_novelobject_adult_lm0)
drop1(diffscore_socioecology_novelobject_adult_lm0, test = "Chisq")


## Novel food ----------------------------------------------------------------------

# LMM 5, difference score data, socio-ecological factors
diffscore_socioecology_novelfood_lm0 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | individual), data = novelfood_control_diff)
diffscore_socioecology_novelfood_lm1 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | lab:individual), data = novelfood_control_diff)
diffscore_socioecology_novelfood_lm2 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | species:individual), data = novelfood_control_diff)
diffscore_socioecology_novelfood_lm3 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | species:lab:individual), data = novelfood_control_diff)
anova(diffscore_socioecology_novelfood_lm0, diffscore_socioecology_novelfood_lm1, diffscore_socioecology_novelfood_lm2, diffscore_socioecology_novelfood_lm3)

summary(diffscore_socioecology_novelfood_lm0)
drop1(diffscore_socioecology_novelfood_lm0, test = "Chisq")

# Adult data

diffscore_socioecology_novelfood_adult_lm0 <- lmer(diff_score_log ~ range + urban_habitat + adult_sociality + flock_size + caching + live_hunting + genus + (1 | individual), data = novelfood_control_diff_adult)

summary(diffscore_socioecology_novelfood_adult_lm0)
drop1(diffscore_socioecology_novelfood_adult_lm0, test = "Chisq")


# Q3 Temporal and contextual repeatability ----------------------------------------------------------------------

# Condition

condition_raw_rep <- tibble(.rows = 0)

for(this_condition in unique(all_data$condition)) {
  condition_long <- all_data %>% 
    filter(condition == this_condition) %>% 
    select(individual, species, round, loglatency)
  this_rpt <- rptR::rpt(loglatency ~ species + round + (1 | individual), data = condition_long, grname = "individual", parallel = TRUE)
  this_row <- c(this_condition, this_rpt$ngroups, round(this_rpt$R$individual, 3), round(this_rpt$CI_emp$`2.5%`, 3), round(this_rpt$CI_emp$`97.5%`, 3), round(this_rpt$P$LRT_P, 3))
  names(this_row) <- c("condition", "n", "icc", "lowerCI", "upperCI", "p_value")
  condition_raw_rep <- bind_rows(condition_raw_rep, this_row)
}

# Novel object & food

novel_data <- all_data %>% 
  filter(condition != "Control")
novel_rpt <- rptR::rpt(loglatency ~ condition + species + round + (1 | individual), data = novel_data, grname = "individual", parallel = TRUE)
novel_row <- c("Overall repeatability", novel_rpt$ngroups, round(novel_rpt$R$individual, 3), round(novel_rpt$CI_emp$`2.5%`, 3), round(novel_rpt$CI_emp$`97.5%`, 3), round(novel_rpt$P$LRT_P, 3))
names(novel_row) <- c("condition", "n", "icc", "lowerCI", "upperCI", "p_value")

# Overall

raw_rpt <- rptR::rpt(loglatency ~ condition + species + round + (1 | individual), data = all_data, grname = "individual", parallel = TRUE)
raw_row <- c("Overall repeatability", raw_rpt$ngroups, round(raw_rpt$R$individual, 3), round(raw_rpt$CI_emp$`2.5%`, 3), round(raw_rpt$CI_emp$`97.5%`, 3), round(raw_rpt$P$LRT_P, 3))
names(raw_row) <- c("condition", "n", "icc", "lowerCI", "upperCI", "p_value")

# Species by condition 

set.seed(2021)

condition_species_raw_rep <- tibble(.rows = 0)

for (this_species in unique(all_data$species)) {
  for (this_condition in unique(all_data$condition)) {
    if (this_species == "'Alalā" & this_condition == "Novel food") {}
    else {
      condition_species_long <- all_data %>% 
        filter(condition == this_condition & species == this_species) %>% 
        select(individual, round, loglatency)
      this_rpt <- rptR::rpt(loglatency ~ round + (1 | individual), data = condition_species_long, grname = "individual", parallel = TRUE)
      this_row <- c(this_species, this_condition, this_rpt$ngroups, round(this_rpt$R$individual, 3), round(this_rpt$CI_emp$`2.5%`, 3), round(this_rpt$CI_emp$`97.5%`, 3), round(this_rpt$P$LRT_P, 3))
      names(this_row) <- c("species", "condition", "n", "icc", "lowerCI", "upperCI", "p_value")
      condition_species_raw_rep <- bind_rows(condition_species_raw_rep, this_row)
    }
  }
}

# Species

species_raw_rep <- tibble(.rows = 0)

for(this_species in unique(all_data$species)) {
  species_long <- all_data %>% 
    filter(species == this_species) %>% 
    select(individual, condition, round, loglatency) %>% 
    group_by(individual, condition) %>% 
    summarise(mean_latency = mean(loglatency, na.rm = TRUE), .groups = "drop")
  this_rpt <- rptR::rpt(mean_latency ~ condition + (1 | individual), data = species_long, grname = "individual", parallel = TRUE)
  this_row <- c(this_species, this_rpt$ngroups, round(this_rpt$R$individual, 3), round(this_rpt$CI_emp$`2.5%`, 3), round(this_rpt$CI_emp$`97.5%`, 3), round(this_rpt$P$LRT_P, 3))
  names(this_row) <- c("species", "n", "icc", "lowerCI", "upperCI", "p_value")
  species_raw_rep <- bind_rows(species_raw_rep, this_row)
}


# Plots -------------------------------------------------------------------

## Condition plot ----------------------------------------------------------

individual_condition_data %>% 
  count(condition)

individual_condition_plot <- individual_condition_data %>% 
  ggplot(aes(x = condition, y = mean_latency)) +
  geom_beeswarm(cex = 0.5, size = 1, shape = 21, fill = 'white', color = "grey40", alpha = 0.5) +
  stat_summary(fun = median, shape = "-", size = 4) +
  labs(x = "Condition", y = "Latency to touch food (s)") +
  geom_signif(comparisons = list(c("Novel food", "Control"), c("Control", "Novel object"), c("Novel food", "Novel object")), step_increase = c(0.1, 0.1, 0), annotation = "***", textsize = 6, vjust = 0.5) +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 18))


## Round plot --------------------------------------------------------------

individual_round_plot <- individual_round_data %>% 
  ggplot(aes(x = round, y = mean_latency)) +
  geom_beeswarm(cex = 0.75, size = 1, shape = 21, fill = 'white', color = "grey40", alpha = 0.5) +
  stat_summary(fun = median, shape = "-", size = 4) +
  labs(x = "Round", y = "Latency to touch food (s)") +
  geom_signif(comparisons = list(c("2", "3"), c("1", "3")), step_increase = 0.1, annotation = "***", textsize = 6, vjust = 0.5) +
  geom_signif(y_position = 680, comparisons = list(), annotation = "***", textsize = 6, vjust = 0.5, tip_length = 0) +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 18))

individual_condition_plot + individual_round_plot + plot_annotation(tag_levels = "A")
ggsave("figures/latency_condition_round.png", width = 10, height = 5)


## Species plot ------------------------------------------------------------

individual_condition_data %>% 
  mutate(species = fct_relevel(species, "Eurasian jay", "Blue jay", "Pinyon jay", "Clark's nutcracker", "Common raven", "'Alalā", "Large-billed crow", "New Caledonian crow", "Carrion crow", "Azure-winged magpie")) %>% 
  ggplot(aes(y = mean_latency, x = species, group = condition, color = condition, shape = condition)) +
  geom_point(alpha = 0.25, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, position = position_dodge(width = 0.5), size = 0.35) +
  labs(y = "Latency to touch food (s)", x = "Species") +
  coord_flip() +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 18),
        legend.position = c(0.85, 0.79),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))
ggsave("figures/latency_species_phylogenetic.png", width = 8, height = 5)

diffscore_species_novelobject_plot <- individual_novelobjectdiff_data %>% 
  mutate(species = fct_relevel(species, "Eurasian jay", "Blue jay", "Pinyon jay", "Clark's nutcracker", "Common raven", "'Alalā", "Large-billed crow", "New Caledonian crow", "Carrion crow", "Azure-winged magpie")) %>% 
  ggplot(aes(y = mean_latency, x = species)) +
  geom_jitter(width = 0.1, color = "grey", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  stat_summary(fun.data = mean_cl_normal) +
  labs(y = "Mean latency difference score (s)", x = "Species") +
  coord_flip() +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 18))

diffscore_species_novelfood_plot <- individual_novelfooddiff_data %>% 
  mutate(species = fct_relevel(species, "Eurasian jay", "Blue jay", "Pinyon jay", "Clark's nutcracker", "Common raven", "'Alalā", "Large-billed crow", "New Caledonian crow", "Carrion crow", "Azure-winged magpie")) %>% 
  ggplot(aes(y = mean_latency, x = species)) +
  geom_jitter(width = 0.1, color = "grey", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  stat_summary(fun.data = mean_cl_normal) +
  labs(y = "Mean latency difference score (s)", x = "Species") +
  coord_flip() +
  theme_classic(base_family = "Arial") +
  theme(text = element_text(size = 18))

diffscore_species_novelobject_plot + diffscore_species_novelfood_plot + plot_layout(nrow = 2) + plot_annotation(tag_levels = "A")
ggsave("figures/differencescore_species_phylogenetic.png", width = 8, height = 8)



## Socio-ecology novel object -----------------------------------------------------

individual_novelobjectdiff_data_long <- individual_novelobjectdiff_data %>% 
  select(-source) %>% 
  pivot_longer(genus:flock_size, names_to = "socioecology", values_to = "values") %>% 
  mutate(socioecology = fct_recode(socioecology, "Adult sociality" = "adult_sociality", "Caching" = "caching", "Genus" = "genus", "Flock size" = "flock_size", "Live hunting" = "live_hunting", "Range" = "range", "Urban habitat" = "urban_habitat"),
         socioecology = fct_relevel(socioecology, c("Range", "Urban habitat", "Adult sociality", "Flock size", "Caching", "Live hunting", "Genus")))

diffscore_socioecology_novelobject_signif <- data.frame(x = rep(1, 4), y = rep(650, 4), xend = rep(2, 4), yend = rep(650, 4), socioecology = c("Adult sociality", "Caching", "Flock size", "Urban habitat"), label = c("***", "*", "*", "**")) %>% 
  mutate(socioecology = factor(socioecology, levels = c("Range", "Urban habitat", "Adult sociality", "Flock size", "Caching", "Live hunting", "Genus")))

diffscore_socioecology_novelobject_plot <- individual_novelobjectdiff_data_long %>% 
  ggplot(aes(x = values, y = mean_latency)) +
  geom_beeswarm(color = "grey", shape = "circle open") +
  stat_summary(fun = median, shape = "-", size = 5) +
  facet_wrap(~ socioecology, nrow = 2, scales = "free_x") +
  geom_segment(data = diffscore_socioecology_novelobject_signif,
               aes(x = x, y = y, yend = yend, xend = xend), inherit.aes = FALSE) +
  geom_text(data = diffscore_socioecology_novelobject_signif,
            aes(x = 1.5, y = 675, label = label), inherit.aes = FALSE, size = 5) +
  ylim(-200, 700) +
  labs(x = "", y = "Mean latency difference score (s)") +
  theme_bw(base_family = "Arial") +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())


## Socio-ecology novel food -----------------------------------------------------

individual_novelfooddiff_data_long <- individual_novelfooddiff_data %>% 
  select(-source) %>% 
  pivot_longer(genus:flock_size, names_to = "socioecology", values_to = "values") %>% 
  mutate(socioecology = fct_recode(socioecology, "Adult sociality" = "adult_sociality", "Caching" = "caching", "Genus" = "genus", "Flock size" = "flock_size", "Live hunting" = "live_hunting", "Range" = "range", "Urban habitat" = "urban_habitat"),
         socioecology = fct_relevel(socioecology, c("Range", "Urban habitat", "Adult sociality", "Flock size", "Caching", "Live hunting", "Genus")))

diffscore_socioecology_novelfood_signif <- data.frame(x = rep(1, 1), y = rep(650, 1), xend = rep(2, 1), yend = rep(650, 1), socioecology = c("Flock size"), label = c("**")) %>% 
  mutate(socioecology = factor(socioecology, levels = c("Range", "Urban habitat", "Adult sociality", "Flock size", "Caching", "Live hunting", "Genus")))

diffscore_socioecology_novelfood_plot <- individual_novelfooddiff_data_long %>% 
  ggplot(aes(x = values, y = mean_latency)) +
  geom_beeswarm(color = "grey", shape = "circle open") +
  stat_summary(fun = median, shape = "-", size = 5) +
  facet_wrap(~ socioecology, nrow = 2, scales = "free_x")  +
  geom_segment(data = diffscore_socioecology_novelfood_signif,
               aes(x = x, y = y, yend = yend, xend = xend), inherit.aes = FALSE) +
  geom_text(data = diffscore_socioecology_novelfood_signif,
            aes(x = 1.5, y = 675, label = label), inherit.aes = FALSE, size = 5) +
  ylim(-400, 700) +
  labs(x = "", y = "Mean latency difference score (s)") +
  theme_bw(base_family = "Arial") +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

diffscore_socioecology_novelobject_plot  + diffscore_socioecology_novelfood_plot + plot_annotation(tag_levels = "A") + plot_layout(nrow = 2)
ggsave("figures/differencescore_socioecology.png", width = 10, height = 12)



## Lab effects ------------------------------------------------------------

lab_data <- individual_data %>% 
  filter(species %in% c("Pinyon jay", "Carrion crow", "Azure-winged magpie")) %>% 
  mutate(lab = factor(lab),
         lab = fct_relevel(lab, "1", "6", "7", "8", "9")) 

lab_signif <- data.frame(x = rep(1, 1), y = rep(650, 1), xend = rep(2, 1), yend = rep(650, 1), species = c("Pinyon jay"), label = c("*"))

lab_data %>% 
  ggplot(aes(x = lab, y = mean_latency)) +
  geom_beeswarm(cex = 2, color = "grey") +
  stat_summary(fun = median, shape = "-", size = 5) +
  facet_wrap(~ species, scales = "free_x") +
  geom_segment(data = lab_signif,
               aes(x = x,y = y,yend = yend,xend = xend), inherit.aes = FALSE) +
  geom_text(data = lab_signif,
            aes(x = 1.5,y = 665,label = label), inherit.aes = FALSE, size = 5) +
  ylim(0, 670) +
  labs(x = "Lab", y = "Latency to touch food (s)") +
  theme_bw(base_family = "Arial") +
  theme(text = element_text(size = 18),
        panel.grid.major.x = element_blank())
ggsave("figures/latency_lab.png", width = 8, height = 5)


