#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for  analyzing ABCD data release 5.1 for the paper entitled                  #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Nov 19th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, summarytools)

options(scipen = 999)

#2# Set paths

datadir <- '/Users/jongonzalez/Desktop/abcd.nosync/'

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

# table for report
vtable::st(table1 <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_original.csv'), stringsAsFactors = T) %>% 
             mutate_at(vars(alcohol, tobacco, marihuana, opioids, opiates, cocaine), 
                       ~factor(ifelse(. == 1, 'Yes', ifelse(. == 0, 'No', ifelse(. == 999, 'Do not know', 'Refuse to answer'))))) %>% 
             mutate(bfq_duration = factor(case_when(bfq_duration == 0 ~ 'No breastfed',
                                                    bfq_duration == 1 ~ '1-6 months',
                                                    bfq_duration == 2 ~ '7-12 months',
                                                    bfq_duration == 3 ~ '>12 months'), 
                                          levels = c('No breastfed', '1-6 months', '7-12 months', '>12 months')),
                    pregnancy_complications = factor(ifelse(pregnancy_complications == 0, 'No', ifelse(pregnancy_complications == 1, 'Yes', 'Do not know'))),
                    birth_complications = factor(ifelse(birth_complications == 0, 'No', ifelse(birth_complications == 1, 'Yes', 'Do not know'))),
                    premature_wk = na_if(premature_wk, 0),
                    eventname = relevel(eventname, ref = 'baseline_year_1_arm_1')),
           digits = 3,
           vars = c('child_age', 'child_sex', 'pds', 'handedness', 'race', 'ethnicity', 'mode_of_delivery', 'birth_complications', 'pregnancy_complications', 
           'premature', 'premature_wk', 'education3', 'annual_income', 'alcohol', 'tobacco', 'marihuana', 'opioids', 'opiates', 'cocaine', 'bfq_duration'), group = 'eventname')

####################################
####################################
#4# Load CT data for main analysis #
####################################
####################################

#3# Set common covariates

(covs2scale <- c('bfq_duration', 'child_age', 'education2', 'smri_vol_scs_intracranialv')) # 'bfq_duration' --> watch out with mean-centering the main predictor; you want the intercept to reflect when bf is 0, if you mean-center then you shift the 0 to the mean of bf, which is 1.43

brain <- 
  # fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_original.csv'), stringsAsFactors = T) %>%               # ORIGINAL data (only combatt'd)
  fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_no_outliers.csv'), stringsAsFactors = T) %>%              # no outliers (NAs) and those with NAs at baselien but not follow-up set to NAs at both timepoints
  # fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_NO_COMBAT.csv'), stringsAsFactors = T) %>%             # regular data (no combatt'd)
  mutate_at(vars(all_of(covs2scale)), ~scale(., scale = T)) %>%  # mean-center continous covariates
  mutate(breastfed = as.factor(ifelse(bfq_duration > 0, 'Yes', 'No')),
         handedness = factor(handedness),
         child_sex = ifelse(child_sex == 'Male', 1, -1),   # effects code sex
         premature = ifelse(premature == 'Yes', 1, -1),    # effects code premature
         interaction = child_age * bfq_duration)

contrasts(brain$handedness) <- contr.sum(levels(brain$handedness))

ct.df <- brain %>% 
  select(., !ends_with(c('_myelin', '_area', '_ratio'))) %>% 
  as.data.table()

# check multicollinearity function

calculate_vif <- function(data) {
  # Ensure the data is numeric
  numeric_data <- as.data.frame(lapply(data, as.numeric))
  # Compute the correlation matrix
  cor_matrix <- cor(numeric_data)
  # Compute the inverse of the correlation matrix
  inv_cor_matrix <- solve(cor_matrix)
  # Calculate VIF for each predictor
  vif_values <- diag(inv_cor_matrix)
  # Name the VIF values based on predictor names
  names(vif_values) <- colnames(numeric_data)
  return(vif_values)
}

# Select only the predictors you want to analyze for VIF
predictors <- ct.df[, c('child_age', 'child_sex', 'premature', 'bfq_duration', 'education2', 'interaction', 'smri_vol_scs_intracranialv')]

# Calculate VIFs
(vif_results <- calculate_vif(predictors))

# Set the DVs
(featurenames <- names(ct.df %>% select(starts_with('mrisdp'))))

cth.model_list <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + premature + education2 + child_sex + (1 | src_subject_id)")  
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = ct.df)  
  cth.model_list[[feature]] <- model
}

posthoc1 <- as.data.frame(matrix(NA, nrow = length(featurenames), ncol = 8))    # we will store all estimates and important values in this table

for (i in 1:length(cth.model_list)) {
  posthoc1[,1][[i]] <- names(cth.model_list[i])
  posthoc1[,2][[i]] <- length(!is.na(resid(cth.model_list[[i]])))
  # calculate beta estimate and CI (standardized)
  ci.tmp <- effectsize::standardize_parameters(cth.model_list[[i]], method = 'basic')
  # extract beta and CI
  posthoc1[,3][[i]] <- ci.tmp[[2]][3]    # 3 main effect, 7 interaction
  posthoc1[,4][[i]] <- ci.tmp[[4]][3]
  posthoc1[,5][[i]] <- ci.tmp[[5]][3]
  # raw p-value
  posthoc1[,6][[i]] <- round(summary(cth.model_list[[i]])[[10]][3,5], 5)
  if (sum(is.na(posthoc1[,6]) < length(posthoc1[,6]))) {
    names(posthoc1) <- c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr', 'fdr2')
    posthoc1[,7] <- round(p.adjust(posthoc1[,6], method = 'BH'), 5) 
    posthoc1[,8] <- round(mutoss::adaptiveBH(posthoc1[,6], 0.05, silent = T)$adjPValues, 5) }
    # posthoc1[,8] <- round(cp4p::adjust.p(posthoc1[,6], pi0.method = 'bky', alpha = 0.05)$adjp[,2], 5) }
}

print(posthoc1)

print(subset(posthoc1, posthoc1$`p-value` < 0.05))    # regions p < 0.05

sum(posthoc1$`p-value` < 0.05)                       # number of regions p < 0.05

print(subset(posthoc1, posthoc1$fdr < 0.05))        # regions fdr < 0.05

sum(posthoc1$fdr < 0.05)                          # number of regions fdr < 0.05

(passed.fdr.cth <- subset(posthoc1$ROI, posthoc1$fdr < 0.05))    # save the regions that fdr < 0.05

# this is to calculate the percentage change and for that we need the unstandardized coefficients

significant.models.list <- lapply(passed.fdr.cth, function(model_name) {
  cth.model_list[[model_name]]
})

dv <- c()
raw.beta <- c()
intercept <- c() 

for (i in 1:length(passed.fdr.cth)) {
  dv[i] <- passed.fdr.cth[i]
  raw.beta[i] <- summary(significant.models.list[[i]])[[10]][3,1]
  intercept[i] <- summary(significant.models.list[[i]])[[10]][1,1]
}

(ct.change <- cbind(raw.beta, intercept) %>% 
  as.data.frame() %>% 
  mutate(change = (raw.beta / intercept) * 100))

mean(ct.change$raw.beta)     # for each increase in breastfeeding duration there was a +0.006mm increase in thickness

mean(ct.change$change)       # this equals a +0.027% change

# Look interactions now

interaction.cth <- as.data.frame(matrix(NA, nrow = length(featurenames), ncol = 8))

for (i in 1:length(cth.model_list)) {
  interaction.cth[,1][[i]] <- names(cth.model_list[i])
  interaction.cth[,2][[i]] <- length(!is.na(resid(cth.model_list[[i]])))
  # confidence intervals
  ci.tmp <- effectsize::standardize_parameters(cth.model_list[[i]], method = 'basic')
  interaction.cth[,3][[i]] <- ci.tmp[[2]][7]    # 3 main effect, 7 interaction
  interaction.cth[,4][[i]] <- ci.tmp[[4]][7]
  interaction.cth[,5][[i]] <- ci.tmp[[5]][7]
  # raw p-value
  interaction.cth[,6][[i]] <- round(summary(cth.model_list[[i]])[[10]][7,5], 5)
  if (sum(is.na(interaction.cth[,6]) < length(interaction.cth[,6]))) {
    names(interaction.cth) <- c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr', 'fdr2')
    interaction.cth[,7] <- round(p.adjust(interaction.cth[,6]), 5)
    interaction.cth[,8] <- round(mutoss::adaptiveBH(interaction.cth[,6], 0.05, silent = T)$adjPValues, 5) }
}

print(interaction.cth)

print(subset(interaction.cth, interaction.cth$`p-value` < 0.05))

sum(interaction.cth$`p-value` < 0.05)

(passed.p.int <- subset(interaction.cth$ROI, interaction.cth$`p-value` < 0.05))

print(subset(interaction.cth, interaction.cth$fdr < 0.05))

sum(interaction.cth$fdr < 0.05)

(cth.rm <- intersect(passed.p.int, passed.fdr.cth))       # IMPORTANT: is there any of the regions with a significant interaction (p < 0.05) also show a main effect?

(passed.fdr.cth <- setdiff(passed.fdr.cth, cth.rm))       # if so, remove the main effect from the original main effect results list 

length(passed.fdr.cth)

# Update raw unstandardized coefficients and change with 27 regions instead of 28 (if there was a overlap between main/interaction, if not the results is the same as above)

significant.models.list <- lapply(passed.fdr.cth, function(model_name) {
  cth.model_list[[model_name]]
})

dv <- c()
raw.beta <- c()
intercept <- c() 

for (i in 1:length(passed.fdr.cth)) {
  dv[i] <- passed.fdr.cth[i]
  raw.beta[i] <- summary(significant.models.list[[i]])[[10]][3,1]
  intercept[i] <- summary(significant.models.list[[i]])[[10]][1,1]
}

(ct.change <- cbind(raw.beta, intercept) %>% 
    as.data.frame() %>% 
    mutate(change = (raw.beta / intercept) * 100))

mean(ct.change$raw.beta)

mean(ct.change$change)

#~~~~~~~~~~~~~#
# Save models #
#~~~~~~~~~~~~~#

for (i in passed.fdr.cth) {
  model <- cth.model_list[[i]]
  save(model, file = paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/models/main_roi_cth_', i, '.Rdata'))
}

posthoc1 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc1$ROI, '_cth')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc1$ROI, '_cth')) %>% 
           select(label2) %>% 
           pull()) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results.csv'), row.names = F)

posthoc1 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc1$ROI, '_cth')) %>% 
           select(label2) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc1$ROI, '_cth')) %>% 
           select(label2) %>% 
           pull()) %>% 
  filter(`p-value` < 0.05) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results_p.csv'), row.names = F)

posthoc1 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc1$ROI, '_cth')) %>% 
           select(label2) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc1$ROI, '_cth')) %>% 
           select(label2) %>% 
           pull()) %>% 
  filter(ROI %in% passed.fdr.cth) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results_fdr.csv'), row.names = F)

#####################################
#####################################
#5# Load AREA data for main analysis
#####################################
#####################################

area.df <- brain %>% 
  select(., !ends_with(c('_myelin', '_cth', '_ratio'))) %>% 
  as.data.table()

# Select only the predictors you want to analyze for VIF
predictors <- area.df[, c('child_age', 'child_sex', 'premature', 'bfq_duration', 'education2', 'smri_vol_scs_intracranialv')]

# Calculate VIFs
(vif_results <- calculate_vif(predictors))

(featurenames <- names(area.df %>% select(starts_with('mrisdp'))))

area.model_list <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + premature + education2 + child_sex + smri_vol_scs_intracranialv + (1 | src_subject_id)")  
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = area.df)
  area.model_list[[feature]] <- model
}

posthoc2 <- as.data.frame(matrix(NA, nrow = length(featurenames), ncol = 8))

for (i in 1:length(area.model_list)) {
  posthoc2[,1][[i]] <- names(area.model_list[i])
  posthoc2[,2][[i]] <- length(!is.na(resid(area.model_list[[i]])))
  ci.tmp <- effectsize::standardize_parameters(area.model_list[[i]], method = 'basic')
  posthoc2[,3][[i]] <- ci.tmp[[2]][3]        # 3 for main, 8 for interaction --> why 8 and not 7? because in this model we have an extra variable (smri_vol_scs_intracranialv)
  posthoc2[,4][[i]] <- ci.tmp[[4]][3]
  posthoc2[,5][[i]] <- ci.tmp[[5]][3]
  # raw p-value
  posthoc2[,6][[i]] <- round(summary(area.model_list[[i]])[[10]][3,5], 5)
  if (sum(is.na(posthoc2[,6]) < length(posthoc2[,6]))) {
    names(posthoc2) <- c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr', 'fdr2')
    posthoc2[,7] <- round(p.adjust(posthoc2[,6], method = 'fdr'), 5)
    posthoc2[,8] <- round(mutoss::adaptiveBH(posthoc2[,6], 0.05, silent = T)$adjPValues, 5) }
    # posthoc2[,8] <- round(cp4p::adjust.p(posthoc2[,6], pi0.method = 'bky', alpha = 0.05)$adjp[,2], 5) }
}

print(posthoc2)

print(subset(posthoc2, posthoc2$`p-value` < 0.05))

sum(posthoc2$`p-value` < 0.05)

print(subset(posthoc2, posthoc2$fdr2 < 0.05))

sum(posthoc2$fdr < 0.05)

(passed.fdr.area <- subset(posthoc2$ROI, posthoc2$fdr < 0.05))

significant.models.list <- lapply(passed.fdr.area, function(model_name) {
  area.model_list[[model_name]]
})

dv <- c()
raw.beta <- c()
intercept <- c()  

for (i in 1:length(passed.fdr.area)) {
  dv[i] <- passed.fdr.area[i]
  raw.beta[i] <- summary(significant.models.list[[i]])[[10]][3,1]
  intercept[i] <- summary(significant.models.list[[i]])[[10]][1,1]
}

(area.change <- cbind(dv, raw.beta, intercept) %>% 
  as.data.frame() %>% 
  mutate(raw.beta = as.numeric(raw.beta),
         intercept = as.numeric(intercept),
          perc_change = (raw.beta / intercept) * 100))

mean(area.change$raw.beta)

mean(area.change$perc_change)

area.change$dv[which.max(area.change$perc_change)]

# Look interactions

interaction.area <- as.data.frame(matrix(NA, nrow = length(featurenames), ncol = 7))

for (i in 1:length(area.model_list)) {
  interaction.area[,1][[i]] <- names(area.model_list[i])
  interaction.area[,2][[i]] <- length(!is.na(resid(area.model_list[[i]])))
  ci.tmp <- effectsize::standardize_parameters(area.model_list[[i]], method = 'basic')
  interaction.area[,3][[i]] <- ci.tmp[[2]][8]    # 3 main effect, 8 interaction
  interaction.area[,4][[i]] <- ci.tmp[[4]][8]
  interaction.area[,5][[i]] <- ci.tmp[[5]][8]
  # raw p-value
  interaction.area[,6][[i]] <- round(summary(area.model_list[[i]])[[10]][8,5], 5)
  if (sum(is.na(interaction.area[,6]) < length(interaction.area[,6]))) {
    names(interaction.area) <- c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
    interaction.area[,7] <- round(p.adjust(interaction.area[,6]), 5) }
}

print(interaction.area)

print(subset(interaction.area, interaction.area$`p-value` < 0.05))

sum(interaction.area$`p-value` < 0.05)

(passed.p.int <- subset(interaction.area$ROI, interaction.area$`p-value` < 0.05))

print(subset(interaction.area, interaction.area$fdr < 0.05))

intersect(passed.p.int, passed.fdr.area)

(area.rm <- intersect(passed.p.int, passed.fdr.area))

(passed.fdr.area <- setdiff(passed.fdr.area, area.rm))

length(passed.fdr.area)    # 51 instead of 52 because 1 region (mrisdp_146_area) showed both main and interaction effects and we favor the higher-order term (interaction)

# recalculate

dv <- c()
raw.beta <- c()
intercept <- c()  

for (i in 1:length(passed.fdr.area)) {
  dv[i] <- passed.fdr.area[i]
  raw.beta[i] <- summary(significant.models.list[[i]])[[10]][3,1]
  intercept[i] <- summary(significant.models.list[[i]])[[10]][1,1]
}

(area.change <- cbind(dv, raw.beta, intercept) %>% 
    as.data.frame() %>% 
    mutate(raw.beta = as.numeric(raw.beta),
           intercept = as.numeric(intercept),
           perc_change = (raw.beta / intercept) * 100))

mean(area.change$raw.beta)

mean(area.change$perc_change)

#~~~~~~~~~~~~~#
# Save models #
#~~~~~~~~~~~~~#

for (i in passed.fdr.area) {
  model <- area.model_list[[i]]
  save(model, file = paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/models/main_roi_area_', i, '.Rdata'))
}

posthoc2 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc2$ROI, '_area')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc2$ROI, '_area')) %>% 
           select(label2) %>% 
           pull()) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/AREA_results.csv'), row.names = F)

posthoc2 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc2$ROI, '_area')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc2$ROI, '_area')) %>% 
           select(label2) %>% 
           pull()) %>% 
  filter(`p-value` < 0.05) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/AREA_results_p.csv'), row.names = F)

posthoc2 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc2$ROI, '_area')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc2$ROI, '_area')) %>% 
           select(label2) %>% 
           pull()) %>% 
  filter(ROI %in% passed.fdr.area) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/AREA_results_fdr.csv'), row.names = F)

# Overlap between CT and SA?

passed.fdr.cth
passed.fdr.area
length(common <- intersect(str_remove(passed.fdr.cth, '_cth'), str_remove(passed.fdr.area, '_area'))) # how many overlapping CT and SA regions
common # overlaping CT and SA regions with a main effect for breastfeeding


