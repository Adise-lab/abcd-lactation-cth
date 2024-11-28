#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for  analyzing ABCD data release 5.1 for the paper entitled                  #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Aug 19th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, summarytools)

options(scipen = 999)

#2# Set paths

datadir <- '/Users/jongonzalez/Desktop/abcd.nosync/'

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

################################################################
################################################################
#6# Conduct follow-up analysis on CM on significant CT regions #
################################################################
################################################################

brain <- 
  # fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_original.csv'), stringsAsFactors = T) %>%               # ORIGINAL data (only combatt'd)
  fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_no_outliers.csv'), stringsAsFactors = T) %>%              # no outliers (NAs) and those with NAs at baselien but not follow-up set to NAs at both timepoints
  # fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_NO_COMBAT.csv'), stringsAsFactors = T) %>%             # regular data (no combatt'd)
  mutate_at(vars(all_of(covs2scale)), ~scale(., scale = F)) %>%  # mean-center continous covariates
  mutate_at(vars(all_of(covs2scale)), ~scale(., scale = F)) %>%  # mean-center continous covariates
  mutate(handedness = factor(handedness),
         child_sex = ifelse(child_sex == 'Male', 1, -1),   # effects code sex
         premature = ifelse(premature == 'Yes', 1, -1),    # effects code premature
         interaction = child_age * bfq_duration)

#3# Set common covariates

(covs2scale <- c('child_age', 'bfq_duration', 'education2', 'smri_vol_scs_intracranialv'))

# generate a list of features to be tested --> from the CT regions that passed fdr correction, remove _cth suffix and replace by _t1t2_ratio

(passed.fdr.cth.myel <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results_fdr.csv')) %>% 
   mutate(ROI = str_replace(ROI, '_cth', '_t1t2_ratio')) %>% 
   select(ROI) %>% 
   pull)

length(passed.fdr.cth.myel)    # confirm is the same number as the main analysis on CT

# now extract from the brain dataframe the regions that do not match cth or area (= select the ones with _t1t2_ratio ending)

myel.df <- brain %>% 
  select(., !ends_with(c('_cth', '_area'))) %>% 
  as.data.table()

contrasts(myel.df$handedness) <- contr.sum(levels(myel.df$handedness))

(featurenames2 <- names(myel.df %>% select(all_of(passed.fdr.cth.myel))))

model_list.myel <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames2) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + premature + education2 + child_sex + (1 | src_subject_id)")  
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = myel.df)
  model_list.myel[[feature]] <- model
}

posthoc3 <- as.data.frame(matrix(NA, nrow = length(featurenames2), ncol = 10))

names(posthoc3) <-  c('ROI', 'n', 'beta.main', 'ci.min.main', 'ci.max.main', 'p-value.main', 'beta.int', 'ci.min.int', 'ci.max.int', 'p-value.int')

for (i in 1:length(model_list.myel)) {
  # name ROIs
  posthoc3[,1][[i]] <- names(model_list.myel[i])
  # N models 
  posthoc3[,2][[i]] <- length(!is.na(resid(model_list.myel[[i]])))
  # CI tmp
  ci.tmp <- effectsize::standardize_parameters(model_list.myel[[i]], method = 'basic')
  # beta and CI main
  posthoc3[,3][[i]] <-ci.tmp[[2]][3]
  posthoc3[,4][[i]] <- ci.tmp[[4]][3]
  posthoc3[,5][[i]] <- ci.tmp[[5]][3]
  # raw p-value main
  posthoc3[,6][[i]] <- summary(model_list.myel[[i]])[[10]][3,5]
  # beta and CI interaction
  posthoc3[,7][[i]] <- effectsize::standardize_parameters(model_list.myel[[i]], method = 'basic')[[2]][7]
  posthoc3[,8][[i]] <- ci.tmp[[4]][7]
  posthoc3[,9][[i]] <- ci.tmp[[5]][7]
  # raw p-value interaction
  posthoc3[,10][[i]] <- summary(model_list.myel[[i]])[[10]][7,5]
}

print(posthoc3)

print(subset(posthoc3[1:6], posthoc3$`p-value.main` < 0.05))

sum(posthoc3$`p-value.main` < 0.05)

(passed.p.myel.cth.main <- subset(posthoc3$ROI, posthoc3$`p-value.main` < 0.05))

posthoc3 %>% 
  filter(`p-value.int` < 0.05) %>% 
  select(ROI, beta.int, ci.min.int, ci.max.int, `p-value.int`)

sum(posthoc3$`p-value.int` < 0.05)

(passed.p.myel.cth.int <- subset(posthoc3$ROI, posthoc3$`p-value.int` < 0.05))

# Review the plots to understand how to interpret
for (i in passed.p.myel.cth.int) {
  plot_data <- ggeffects::ggpredict(model_list.myel[[i]], terms = c('child_age', 'bfq_duration')) %>% plot()
  plot(plot_data)
  cat('Press [Enter] to continue to the next plot...')
  invisible(readline())
}

#~~~~~~~~~~~~~~#
# Saved models #
#~~~~~~~~~~~~~~#

passed.p.myel.cth <- union(passed.p.myel.cth.main, passed.p.myel.cth.int)

for (i in passed.p.myel.cth) {
  model <- model_list.myel[[i]]
  save(model, file = paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/models/myel_ct_roi_', i, '.Rdata'))
}

posthoc3 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc3$ROI, '_t1t2_ratio')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc3$ROI, '_t1t2_ratio')) %>% 
           select(label2) %>% 
           pull()) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/MYEL_CT_results.csv'), row.names = F)

################################################################
################################################################
#7# Conduct follow-up analysis on CM on SA significant regions #
################################################################
################################################################

(passed.fdr.area.myel <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Area_results_fdr.csv')) %>% 
   mutate(ROI = str_replace(ROI, '_area', '_t1t2_ratio')) %>% 
   select(ROI) %>% 
   pull())

myel.df <- brain %>% 
  select(., !ends_with(c('_cth', '_area'))) %>% 
  as.data.table()

(featurenames2 <- names(myel.df %>% select(all_of(passed.fdr.area.myel))))

model_list.myel <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames2) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + premature + education2 + child_sex + (1 | src_subject_id)")  
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = myel.df)
  model_list.myel[[feature]] <- model
}

posthoc4 <- as.data.frame(matrix(NA, nrow = length(featurenames2), ncol = 10))

names(posthoc4) <-  c('ROI', 'n', 'beta.main', 'ci.min.main', 'ci.max.main', 'p-value.main', 'beta.int', 'ci.min.int', 'ci.max.int', 'p-value.int')

for (i in 1:length(model_list.myel)) {
  # name ROIs
  posthoc4[,1][[i]] <- names(model_list.myel[i])
  # N models 
  posthoc4[,2][[i]] <- length(!is.na(resid(model_list.myel[[i]])))
  # CI tmp
  ci.tmp <- effectsize::standardize_parameters(model_list.myel[[i]], method = 'basic')
  # beta and CI main
  posthoc4[,3][[i]] <- effectsize::standardize_parameters(model_list.myel[[i]], method = 'basic')[[2]][3]
  posthoc4[,4][[i]] <- ci.tmp[[4]][3]
  posthoc4[,5][[i]] <- ci.tmp[[5]][3]
  # raw p-value main
  posthoc4[,6][[i]] <- summary(model_list.myel[[i]])[[10]][3,5]
  # beta and CI interaction
  posthoc4[,7][[i]] <- effectsize::standardize_parameters(model_list.myel[[i]], method = 'basic')[[2]][7]
  posthoc4[,8][[i]] <- ci.tmp[[4]][7]
  posthoc4[,9][[i]] <- ci.tmp[[5]][7]
  # raw p-value interaction
  posthoc4[,10][[i]] <- summary(model_list.myel[[i]])[[10]][7,5]
}

print(posthoc4)

print(subset(posthoc4[1:6], posthoc4$`p-value.main` < 0.05))

sum(posthoc4$`p-value.main` < 0.05)

(passed.p.myel.area.main <- subset(posthoc4$ROI, posthoc4$`p-value.main` < 0.05))

posthoc4 %>% 
  filter(`p-value.int` < 0.05) %>% 
  select(ROI, beta.int, ci.min.int, ci.max.int, `p-value.int`)

sum(posthoc4$`p-value.int` < 0.05)

(passed.p.myel.area.int <- subset(posthoc4$ROI, posthoc4$`p-value.int` < 0.05))

# Review the plots to understand how to interpret
for (i in passed.p.myel.area.int) {
  plot_data <- ggeffects::ggpredict(model_list.myel[[i]], terms = c('child_age', 'bfq_duration')) %>% plot()
  plot(plot_data)
  cat('Press [Enter] to continue to the next plot...')
  invisible(readline())
}

#~~~~~~~~~~~~~~#
# Saved models #
#~~~~~~~~~~~~~~#

passed.p.myel.area <- union(passed.p.myel.area.main, passed.p.myel.area.int)

for (i in passed.p.myel.area) {
  model <- model_list.myel[[i]]
  save(model, file = paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/models/myel_area_roi_', i, '.Rdata'))
}

posthoc4 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc4$ROI, '_t1t2_ratio')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc4$ROI, '_t1t2_ratio')) %>% 
           select(label2) %>% 
           pull()) %>% 
  select(ROI, label, label2, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/MYEL_AREA_results.csv'), row.names = F)

# Overlap CT, SA and CM?

intersect(passed.p.myel.cth.main, passed.p.myel.cth.int)      # CAUTION: main and interaction effect in this region warrant further examination

intersect(passed.p.myel.area.main, passed.p.myel.area.int)    # CAUTION: main and interaction effect in this region warrant further examination

# Since the effects were not truly due to breastfeeding alone, will remove from passed.p.myel.area.main but keep at passed.p.myel.area.int

(passed.p.myel.cth.main <- setdiff(passed.p.myel.cth.main, intersect(passed.p.myel.cth.main, passed.p.myel.cth.int)))

(passed.p.myel.area.main <- setdiff(passed.p.myel.area.main, intersect(passed.p.myel.area.main, passed.p.myel.area.int)))

# Confirm all regions (CT, SA, myelin) are associated with BF at baseline

(passed.fdr.cth <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results_fdr.csv')) %>% 
  pull(ROI))

(passed.fdr.area <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Area_results_fdr.csv')) %>% 
    pull(ROI))

length(featurenames2 <- union(union(union(passed.p.myel.cth.main, passed.p.myel.area.main), passed.fdr.area), passed.fdr.cth))

baseline_list <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames2) {
  # Check if the feature name ends with '_area'
  if (grepl("_area$", feature)) {
    # Use the formula with 'total_volume' for features ending with '_area'
    formula_str <- paste0(feature, " ~ child_age + bfq_duration + premature + education2 + child_sex + smri_vol_scs_intracranialv")
  } else {
    # Use the standard formula for other features
    formula_str <- paste0(feature, " ~ child_age + bfq_duration + premature + education2 + child_sex")
  }
  
  # Fit the model using lm and store in baseline_list
  model <- lm(as.formula(formula_str), subset(brain, brain$eventname == 'baseline_year_1_arm_1'))
  baseline_list[[feature]] <- model
}

posthoc4.5 <- as.data.frame(matrix(NA, nrow = length(featurenames2), ncol = 6))

names(posthoc4.5) <-  c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value')

for (i in 1:length(baseline_list)) {
  # name ROIs
  posthoc4.5[,1][[i]] <- names(baseline_list[i])
  # N models 
  posthoc4.5[,2][[i]] <- length(!is.na(resid(baseline_list[[i]])))
  # CI tmp
  ci.tmp <- effectsize::standardize_parameters(baseline_list[[i]], method = 'basic')
  # beta and CI main
  posthoc4.5[,3][[i]] <- effectsize::standardize_parameters(baseline_list[[i]], method = 'basic')[[2]][3]
  posthoc4.5[,4][[i]] <- ci.tmp[[4]][3]
  posthoc4.5[,5][[i]] <- ci.tmp[[5]][3]
  # raw p-value main
  posthoc4.5[,6][[i]] <- summary(baseline_list[[i]])[[4]][3,4]
}

print(posthoc4.5)

length(posthoc4.5$ROI)

subset(posthoc4.5, posthoc4.5$`p-value` < 0.05)

sum(posthoc4.5$`p-value` < 0.05)

length(posthoc4.5 %>% 
  filter(`p-value` < 0.05 &
           str_ends(ROI, 'cth')) %>% 
  pull(ROI) -> new.cth)

length(posthoc4.5 %>% 
  filter(`p-value` < 0.05 &
           str_ends(ROI, 'area')) %>% 
  pull(ROI) -> new.area)

length(posthoc4.5 %>% 
  filter(`p-value` < 0.05 &
           str_ends(ROI, 'ratio')) %>% 
  pull(ROI) -> new.myelin)

#################################################
#################################################
#8# Conduct supplementary analysis on cognition #
#################################################
#################################################

(passed.p.main <- union(union(new.cth, new.area), new.myelin))

(passed.p.myel_interaction <- union(passed.p.myel.cth.int, passed.p.myel.area.int))

(covs2scale2 <- c(covs2scale, 'vision'))

cog <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv'), stringsAsFactors = T) %>% 
  inner_join(., brain %>% 
               select(., src_subject_id, eventname, matches(c(passed.p.main, passed.p.myel_interaction))), 
             by = c('src_subject_id', 'eventname'))

cog.df <- cog %>% 
  # filter for baseline as Shana asked
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  mutate_at(vars(all_of(covs2scale2)), ~as.numeric(scale(., scale = F))) %>% 
  mutate(handedness = as.factor(handedness),
         child_sex = ifelse(child_sex == 'Male', 1, -1),
         premature = ifelse(premature == 'Yes', 1, -1)) %>% 
  as.data.table()

# table for report
vtable::st(fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv'), stringsAsFactors = T) %>% 
             mutate_at(vars(alcohol, tobacco, marihuana, opioids, opiates, cocaine), ~factor(ifelse(. == 1, 'Yes', ifelse(. == 0, 'No', 'Do not know')))) %>% 
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
           vars = c('child_age', 'child_sex', 'pds', 'handedness', 'vision', 'race', 'ethnicity', 'mode_of_delivery', 'birth_complications', 'pregnancy_complications', 
                    'premature', 'premature_wk', 'education3', 'annual_income', 'alcohol', 'tobacco', 'marihuana', 'opioids', 'opiates', 'cocaine', 'bfq_duration'), 
           group = 'eventname')

summary(model.cog <- lm(nih_fluid ~ bfq_duration + child_age + premature + education2 + child_sex + vision, cog.df))

# summary(model.cog <- lmerTest::lmer(nih_fluid ~ bfq_duration * child_age + premature + education2 + child_sex + vision + (1 | src_subject_id), cog.df))

effectsize::standardize_parameters(model.cog, method = 'basic')

##################################################
##################################################
# Now test what regions are related to cognition #
##################################################
##################################################

length(baseline.brain <- names(cog.df %>% select(all_of(passed.p.main))))

baseline.brain

model_list <- list()

# Define main predictors
(main_predictors <- names(cog.df %>% select(., all_of(baseline.brain))))     # passed.fdr.cth, passed.fdr.area, passed.fdr.cth.myel, passed.fdr.area.myel

# only run in baseline because bf effects on brain and cognition are stable
cog.df2 <- cog %>% 
  filter(eventname == 'baseline_year_1_arm_1') %>%
  mutate_at(vars(all_of(main_predictors)), ~scale(., scale = F))

for (main_pred in main_predictors) {
  # Construct the formula string for the current feature and main predictor
  formula_str <- paste0("nih_fluid ~", main_pred, "+ child_age + premature + education2 + child_sex + vision")
  # Fit the model using lmerTest::lmer and store in model_list
  model_list[[main_pred]] <- lm(as.formula(formula_str), data = cog.df2)
}

posthoc6 <- as.data.frame(matrix(NA, nrow = length(main_predictors), ncol = 7))

for (i in 1:length(model_list)) {
  posthoc6[,1][[i]] <- names(model_list[i])
  posthoc6[,2][[i]] <- length(!is.na(resid(model_list[[i]])))
  # beta estimate
  posthoc6[,3][[i]] <- effectsize::standardize_parameters(model_list[[i]], method = 'basic')[[2]][2]
  # confidence intervals
  ci.tmp <- effectsize::standardize_parameters(model_list[[i]], method = 'basic')
  posthoc6[,4][[i]] <- ci.tmp[[4]][2]
  posthoc6[,5][[i]] <- ci.tmp[[5]][2]
  # raw p-value
  posthoc6[,6][[i]] <- round(summary(model_list[[i]])[[4]][2,4], 5)
  # posthoc6[,5][[i]] <- round(summary(model_list[[i]])[[10]][3,5], 5)
  if (sum(is.na(posthoc6[,6]) < length(posthoc6[,6]))) {
    names(posthoc6) <- c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
    posthoc6[,7] <- round(p.adjust(posthoc6[,6], method = 'fdr'), 3) }
}

print(subset(posthoc6, posthoc6$`p-value` < 0.05))

print(subset(posthoc6$ROI, posthoc6$`p-value` < 0.05))

sum(posthoc6$`p-value` < 0.05)

(passed.p2 <- subset(posthoc6$ROI, posthoc6$`p-value` < 0.05))

posthoc6 %>% 
  filter(str_ends(ROI, '_cth')) %>% 
  mutate(ROI = str_remove(ROI, '_cth')) %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc6$ROI, '_cth')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc6$ROI, '_cth')) %>% 
           select(label2) %>% 
           pull()) %>% 
  mutate(ROI = paste0(ROI, '_cth')) %>% 
  select(ROI, label, label2, everything()) -> cog.cth

cog.cth %>% 
  filter(`p-value` < 0.05) %>% 
  nrow()

posthoc6 %>% 
  filter(str_ends(ROI, '_area')) %>% 
  mutate(ROI = str_remove(ROI, '_area')) %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc6$ROI, '_area')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc6$ROI, '_area')) %>% 
           select(label2) %>% 
           pull()) %>% 
  mutate(ROI = paste0(ROI, '_area')) %>% 
  select(ROI, label, label2, everything()) -> cog.area

cog.area %>% 
  filter(`p-value` < 0.05) %>% 
  nrow()

posthoc6 %>% 
  filter(str_ends(ROI, '_t1t2_ratio')) %>% 
  mutate(ROI = str_remove(ROI, '_t1t2_ratio')) %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc6$ROI, '_t1t2_ratio')) %>% 
           select(ggseg) %>% 
           pull()) %>% 
  mutate(label2 = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% str_remove(posthoc6$ROI, '_t1t2_ratio')) %>% 
           select(label2) %>% 
           pull()) %>% 
  mutate(ROI = paste0(ROI, '_t1t2_ratio')) %>% 
  select(ROI, label, everything()) -> cog.myelin

cog.myelin %>% 
  filter(`p-value` < 0.05) %>% 
  nrow()

brain.cog <- rbind(cog.cth, cog.area, cog.myelin) %>% 
  mutate(measure = ifelse(str_ends(ROI, '_cth'), 'Thickness',
                          ifelse(str_ends(ROI, '_area'), 'Area', 'Myelin'))) %>% 
  select(ROI, label, measure, everything())

write.csv(brain.cog, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/brain_cog_results.csv'), row.names = F)

# Make sure no CM interaction-regions predict fluid cognition at y2

# length(year2.brain <- names(cog.df %>% select(all_of(passed.p.myel_interaction))))
# 
# year2.brain
# 
# model_list <- list()
# 
# # Define main predictors
# (main_predictors <- names(cog.df %>% select(., all_of(year2.brain))))     # passed.fdr.cth, passed.fdr.area, passed.fdr.cth.myel, passed.fdr.area.myel
# 
# # only run in baseline because bf effects on brain and cognition are stable
# cog.df3 <- cog.df %>%
#   filter(eventname != 'baseline_year_1_arm_1') %>%
#   mutate_at(vars(all_of(main_predictors)), ~scale(., scale = F))
# 
# for (main_pred in main_predictors) {
#   # Construct the formula string for the current feature and main predictor
#   formula_str <- paste0("nih_fluid ~", main_pred, "+ child_age + premature + education2 + child_sex + vision")
#   # Fit the model using lmerTest::lmer and store in model_list
#   model_list[[main_pred]] <- lm(as.formula(formula_str), data = cog.df3)
# }
# 
# posthoc7 <- as.data.frame(matrix(NA, nrow = length(main_predictors), ncol = 7))
# 
# for (i in 1:length(model_list)) {
#   posthoc7[,1][[i]] <- names(model_list[i])
#   posthoc7[,2][[i]] <- length(!is.na(resid(model_list[[i]])))
#   # beta estimate
#   posthoc7[,3][[i]] <- effectsize::standardize_parameters(model_list[[i]], method = 'basic')[[2]][2]
#   # confidence intervals
#   ci.tmp <- effectsize::standardize_parameters(model_list[[i]], method = 'basic')
#   posthoc7[,4][[i]] <- ci.tmp[[4]][2]
#   posthoc7[,5][[i]] <- ci.tmp[[5]][2]
#   # raw p-value
#   posthoc7[,6][[i]] <- round(summary(model_list[[i]])[[4]][2,4], 5)
#   # posthoc6[,5][[i]] <- round(summary(model_list[[i]])[[10]][3,5], 5)
#   if (sum(is.na(posthoc7[,6]) < length(posthoc7[,6]))) {
#     names(posthoc7) <- c('ROI', 'n', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
#     posthoc7[,7] <- round(p.adjust(posthoc7[,6], method = 'fdr'), 3) }
# }
# 
# print(posthoc7)
# 
# print(subset(posthoc7, posthoc7$`p-value` < 0.05))
# 
# sum(posthoc7$`p-value` < 0.05)

# Continue

covs <- c('child_age', 'vision', 'bfq_duration', 'education2')

(neg <- subset(posthoc6$ROI, posthoc6$`p-value` < 0.05 & posthoc6$beta < 0))

(cth.composite <- setdiff(subset(cog.cth$ROI, cog.cth$`p-value` < 0.05), neg))

(area.composite <- setdiff(subset(cog.area$ROI, cog.area$`p-value` < 0.05), neg))

(myelin.composite <- subset(cog.myelin$ROI, cog.myelin$`p-value` < 0.05))

med.y0 <- cog %>% 
  filter(eventname == 'baseline_year_1_arm_1') %>%
  mutate_at(vars(all_of(covs)), ~as.numeric(scale(., scale = F))) %>% 
  # it is important to set na.rm = T because if there's any NA in the brain variables it will return a mean of NA (it is not the case but...)
  mutate(brain.cth = scale(rowMeans(select(., all_of(cth.composite)), na.rm = T), scale = F),
         # brain.myel = scale(rowMeans(select(., all_of(myelin.composite)), na.rm = T), scale = F),
         brain.area = scale(rowMeans(select(., all_of(area.composite)), na.rm = T), scale = F),
         cognition = nih_fluid,
         child_sex = ifelse(child_sex == 'Male', 1, -1),
         premature = ifelse(premature == 'Yes', 1, -1))

# PROCESS 

source(paste0(drivedir, 'Lab_member_files/Jonatan_files/processv43/PROCESS v4.3 for R/process.R'))

process(as.data.frame(med.y0), 
        y = 'nih_fluid', 
        x = 'bfq_duration', 
        m = c('brain.cth',
              # 'brain.myel',    
              'brain.area'), 
        model = 4,        # 4 parallel mediation, 6 serial mediation
        normal = 1,       # sobel test
        effsize = 1,      # effect sizes
        total = 1,        # total effect (c-path)
        stand = 1,        # standardized coefficients
        modelbt = 1,      # bootstrap all effects (not only indirect)
        contrast = 1,     # pairwise comparison between indirect effects
        boot = 10000,
        cov = c('child_age', 
                'premature',
                'education2', 
                'child_sex', 
                'vision'), 
        seed = 1904)
