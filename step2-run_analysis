#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for  analyzing ABCD data release 5.1 for the paper entitled                  #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Aug 14th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, summarytools)

#2# Set paths

datadir <- '/Users/jongonzalez/Desktop/abcd.nosync/'

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

#3# Load MRI data (CT) for main analysis

(covs2scale <- c('child_age', 'bfq_duration', 'premature_wk', 'education_max', 'handedness_c'))

ct.df <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data.csv')) %>% 
  mutate_at(vars(all_of(covs2scale)), ~scale(., scale = F)) %>% 
  select(., !ends_with('_myelin'))

(featurenames <- names(ct.df %>% select(starts_with('mrisdp'))))

model_list <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + premature_wk + education_max + handedness_c + child_sex + (1 | src_subject_id)")  
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = ct.df)
  model_list[[feature]] <- model
}

posthoc1 <- as.data.frame(matrix(NA, nrow = length(featurenames), ncol = 6))

for (i in 1:length(model_list)) {
  posthoc1[,1][[i]] <- names(model_list[i])
  # beta estimate
  posthoc1[,2][[i]] <- summary(model_list[[i]])[[10]][3,1]     # 3 for main, 11 for interaction
  # confidence intervals min [1] and max [2]
  ci.tmp <- as.numeric(confint(model_list[[i]], 'bfq_duration'))
  posthoc1[,3][[i]] <- ci.tmp[1]
  posthoc1[,4][[i]] <- ci.tmp[2]
  # raw p-value
  posthoc1[,5][[i]] <- summary(model_list[[i]])[[10]][3,5]
  if (sum(is.na(posthoc1[,5]) < length(posthoc1[,5]))) {
    names(posthoc1) <- c('ROI', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
    posthoc1[,6] <- round(p.adjust(posthoc1[,5], method = 'fdr'), 5) }
}

print(posthoc1)

print(subset(posthoc1, posthoc1$`p-value` < 0.05))

sum(posthoc1$`p-value` < 0.05)

print(subset(posthoc1, posthoc1$fdr < 0.05))

(passed.fdr <- subset(posthoc1$ROI, posthoc1$fdr < 0.05))

sum(posthoc1$fdr < 0.05)

posthoc1 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% posthoc1$ROI) %>% 
           select(ggseg) %>% pull()) %>% 
  select(ROI, label, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/CT_results.csv'), row.names = F)

posthoc1 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% posthoc1$ROI) %>% 
           select(ggseg) %>% pull()) %>% 
  filter(`p-value` < 0.05) %>% 
  select(ROI, label, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/CT_results_p.csv'), row.names = F)

posthoc1 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% posthoc1$ROI) %>% 
           select(ggseg) %>% pull()) %>% 
  filter(fdr < 0.05) %>% 
  select(ROI, label, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/CT_results_fdr.csv'), row.names = F)

#4# Conduct follow-up analysis on Intracortical Myelin Content

myel.df <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data.csv')) %>% 
  mutate_at(vars(all_of(covs2scale)), ~scale(., scale = F)) %>% 
  rename_with(~ ifelse(grepl('^mrisdp_\\d+$', .) & !grepl('_myelin$', .), paste0(., '_ct'),
                       ifelse(grepl('_myelin$', .), sub('_myelin$', '', .), .))) %>% 
  select(-ends_with('_ct'))

(featurenames2 <- names(myel.df %>% select(all_of(passed.fdr))))

model_list <- list()        # Initialize an empty list to store models

# Iterate over each feature name
for (feature in featurenames2) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + premature_wk + education_max + handedness_c + child_sex + (1 | src_subject_id)")  
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = myel.df)
  model_list[[feature]] <- model
}

posthoc2 <- as.data.frame(matrix(NA, nrow = length(featurenames2), ncol = 6))

for (i in 1:length(model_list)) {
  posthoc2[,1][[i]] <- names(model_list[i])
  # beta estimate
  posthoc2[,2][[i]] <- summary(model_list[[i]])[[10]][3,1]     # 3 for main, 11 for interaction
  # confidence intervals min [1] and max [2]
  ci.tmp <- as.numeric(confint(model_list[[i]], 'bfq_duration'))
  posthoc2[,3][[i]] <- ci.tmp[1]
  posthoc2[,4][[i]] <- ci.tmp[2]
  # raw p-value
  posthoc2[,5][[i]] <- summary(model_list[[i]])[[10]][3,5]
  if (sum(is.na(posthoc2[,5]) < length(posthoc2[,5]))) {
    names(posthoc2) <- c('ROI', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
    posthoc2[,6] <- round(p.adjust(posthoc2[,5], method = 'fdr'), 5) }
}
  
print(posthoc2)

print(subset(posthoc2, posthoc2$`p-value` < 0.05))

sum(posthoc2$`p-value` < 0.05)

posthoc2 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% posthoc2$ROI) %>% 
           select(ggseg) %>% pull()) %>% 
  select(ROI, label, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/Myel_results.csv'), row.names = F)

posthoc2 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% posthoc2$ROI) %>% 
           select(ggseg) %>% pull()) %>% 
  filter(`p-value` < 0.05) %>% 
  select(ROI, label, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/Myel_results_p.csv'), row.names = F)

#5# Conduct supplementary analysis on cognition

(covs2scale2 <- c(covs2scale, 'vision'))

cog.df <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv')) %>% 
  mutate_at(vars(all_of(covs2scale2)), ~scale(., scale = F)) %>% 
  select(., !ends_with('_myelin'))

(featurenames3 <- names(cog.df %>% select(starts_with(c('nih', 'pea', 'lmt')))))

cog.df$cognition <- cog.df %>% group_by(eventname) %>% mutate_at(vars(all_of(featurenames3)), ~scale(., scale = T)) %>% ungroup() %>% mutate(cognition = rowMeans(select(., all_of(featurenames3)))) %>% select(cognition) 

(featurenames3 <- c('cognition', featurenames3))

# Initialize an empty list to store models
model_list <- list()

# Iterate over each feature name
for (feature in featurenames3) {
  # Construct the formula string for the current feature
  formula_str <- paste0(feature, " ~ child_age * bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + premature_wk + education_max + handedness_c + child_sex + vision + (1 | src_subject_id)")
  # Fit the model using lmerTest::lmer and store in model_list
  model <- lmerTest::lmer(as.formula(formula_str), data = cog.df)
  model_list[[feature]] <- model
}

posthoc3 <- as.data.frame(matrix(NA, nrow = length(featurenames3), ncol = 6))

for (i in 1:length(model_list)) {
  posthoc3[,1][[i]] <- names(model_list[i])
  # beta estimate
  posthoc3[,2][[i]] <- round(summary(model_list[[i]])[[10]][3,1], 3)     # 3 for main, 13 for interaction
  # confidence intervals min [1] and max [2]
  ci.tmp <- confint(model_list[[i]], 'bfq_duration')
  posthoc3[,3][[i]] <- ci.tmp[1]
  posthoc3[,4][[i]] <- ci.tmp[2]
  # raw p-value
  posthoc3[,5][[i]] <- round(summary(model_list[[i]])[[10]][3,5], 5)
  if (sum(is.na(posthoc3[,5]) < length(posthoc3[,5]))) {
    names(posthoc3) <- c('ROI', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
    posthoc3[,6] <- round(p.adjust(posthoc3[,5], method = 'fdr'), 5) }
}

print(posthoc3)

print(subset(posthoc3, posthoc3$`p-value` < 0.05))

sum(posthoc3$`p-value` < 0.05)

# Now test what regions are related to cognition

model_list <- list()

# Define main predictors
(main_predictors <- names(cog.df %>% select(., all_of(passed.fdr))))

for (main_pred in main_predictors) {
  # Construct the formula string for the current feature and main predictor
  formula_str <- paste0("cognition ~ child_age * ", main_pred, "+ mode_of_delivery + birth_complications_n + pregnancy_complications + premature_wk + education_max + handedness_c + child_sex + vision + (1 | src_subject_id)")
  # Fit the model using lmerTest::lmer and store in model_list
  model_list[[main_pred]] <- lmerTest::lmer(as.formula(formula_str), data = cog.df)
}

posthoc4 <- as.data.frame(matrix(NA, nrow = length(main_predictors), ncol = 6))

for (i in 1:length(model_list)) {
  posthoc4[,1][[i]] <- names(model_list[i])
  # beta estimate
  posthoc4[,2][[i]] <- round(summary(model_list[[i]])[[10]][3,1], 3)     # 3 for main, 13 for interaction
  # confidence intervals min [1] and max [2]
  ci.tmp <- confint(model_list[[i]], names(model_list[i]))
  posthoc4[,3][[i]] <- ci.tmp[1]
  posthoc4[,4][[i]] <- ci.tmp[2]
  # raw p-value
  posthoc4[,5][[i]] <- round(summary(model_list[[i]])[[10]][3,5], 5)
  if (sum(is.na(posthoc4[,5]) < length(posthoc4[,5]))) {
    names(posthoc4) <- c('ROI', 'beta', 'ci.min', 'ci.max', 'p-value', 'fdr')
    posthoc4[,6] <- round(p.adjust(posthoc4[,5], method = 'fdr'), 3) }
}

print(posthoc4)

print(subset(posthoc4, posthoc4$`p-value` < 0.05))

sum(posthoc4$`p-value` < 0.05)

(passed.p <- subset(posthoc4$ROI, posthoc4$`p-value` < 0.05))

posthoc4 %>% 
  mutate(label = readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
           filter(roi %in% posthoc4$ROI) %>% 
           select(ggseg) %>% pull()) %>% 
  select(ROI, label, everything()) %>% 
  write.csv(., paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/CogCT_results.csv'), row.names = F)

(featurenames3 <- setdiff(featurenames3, 'cognition'))

covs <- c('mode_of_delivery', 'birth_complications_n', 'pregnancy_complications', 'premature_wk', 'handedness_c', 'education_max', 'child_sex', 'child_age', 'vision')

# lav.df <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv')) %>%
#   mutate(eventname = ifelse(eventname == 'baseline_year_1_arm_1', 'y0', 'y2')) %>%
#   select(., src_subject_id, eventname, child_age, bfq_duration, vision, all_of(c(passed.p, covs, featurenames3))) %>%
#   group_by(eventname) %>%
#   mutate_at(vars(c(all_of(c(featurenames3, passed.p)), child_age, bfq_duration, vision)), ~as.numeric(scale(., scale = T))) %>%
#   ungroup() %>%
#   mutate(brain = rowMeans(select(., all_of(passed.p))),
#          cognition = rowMeans(select(., all_of(featurenames3)))) %>%
#   select(-all_of(passed.p), -starts_with(c('nih', 'pea', 'lmt'))) %>%
#   pivot_wider(names_from = 'eventname', values_from = c(brain, cognition, child_age, vision))

extract_mediation_summary <- function (x) { 
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}

med.y0 <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv')) %>%
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(., src_subject_id, eventname, child_age, bfq_duration, vision, all_of(c(passed.p, covs)), starts_with(c('nih', 'pea', 'lmt'))) %>%
  mutate_at(vars(starts_with(c('nih', 'pea', 'lmt'))), ~as.numeric(scale(., scale = T))) %>%
  mutate_at(vars(c(pregnancy_complications, birth_complications_n)), ~as.factor(.)) %>% 
  mutate(mode_of_delivery = as.factor(ifelse(mode_of_delivery == 'Vaginal', 0, 1)),
         child_sex = as.factor(ifelse(child_sex == 'Male', 1, 0)),
         brain = rowMeans(select(., all_of(passed.p))),
         cognition = rowMeans(select(., starts_with(c('nih', 'pea', 'lmt'))))) %>%
  select(-all_of(passed.p), -starts_with(c('nih', 'pea', 'lmt'))) %>% 
  mutate(across(where(is.numeric), ~as.vector(scale(., scale = T))))

# baseline - baseline

covs_y0 <- c('mode_of_delivery', 'birth_complications_n', 'pregnancy_complications', 'premature_wk', 'handedness_c', 'education_max', 'child_sex', 'child_age', 'vision')

summary(mediator_modely0 <- lm(brain ~ bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + handedness_c + premature_wk + education_max + child_sex + child_age + vision, med.y0))
summary(outcome_modely0 <- lm(cognition ~ brain + bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + handedness_c + premature_wk + education_max + child_sex + child_age + vision, med.y0))

confint(outcome_modely0, 'bfq_duration')     # lactation on cognition
confint(mediator_modely0, 'bfq_duration')    # lactation on brain
confint(outcome_modely0, 'brain')            # brain on cognition

set.seed(4994)

mediation_resultsy0 <- mediation::mediate(mediator_modely0, outcome_modely0, treat = 'bfq_duration', mediator = 'brain', covariates = covs_y0, boot = T, sims = 5000)

(medy0 <- extract_mediation_summary(summary(mediation_resultsy0)))

write.csv(medy0, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/ABCD_medy0.csv'))

# follow up - follow up

med.y2 <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv')) %>%
  filter(eventname != 'baseline_year_1_arm_1') %>%
  select(., src_subject_id, eventname, child_age, bfq_duration, vision, all_of(c(passed.p, covs)), starts_with(c('nih', 'pea', 'lmt'))) %>%
  mutate_at(vars(starts_with(c('nih', 'pea', 'lmt'))), ~as.numeric(scale(., scale = T))) %>%
  mutate_at(vars(c(pregnancy_complications, birth_complications_n)), ~as.factor(.)) %>% 
  mutate(mode_of_delivery = as.factor(ifelse(mode_of_delivery == 'Vaginal', 0, 1)),
         child_sex = as.factor(ifelse(child_sex == 'Male', 1, 0)),
         brain = rowMeans(select(., all_of(passed.p))),
         cognition = rowMeans(select(., starts_with(c('nih', 'pea', 'lmt'))))) %>%
  select(-all_of(passed.p), -starts_with(c('nih', 'pea', 'lmt'))) %>% 
  mutate(across(where(is.numeric), ~as.vector(scale(., scale = T))))

covs_y2 <- c('mode_of_delivery', 'birth_complications_n', 'pregnancy_complications', 'premature_wk', 'handedness_c', 'education_max', 'child_sex', 'child_age', 'vision')

summary(mediator_modely2 <- lm(brain ~ bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + handedness_c + premature_wk + education_max + child_sex + child_age + vision, med.y2))
summary(outcome_modely2 <- lm(cognition ~ brain + bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + handedness_c + premature_wk + education_max + child_sex + child_age + vision, med.y2))

confint(outcome_modely2, 'bfq_duration')     # lactation on cognition
confint(mediator_modely2, 'bfq_duration')    # lactation on brain
confint(outcome_modely2, 'brain')    

set.seed(4994)

mediation_resultsy2 <- mediation::mediate(mediator_modely2, outcome_modely2, treat = 'bfq_duration', mediator = 'brain', covariates = covs_y2, boot = T, sims = 5000)

(medy2 <- extract_mediation_summary(summary(mediation_resultsy2)))

write.csv(medy2, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/ABCD_medy2.csv'))

