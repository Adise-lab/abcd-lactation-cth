#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for preparing the data from the ABCD data release 5.1 for the paper entitled #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Aug 14th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, summarytools)

#2# Set paths

datadir <- '/Users/jongonzalez/Desktop/abcd.nosync/'

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

(original.sample.n <- length(abcd.excl <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_screen.csv'))$src_subject_id))      # initial sample size

# Create a character vector of subject IDs (src_subject_id) that passed ABCD minimal exclusion criteria

neurol <- c('scrn_cpalsy', 'scrn_tumor', 'scrn_stroke', 'scrn_aneurysm', 'scrn_hemorrhage', 'scrn_hemotoma', 'scrn_tbi')   # vector of variables that will be later added for exclusion

abcd.excl <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_screen.csv')) %>%                                
  filter(scrn_gestage == 0,             # no extremely premature (>28 weeks born early) --> 69 excluded (3 NA)
         scrn_schiz == 0,               # exclude those with schizophrenia --> 3 excluded (37 NA)
         scrn_asd == 0,                 # exclude those with autism spectrum disorder --> 201 excluded (39 NA)
         scrn_intdisab == 0,            # exclude those with intellectual dissabilities --> 13 excluded (38 NA)
         scrn_birthcomp == 0,           # exclude those with birth complications that required >30d hospitalization
         scrn_birthwt == 0) %>%         # exclude those with low birth weight (<1,200 gr)
  mutate(scrn_tbi_scan_excl = ifelse(scrn_tbi_scan_excl > 0, 0, 1),
         scrn_tbi = ifelse(rowSums(select(., c(scrn_tbi_loc, scrn_tbi_mem, scrn_tbi_scan_excl)), na.rm = T) > 0, 0, 1),    # will be reversed later (originally, 0 = yes); na.rm = T bc when inverted 0 transforms to NA
         scrn_epls = ifelse(is.na(scrn_epls), 0, scrn_epls),                                                               # replace NAs by 0, don't know why these are not filled
         scrn_abn_scan = ifelse(is.na(scrn_abn_scan), 0, scrn_abn_scan)) %>%                                               # if NA means there was no finding (1 is yes, 2 is unsure; the ABCD wiki is wrong -> 0 = yes)
  mutate_at(vars(all_of(neurol)), ~as.numeric(!.)) %>%                                                                     # reverse all neurological scores so now yes is 1 (and not 0)
  select(src_subject_id, eventname, scrn_abn_scan, all_of(neurol)) %>% 
  mutate(neurol_excl = rowSums(select(., all_of(neurol)))) %>%    # now sum all neurological scores                                                       
  filter(neurol_excl == 0,              # exclude those with ANY neurological disorder (neurol scores > 0)
         scrn_abn_scan == 0) %>%        # exclude those with abnormal scan findings
  select(src_subject_id) %>% 
  pull()

(n.excluded <- original.sample.n - length(abcd.excl))   # this is the amount of individuals excluded with ABCD minimal criteria

length(abcd.excl)     # this is your sample after ABCD minimal exclusion

# Generate summary numbers of how many subjects met study-specific exclusion criteria -> 1. Non-biological mothers and 2. Drug use during pregnancy

fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>%
  filter(eventname == 'baseline_year_1_arm_1', 
         src_subject_id %in% abcd.excl) %>%
  mutate_at(vars(devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack), ~replace_na(., 1)) %>% 
  summarise(non_biological = sum(devhx_1_p == 0),
            alcohol_during = sum(devhx_9_alcohol == 1),
            tobacco_during = sum(devhx_9_tobacco == 1),
            marijuana_during = sum(devhx_9_marijuana == 1),
            opioids_during = sum(devhx_9_oxycont == 1),
            opiates_during = sum(devhx_9_her_morph == 1),
            coc_crack_during = sum(devhx_9_coc_crack == 1))

# Now extract the subjects who DID NOT MEET this exclusion criteria

study.excl <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>%
  filter(src_subject_id %in% abcd.excl) %>%
  mutate_at(vars(devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack), ~replace_na(., 1)) %>% 
  filter(eventname == 'baseline_year_1_arm_1',     # focus on baseline data (not 4y assessment) to exclude
         devhx_1_p == 1,                           # biological mother only (1082 non-biological, 0 NA)
         devhx_9_alcohol == 0,                     # 0 use of alcohol AFTER knowing pregnancy (238 used)
         devhx_9_tobacco == 0,                     # 0 use of tobacco AFTER knowing pregnancy (404 used)
         devhx_9_marijuana == 0,                   # 0 se of marihuana AFTER knowing pregnancy (147 used)
         devhx_9_oxycont == 0,                     # 0 use of opiods AFTER knowing pregnancy (15 used)
         devhx_9_her_morph == 0,                   # 0 use of opiates AFTER knowing pregnancy (8 used)
         devhx_9_coc_crack == 0) %>%               # 0 use of crack/cocaine AFTER knowing pregnancy (29 used)
  select(src_subject_id) %>% 
  pull()

(n.excluded2 <- length(abcd.excl) - length(study.excl))    # this is the amount of individuals excluded with study-specific criteria (ABCD already excluded)

length(study.excl)    # this is your sample

#3# Load behavioral/demographic data -- note: additional exclusion applied (missing data, maternal drug use while pregnant, etc)

df <-                              
  # append parent education
  fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_demo.csv')) %>%
  #   # only retain baseline and 2y observations & remove subjects that meet exclusion criteria 
  filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1'),
         src_subject_id %in% abcd.excl) %>%
  mutate(child_sex = factor(ifelse(demo_sex_v2 == 1, 'Male', 
                                   ifelse(demo_sex_v2 == 2, 'Female', NA)))) %>%    # 1 = male, 2 = female, 3/4 = intersex, 999 = don't know, 777 = refuse 
  # group by subject to fill not filled observations that should be invariant
  group_by(src_subject_id) %>% 
  # fill assigned sex-at-birth by subject in updown direction (asked at baseline, need to be extended to 2y)
  fill(c(child_sex, race_ethnicity, demo_prnt_prtnr_v2, demo_prnt_prtnr_bio), .direction = 'updown') %>% 
  # fill parent education by subject in downup direction (asked at 2y, need to be extended to baseline)
  fill(c(demo_prnt_ed_v2_2yr_l, demo_prtnr_ed_v2_2yr_l), .direction = 'downup') %>%         # demo_prnt_ed_v2_2yr_l = highest grade acheived; demo_prtnr_ed_v2_2yr_l = highest grade partner achieved 
  ungroup() %>% 
  mutate(race_ethnicity = as.factor(case_when(race_ethnicity == 1 ~ 'White',
                                              race_ethnicity == 2 ~ 'Black',
                                              race_ethnicity == 3 ~ 'Hispanic',
                                              race_ethnicity == 4 ~ 'Asian',
                                              TRUE ~ 'Other')),
         # fix inconsistencies to calculate CONTINOUS education for the parent (prnt) and partner (prtnr) so PhD (21y) is the max
         education_parent1 = ifelse(demo_prnt_ed_v2_2yr_l == 22, 15,                  # 22 = <1 year of college (or less than 10 classes) - should equal to 15 yrs of education
                                    ifelse(demo_prnt_ed_v2_2yr_l == 23, 16,           # 23 = 1<= year of college but no degree - this should equal 16 yrs of education
                                           ifelse(demo_prnt_ed_v2_2yr_l == 999 |      # 999 = Do not know / 777 = Refuse to answer -> treat as NA
                                                    demo_prnt_ed_v2_2yr_l == 777, NA, demo_prnt_ed_v2_2yr_l))),    
         education_parent2 = ifelse(demo_prtnr_ed_v2_2yr_l == 22, 15,
                                    ifelse(demo_prtnr_ed_v2_2yr_l == 23, 16, 
                                           ifelse(demo_prtnr_ed_v2_2yr_l == 999 | 
                                                    demo_prtnr_ed_v2_2yr_l == 777, NA, demo_prtnr_ed_v2_2yr_l))),
         # find the maximum value in education and assign that to the household
         education_max = ifelse(!is.na(education_parent1) & is.na(education_parent2), education_parent1,
                                ifelse(is.na(education_parent1) & !is.na(education_parent2), education_parent2,
                                       ifelse(!is.na(education_parent1) & !is.na(education_parent2), pmax(education_parent1, education_parent2), NA))),
         # Create 5-level factor degree only for reporting according to Adise scripts
         education_max_f = factor(case_when(education_max < 13 ~ 'HS',
                                               between(education_max, 13, 14) ~ 'HS/GD',
                                               between(education_max, 15, 17) ~ 'Some college',
                                               education_max == 18 ~ 'BS degree',
                                               between(education_max, 19, 21) ~ 'Postgraduate degree',
                                               TRUE ~ NA), level = c('HS', 'HS/GD', 'Some college', 'BS degree', 'Postgraduate degree'))) %>%    
  # retain only subjects with data available at both baseline and 2y
  select(src_subject_id, eventname, child_sex, race_ethnicity, education_max, education_max_f) %>% 
  # append child's age at visit (months), site and family ID
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_y_lt.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>% 
               rename(site = site_id_l,
                      child_age = interview_age) %>% 
               group_by(src_subject_id) %>% 
               # extend rel_family_id to 2y follow-up
               fill(rel_family_id, .direction ='updown') %>% 
               ungroup() %>% 
               # recode visit type for interpretation
               mutate(visit_type = as.factor(case_when(visit_type == 1 ~ 'In-person',
                                                       visit_type == 2 ~ 'Remote',
                                                       visit_type == 3 ~ 'Hybrid',
                                                       TRUE ~ NA))) %>% 
               select(src_subject_id, eventname, visit_type, site, rel_family_id, child_age),
             by = c('src_subject_id', 'eventname')) %>% 
  # append handedness 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_ehis.csv')) %>% 
               filter(eventname == 'baseline_year_1_arm_1') %>% 
               mutate(handedness_raw = rowMeans(select(., c('ehi1b', 'ehi2b', 'ehi3b', 'ehi4b')), na.rm = T),   # continous variable designating handedness (-100 left, 0 ambi, 100 right)
                      handedness_c = scales::rescale(handedness_raw, from = c(-100, 100), to = c(-1, 1)),       # rescale to -1 to 1 to avoid large values and potential fit problems
                      handedness_f = factor(case_when(ehi_y_ss_scoreb == 1 ~ 'Right',                           # recode handedness factor variable for interpretation
                                                    ehi_y_ss_scoreb == 2 ~ 'Left', 
                                                    TRUE ~ 'Mixed'))) %>% 
              select(src_subject_id, handedness_c, handedness_f), 
             by = 'src_subject_id') %>% 
  # append data to be excluded and important info
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>%
               filter(eventname == 'baseline_year_1_arm_1',     # focus on baseline data (not 4y assessment) to exclude
                      src_subject_id %in% study.excl) %>% 
                mutate_at(vars(devhx_10a3_p:devhx_10m3_p), ~na_if(., 999)) %>%                                     # if 999, transform to NA -> variables for pregnancy problems
                mutate_at(vars(devhx_14a3_p:devhx_14b3_p), ~na_if(., 999)) %>%                                     # if 999, transform to NA -> variables for birth problems
                mutate_at(vars(devhx_8_tobacco, devhx_8_alcohol, devhx_8_marijuana, devhx_8_oxycont, devhx_8_her_morph, devhx_8_coc_crack), ~na_if(., 999)) %>%         # if 999, transform to NA -> variables for illicit drugs before knowing pregnancy
                mutate(mode_of_delivery = factor(ifelse(devhx_13_3_p == 1, 'Cesarean', ifelse(devhx_13_3_p == 0, 'Vaginal', NA))),           # if not 1 or 0, means it's 999 --> set to NA to further exclude
                       premature = as.factor(ifelse(devhx_12a_p == 999, NA, devhx_12a_p)),                                                   # if don't know (999) set to NA to further exclude
                       # if premature is 0 then premature_wk is 0, if premature is 1 but don't know (999) OR > than 12 wk (extremely premature, should have answered 'yes' in the scnr but some didn't) -> NA to further exclude
                       premature_wk = ifelse(premature == 0, 0, ifelse(premature == 1 & devhx_12_p == 999 | premature == 1 & devhx_12_p > 12, NA, devhx_12_p)),           
                       pregnancy_complications = as.factor(ifelse(rowSums(select(., c(devhx_10a3_p:devhx_10m3_p))) > 0, 1, 0)),  # if ANY complication during pregnancy, 1
                       birth_complications_n = as.factor(ifelse(rowSums(select(., c(devhx_14a3_p:devhx_14b3_p))) > 0, 1, 0)),    # if ANY complication during birth, 1
                       tobacco = as.factor(ifelse(devhx_8_tobacco == 1, 1, ifelse(devhx_8_tobacco == 0, 0, NA))),              # use of tobacco BEFORE pregnancy
                       alcohol = as.factor(ifelse(devhx_8_alcohol == 1, 1, ifelse(devhx_8_alcohol == 0, 0, NA))),              # use of alcohol BEFORE  pregnancy
                       marihuana = as.factor(ifelse(devhx_8_marijuana == 1, 1, ifelse(devhx_8_marijuana == 0, 0, NA))),          # use of cannabis BEFORE  pregnancy
                       other_drugs = as.factor(ifelse(rowSums(select(., c(devhx_8_oxycont, devhx_8_her_morph, devhx_8_coc_crack))) > 0, 1, 0)),          # if use of ANY illicit drug, excluding marihuana, 1
                       illicit_drugs = as.factor(ifelse(rowSums(select(., c(devhx_8_marijuana, devhx_8_oxycont, devhx_8_her_morph, devhx_8_coc_crack))) > 0, 1, 0)),      # if use of ANY illicit drug, including marihuana, 1
                       bfq_months_y0 = ifelse(devhx_18_p > 72, NA, devhx_18_p),                   # remove imposible values (182, 184 months of lactation --> likely 82 and 84 but who knows)
                       bfq_duration_y0 = case_when(bfq_months_y0 == 0 ~ 0,                  # 0 lactation
                                                between(bfq_months_y0, 1, 6) ~ 1,        # 1 to 6 months
                                                between(bfq_months_y0, 7, 12) ~ 2,       # 7 to 12 months
                                                bfq_months_y0 > 12 ~ 3,                  # > 12 months
                                                TRUE ~ NA)) %>%
                select(src_subject_id, bfq_months_y0, bfq_duration_y0, mode_of_delivery, premature, premature_wk, birth_complications_n, pregnancy_complications, tobacco, alcohol, marihuana, other_drugs, illicit_drugs),
             by = 'src_subject_id') %>% 
  # append 3y assessment on lactation to cross-reference
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_bfq.csv')) %>% 
               filter(bfq_breastfeed_p != 999) %>%
               mutate(bfq_breastfeed_p = as.factor(bfq_breastfeed_p),
                      bfq_duration_y3 = ifelse(bfq_breastfeed_p == 0, 0, 
                                               ifelse(bfq_breastfeed_time_p == 999, NA, bfq_breastfeed_time_p))) %>%    # if bfq_breastfeed_p was no (0), replace in bfq_breastfeed_time_p's NA by 0 (0 times)
               select(src_subject_id, bfq_duration_y3), 
             by = 'src_subject_id') %>% 
  # append PDS; take the average of parent and youth report if possible; if not only parent or youth 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_pds.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>% 
               left_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_y_pds.csv')), by = c('src_subject_id', 'eventname')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>% 
               mutate(pds_female = case_when(!is.na(pds_y_ss_female_category_2) & !is.na(pds_p_ss_female_category_2) ~ (pds_y_ss_female_category_2 + pds_p_ss_female_category_2) / 2,
                                             !is.na(pds_y_ss_female_category_2) & is.na(pds_p_ss_female_category_2) ~ pds_y_ss_female_category_2,
                                             is.na(pds_y_ss_female_category_2) & !is.na(pds_p_ss_female_category_2) ~ pds_p_ss_female_category_2,
                                             TRUE ~ NA_real_),
                      pds_male = case_when(!is.na(pds_y_ss_male_cat_2) & !is.na(pds_p_ss_male_category_2) ~ (pds_y_ss_male_cat_2 + pds_p_ss_male_category_2) / 2,
                                           !is.na(pds_y_ss_male_cat_2) & is.na(pds_p_ss_male_category_2) ~ pds_y_ss_male_cat_2,
                                           is.na(pds_y_ss_male_cat_2) & !is.na(pds_p_ss_male_category_2) ~ pds_p_ss_male_category_2,
                                           TRUE ~ NA_real_),
                      pds = coalesce(pds_female, pds_male),
                      eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               select(src_subject_id, eventname, pds),
             by = c('src_subject_id', 'eventname')) %>% 
  mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>%
  droplevels()

df.split <- split(df, df$eventname)

count_missing <- function(data) {
  sapply(data, function(x) sum(is.na(x)))
}
(missing_counts_by_event <- lapply(df.split, count_missing))

before.na.rm <- df %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) %>% 
  select(n) %>% 
  pull()

df <- df %>% 
  drop_na()

after.na.rm <- df %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) %>% 
  select(n) %>% 
  pull()

before.na.rm;after.na.rm;before.na.rm[1:2]-after.na.rm[1:2]

df %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))

# this is the classification of bf duration based on baseline assessment
df %>% 
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  group_by(bfq_duration_y0) %>% 
  freq(bfq_months_y0)

# this cross-reference baseline (devhx_18_p, months) and y3 reports on lactation duration (bfq_duration, 4-level)
df %>% 
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  mutate(flag = case_when(bfq_duration_y3 == 0 & bfq_months_y0 == 0 ~ 'no',                      # if bf 3y is 'No lactation' (0) and bf 0y is 0, DO NOT flag
                          between(bfq_duration_y3, 1, 2) & between(bfq_months_y0, 1, 3) ~ 'no',  # if bf 3y is between 'Several days' (1) and '1-3 months' (2) and bf 0y is between 1 and 3 months, DO NOT flag
                          bfq_duration_y3 == 3 & between(bfq_months_y0, 4, 6) ~ 'no',            # if bf 3y is between '4-6 months' (3) and bf 0y is between 4 and 6 months, DO NOT flag
                          bfq_duration_y3 == 4 & between(bfq_months_y0, 7, 9) ~ 'no',            # if bf 3y is between '7-9 months' (4) and bf 0y is between 7 and 9 months, DO NOT flag
                          bfq_duration_y3 == 5 & between(bfq_months_y0, 10, 12) ~ 'no',          # if bf 3y is between '10-12 months' (5) and bf 0y is between 10 and 12 months, DO NOT flag
                          bfq_duration_y3 == 6 & between(bfq_months_y0, 13, 18) ~ 'no',          # if bf 3y is between '13-18 months' (6) and bf 0y is between 13 and 18 months, DO NOT flag
                          bfq_duration_y3 == 7 & between(bfq_months_y0, 19, 24) ~ 'no',          # if bf 3y is between '19-24 months' (7) and bf 0y is between 19 and 24 months, DO NOT flag
                          bfq_duration_y3 == 8 & bfq_months_y0 > 24 ~ 'no',                      # if bf 3y is between '>24 months' (8) and bf 0y is greater than 24 months, DO NOT flag
                          TRUE ~ 'yes')) %>% 
  filter(flag == 'yes') %>% 
  select(src_subject_id) %>% 
  pull() -> excl2

length(excl2)    # this is the n of people with inconsistent reporting between baseline and y3

# remove those participants with inconsistent findings
df <- df %>%     
  filter(!src_subject_id %in% excl2) %>% 
  mutate(bfq_duration_y3 = case_when(bfq_duration_y3 == 0 ~ 0,
                                     between(bfq_duration_y3, 1, 3) ~ 1,
                                     between(bfq_duration_y3, 4, 5) ~ 2,
                                     bfq_duration_y3 > 5 ~ 3)) %>% 
  rename(bfq_duration = bfq_duration_y0)

# confirm the two metrics are now equal
cor.test(df$bfq_duration[df$eventname == 'baseline_year_1_arm_1'], df$bfq_duration_y3[df$eventname == 'baseline_year_1_arm_1'])
plot(df$bfq_duration[df$eventname == 'baseline_year_1_arm_1'] ~ df$bfq_duration_y3[df$eventname == 'baseline_year_1_arm_1'])

df %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))  

#4# Load and append MRI data -- note: more exclusion criteria applied (bad scans, incidental findings)

df2 <- df %>% 
  # append mri id
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_adm_info.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1'),
                      mri_id = factor(mri_info_deviceserialnumber),
                      mri_date = ymd(mri_info_studydate)) %>% 
               select(src_subject_id, eventname, mri_id, mri_date), by = c('src_subject_id', 'eventname')) %>% 
  # append incidental findings of the images (mrif_score) -- bear in mind that you have to decide whether to drop an entire subject based on bad scan at baseline/follow-up, or just the timepoint with the bad scan
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_qc_clfind.csv')) %>%
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>%
               filter(between(mrif_score, 1, 2)) %>%                                              # filter only those no abnormal (1) or normal variation (2) --> 581 excluded baseline and 627 at follow up (304 also baseline)
               select(src_subject_id, eventname, mrif_score), by = c('src_subject_id', 'eventname')) %>%
  # T1w and T2w passed ALL criteria to be included (see https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html)
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_qc_incl.csv')) %>%
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1'),
                      imgincl_t1w_include == 1,                                                 # filter only usable t1w --> 227 at baseline and 501 follow up excluded (57 also baseline)
                      imgincl_t2w_include == 1) %>%                                           # filter only usable t2w --> 600 at baseline and 1217 follow up excluded (182 from baseline)
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>%
               select(src_subject_id, imgincl_t1w_include, imgincl_t2w_include, eventname), by = c('src_subject_id', 'eventname')) %>%
  # append freesurfer data (CT data, desterieux)
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_thk_dst.csv')) %>%      #mri_y_smr_area_dst.csv mri_y_smr_thk_dst.csv
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               select(src_subject_id, eventname, !ends_with(c('_149', '_150', '_151'))),
             by = c('src_subject_id', 'eventname')) %>%
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_t1_contr_dst.csv')) %>%      # mri_y_smr_t1_contr_dst.csv --> higher t1/t2 ratio fluid (t2) < fat (t1); mri_y_smr_t1_gray_dst.csv or mri_y_smr_t1_white_dst.csv
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               rename_with(.fn = ~paste0("mrisdp_", seq_along(.), '_myelin'), .cols = starts_with("mrisdp_")) %>%    # mrisdp or dmri
               select(src_subject_id, eventname, !ends_with(c('149_myelin', '150_myelin', '151_myelin'))), 
             by = c('src_subject_id', 'eventname')) %>% 
  drop_na() 

# this is your sample after meeting MRI criteria (mrif_score, imgincl_t1w_include, imgincl_t2w_include)
df2 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) 

# create a vector of participants with data at y2 but not at baseline (bc of exclusion criteria applied)
length(not.at.baseline <- setdiff(df2$src_subject_id[df2$eventname == '2_year_follow_up_y_arm_1'], 
                                  df2$src_subject_id[df2$eventname == 'baseline_year_1_arm_1'])) 

# remove subjects who had no baseline data
df2 <- df2 %>% 
  filter(!src_subject_id %in% not.at.baseline)

df2 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) 

#5# Picking one sibling at random

set.seed(4994)       # set seed reproducibility

unique_participants <- df2 %>%
  # Determine the number of unique time points for each participant
  group_by(src_subject_id, rel_family_id) %>%
  summarise(time_points = n_distinct(eventname), .groups = "drop") %>%
  # Flag participants with both baseline and follow-up data
  mutate(has_both = ifelse(time_points == 2, 1, 0)) %>%
  # Group participants by family
  group_by(rel_family_id) %>%
  # Prioritize: randomly select one from those with both time points if available
  filter(has_both == max(has_both)) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  pull(src_subject_id)

length(unique_participants)   # number of unique participants (of 3,970)

# remove random sibling
df3 <- df2 %>%
  filter(src_subject_id %in% unique_participants)

# confirmation
df3 %>%         
  group_by(rel_family_id) %>%
  summarize(count = n()) %>%
  ungroup() %>% 
  filter(count > 2) %>% 
  nrow()                        # this should result in 0; no participants with > 2 rel_family_id (no siblings)

df3 %>%                         # numbers will be the same regardless of seed but want to ensure we pick the same sibling if data is good on both
  group_by(eventname) %>%
  summarize(n = length(src_subject_id))

#6# Remove MRI device effects from cortical thickness and t1w/t2w data (Destrieux)

library(longCombat)

# create a vector of the variable names  (149:151 -> mean and total thickness/myelin)
(featurenames <- names(df3 %>% select(., starts_with('mrisdp'))))

length(featurenames) # destrieux atlas is 148 regions x 2 (CT and myelin) = 296

# make sure to center your data as in the main analysis before running longCombat (it also runs a LME, formula must be the same as the analysis)
covs2scale <- c('child_age', 'bfq_duration', 'premature_wk', 'education_max', 'handedness_c')

tmp <- df3 %>% 
  mutate_at(vars(all_of(covs2scale)), ~scale(., scale = F))

# apply longCombat
simdata_combat <- longCombat(idvar = 'src_subject_id', 
                             timevar = 'child_age',        # time variable is not eventname but (centered) age
                             batchvar = 'mri_id',          # batch effect is MRI serial device ID, not site
                             features = featurenames,      # featurenames = MRI features to remove MRI effects from
                             formula = 'child_age * bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + premature_wk + education_max + handedness_c + child_sex',   
                             ranef= '(1 | src_subject_id)',   # random intercept for subject bc repeated measurements
                             data = tmp)

df3 %>% 
  droplevels() %>% 
  select(., -all_of(featurenames)) %>%       # this will remove the original features from the dataframe 
  cbind(., simdata_combat$data_combat %>%    # this is the dataframe w/features resulting from longCombat 
          rename_at(vars(all_of(paste0(featurenames, '.combat'))), ~str_remove(., '.combat')) %>%    # we need to remove the '.combat' suffix
          select(., all_of(featurenames))) -> df.mri     # update w/combatt'd features (that now match names in featurenames)

# check very extreme values in numeric
df.mri %>% 
  group_by(eventname) %>% 
  mutate(across(where(is.numeric), ~scale(., scale = T))) %>% 
  summarise(across(where(is.numeric), ~sum(abs(.) > 4, na.rm = T))) %>%
  ungroup() %>%   
  View()
 
write.csv(df.mri, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data.csv'), row.names = F)

#7# Append cognition

df4 <- df.mri %>%  
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_svs.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               rename(vision = snellen_va_y) %>% 
               select(src_subject_id, eventname, vision) %>% 
               drop_na(),
             by = c('src_subject_id', 'eventname')) %>% 
  # append raw crystalized and fluid uncorrected scores from NIH toolbox
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_nihtb.csv')) %>%
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               group_by(eventname) %>%
               mutate_at(vars(nihtbx_flanker_uncorrected, nihtbx_pattern_uncorrected, nihtbx_picture_uncorrected), ~scales::rescale(., to = c(0, 19))) %>%
               ungroup() %>%
               mutate(nih_fluid = rowMeans(select(., c(nihtbx_flanker_uncorrected, nihtbx_pattern_uncorrected, nihtbx_picture_uncorrected)), na.rm = T)) %>%
               select(src_subject_id, eventname, nih_fluid) %>%
               drop_na(), by = c('src_subject_id', 'eventname')) %>% 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_ravlt.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               # mutate(pea_ravlt_learning_tc = rowSums(select(., pea_ravlt_sd_trial_i_tc, pea_ravlt_sd_trial_ii_tc, pea_ravlt_sd_trial_iii_tc, pea_ravlt_sd_trial_iv_tc, pea_ravlt_sd_trial_v_tc))) %>% 
               select(., src_subject_id, eventname, pea_ravlt_ld_trial_vii_tc) %>% 
               drop_na(), by = c('src_subject_id', 'eventname')) %>% 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_lmt.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               select(., src_subject_id, eventname, lmt_scr_perc_correct) %>% 
               mutate(lmt_scr_perc_correct = ifelse(lmt_scr_perc_correct > 1, lmt_scr_perc_correct / 100, lmt_scr_perc_correct)) %>%    # there are some inconsistencies that I believe are due to someone reporting percentage instead of 0.
               drop_na(), by = c('src_subject_id', 'eventname')) %>%
  mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) 

df4 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))

df4 <- df4 %>%   # alternatively get rid of remote/hybrid evaluations
  filter(visit_type == 'In-person')

df4 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))


length(not.at.baseline <- setdiff(df4$src_subject_id[df4$eventname == '2_year_follow_up_y_arm_1'], 
                                  df4$src_subject_id[df4$eventname == 'baseline_year_1_arm_1']))

length(both.timepoints <- intersect(df4$src_subject_id[df4$eventname == '2_year_follow_up_y_arm_1'], 
                                    df4$src_subject_id[df4$eventname == 'baseline_year_1_arm_1']))

df4 <- df4 %>%
  filter(!src_subject_id %in% not.at.baseline)
# filter(src_subject_id %in% both.timepoints)

df4 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))

# table for report
vtable::st(df4 %>% 
             mutate(bfq_duration = factor(case_when(bfq_duration == 0 ~ 'No breastfed',
                                                    bfq_duration == 1 ~ '1-6 months',
                                                    bfq_duration == 2 ~ '7-12 months',
                                                    bfq_duration == 3 ~ '>12 months'), 
                                          levels = c('No breastfed', '1-6 months', '7-12 months', '>12 months')),
                    premature_wk = na_if(premature_wk, 0)), 
           vars = c('bfq_duration', 'devhx_18_p', 'mode_of_delivery', 'alcohol', 'tobacco', 'marihuana', 'illicit_drugs', 'loweight',
                    'birth_complications', 'birth_complications_n', 'pregnancy_complications', 'premature', 'premature_wk',
                    'education_max_f', 'education_max', 'handedness_f', 'child_sex', 'pds', 'child_age', 'race_ethnicity', 'vision'), 
           group = 'eventname')

(featurenames2 <- names(df4 %>% select(starts_with(c('nih', 'pea', 'lmt')))))

covs2scale2 <- c('child_age', 'bfq_duration', 'premature_wk', 'handedness_c', 'education_max', 'vision')

tmp2 <- df4 %>% 
  mutate_at(vars(all_of(covs2scale2)), ~scale(., scale = F))

simdata_combat <- longCombat(idvar = 'src_subject_id', 
                             timevar = 'child_age',
                             batchvar = 'site', 
                             features = featurenames2, 
                             formula = 'child_age * bfq_duration + mode_of_delivery + birth_complications_n + pregnancy_complications + premature_wk + education_max + handedness_c + child_sex + vision',
                             ranef= '(1| src_subject_id)',
                             data = tmp2)

df4 %>% 
  droplevels() %>% 
  select(., -matches(featurenames2)) %>%
  cbind(., simdata_combat$data_combat %>% 
          rename_at(vars(matches(featurenames2)), ~str_remove(., '.combat')) %>% 
          select(., matches(featurenames2))) -> df.cog

write.csv(df.cog, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv'), row.names = F)

