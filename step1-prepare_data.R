#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for preparing the data from the ABCD data release 5.1 for the paper entitled #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Nov 19th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, summarytools)

#2# Set paths

datadir <- '/Users/jongonzalez/Desktop/abcd.nosync/'

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

# Save the number of the original sample before any exclusion criteria is applied to help follow the process

(original.sample.n <- length(fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_screen.csv'))$src_subject_id))      # initial sample size

#3# Initial exclusion criteria (ABCD study general)

# Create a character vector of subject IDs (src_subject_id) that passed ABCD minimal exclusion criteria (look up variables for meaning at https://data-dict.abcdstudy.org/?)

reverse.code <- c('scrn_cpalsy', 'scrn_tumor', 'scrn_stroke', 'scrn_aneurysm', 'scrn_hemorrhage', 'scrn_hemotoma', 'scrn_percept_corr')    # these variables are exclusionary

# This will only give you a summary of how many participants there are who met exclusionary criteria 

fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_screen.csv')) %>% 
  # this line will reverse code exclusionary variables so 0 (originally meant STOP evaluation) now is 1
  mutate_at(vars(all_of(reverse.code)), ~as.numeric(!.)) %>%
  # these are different because it is not 0 or 1 but 0, 1, 2 (0 yes and stop, 1 no, and 2 unsure) -- if 1 or 2 (no or unsure) recode as 0 (no)
  mutate(scrn_abn_scan = ifelse(scrn_abn_scan == 0, 1, 0),
         scrn_birthcomp_excl = ifelse(scrn_birthcomp_excl == 0, 1, 0),
         scrn_tbi_scan_excl = ifelse(scrn_tbi_scan_excl == 0, 1, 0)) %>% 
  # transform NAs to 0 just for the summaries
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>% 
  # sum all reversed neurological scores so now yes is 1 (and not 0)
  summarise(scrn_gestage = sum(scrn_gestage == 1),                  # extreme prematurity (>28wk), n = 0
            scrn_schiz = sum(scrn_schiz == 1),                      # schizophrenia, n = 3
            scrn_asd = sum(scrn_asd == 1),                          # autism spectrum disorder, n = 201
            scrn_intdisab = sum(scrn_intdisab == 1),                # intellectual disability, n = 13
            scrn_birthcomp_excl = sum(scrn_birthcomp_excl == 1),    # birth complications with >30d hospitalization, n = 0
            scrn_birthwt = sum(scrn_birthwt == 1),                  # extreme low birth weight (<1200gr), n = 0
            scrn_seizure = sum(scrn_seizure == 1),                  # number of last month's seizures while in medication, n = 0
            scrn_cpalsy = sum(scrn_cpalsy == 1),                    # cerebral palsy, n = 4
            scrn_tumor = sum(scrn_tumor == 1),                      # brain tumor, n = 3
            scrn_stroke = sum(scrn_stroke == 1),                    # stroke, n = 3
            scrn_aneurysm = sum(scrn_aneurysm == 1),                # aneurysm, n = 2
            scrn_hemotoma = sum(scrn_hemotoma == 1),                # brain hematoma, n = 0
            scrn_hemorrhage = sum(scrn_hemorrhage == 1),            # brain hemorrhage, n = 5
            scrn_abn_scan = sum(scrn_abn_scan == 1),                # history of neurological abnormalities in brain, n = 0
            scrn_tbi_loc = sum(scrn_tbi_loc == 1),                  # loss of consciouscness due to traumatic brain injury, n = 1
            scrn_tbi_mem = sum(scrn_tbi_mem == 1),                  # memory loss that lasted >1 day due to traumatic brain injury, n = 0
            scrn_tbi_scan_excl = sum(scrn_tbi_scan_excl == 1),      # abnormalities in brain scan after traumatic brain injury, n = 0
            scrn_percept_corr = sum(scrn_percept_corr == 1))        # visual/auditory problems that cannot be corrected with the use of glasses/hearing aid/others, n = 0

# Now this will save a list of participants (abcd.excl) who met ANY exclusionary criteria stated above to be further excluded

abcd.excl <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_screen.csv')) %>%
  # we have to reverse code again as before
  mutate_at(vars(all_of(reverse.code)), ~as.numeric(!.)) %>%      
  mutate(scrn_abn_scan = ifelse(scrn_abn_scan == 0, 1, 0),
         scrn_birthcomp_excl = ifelse(scrn_birthcomp_excl == 0, 1, 0),
         scrn_tbi_scan_excl = ifelse(scrn_tbi_scan_excl == 0, 1, 0)) %>%
  # 1 = YES means EXCLUDE
  filter(scrn_gestage == 1 |
           scrn_schiz == 1 |
           scrn_asd == 1 |          
           scrn_intdisab == 1 |
           scrn_birthcomp_excl == 1 |
           scrn_birthwt == 1 |
           scrn_seizure == 1 |
           scrn_cpalsy == 1 |
           scrn_tumor == 1 |
           scrn_stroke == 1 |
           scrn_aneurysm == 1 |
           scrn_hemotoma == 1 |
           scrn_hemorrhage == 1 |
           scrn_abn_scan == 1 |
           scrn_tbi_loc == 1 |
           scrn_tbi_mem == 1 |
           scrn_tbi_scan_excl == 1 |
           scrn_percept_corr == 1) %>%
  select(src_subject_id) %>%
  pull()

length(abcd.excl)     # this is the number of people to exclude with ABCD minimal criteria + ASD (n = 231)

(n.excluded <- original.sample.n - length(abcd.excl))   # this is how the sample looks like after applying exclusion criteria (11,867 - 231)

#4# Study-specific exclusionary criteria

# Generate summary numbers of how many subjects met study-specific exclusion criteria -> 1. Non-biological mothers, 2. Drug use during pregnancy, etc.

fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>% 
  filter(eventname == 'baseline_year_1_arm_1',          # focus on baseline
         !src_subject_id %in% abcd.excl) %>%            # only look on subjects that were not excluded in previous step
  mutate(weight_kg1 = birth_weight_lbs * 453.6,         # transform pounds to grams
         weight_kg2 = birth_weight_oz * 28.35,          # transform ounces to grams
         birthweight_gr = ifelse(!is.na(weight_kg1) & !is.na(weight_kg2), weight_kg1 + weight_kg2, ifelse(!is.na(weight_kg1) & is.na(weight_kg2), weight_kg1, NA))) %>%   # if weight2 is NA is bc there was no ounces, keep weight1
  mutate_at(vars(devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack, devhx_12_p, devhx_12a_p), ~na_if(., 999)) %>%      # if response in any is 999 (Do not know/Do not remember) set to NA
  # transform NAs to 1 = if they do not answer is like they 'used' --> it is important that we are sure that kids were not exposed to drugs during pregnancy. If mom's do not remember that's a bad sign
  mutate_at(vars(devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack, devhx_12_p, devhx_12a_p), ~replace_na(., 1)) %>%   
  select(src_subject_id, eventname, devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack, devhx_12a_p, devhx_12_p, birthweight_gr) %>% 
  # count the number of subjects that would meet study-specific exclusionary criteria
  summarise(non_biological = sum(devhx_1_p == 0),                                        # n = 1,585 were not the biological mom - we need to secure breastfeeding reports and unfortunately this means only including biological moms
            alcohol_during = sum(devhx_9_alcohol == 1),                                  # n = 597 used alcohol after knowing of pregnancy
            tobacco_during = sum(devhx_9_tobacco == 1),                                  # n = 862 used tobacco after knowing of pregnancy
            marijuana_during = sum(devhx_9_marijuana == 1),                              # n = 512 used marijuana after knowing of pregnancy
            opioids_during = sum(devhx_9_oxycont == 1),                                  # n = 322 opioids
            opiates_during = sum(devhx_9_her_morph == 1),                                # n = 299 opiates
            coc_crack_during = sum(devhx_9_coc_crack == 1),                              # n = 316 crack/cocaine
            birthweight_miss = sum(is.na(birthweight_gr)),                               # n = 509 have birth weight missing -- we must be sure that they are had not extreme low birth weight
            extreme_low_birthweight = sum(birthweight_gr < 1200, na.rm = T),             # n = 22 had extreme low birth weight not caught by the screener (they answer no)
            extreme_premature = sum(devhx_12_p > 12),                                    # n = 43 were extremely premature (born before 28 week of gestation); not caught by the screener
            premature_yes_wk_miss = sum(devhx_12a_p == 1 & is.na(devhx_12_p)))           # n = 0 were premature and had missing data for weeks premature

# Now extract the subject ids of those who met this exclusion criteria (study.excl)

study.excl <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>%
  filter(eventname == 'baseline_year_1_arm_1',
         !src_subject_id %in% abcd.excl) %>%
  mutate(weight_kg1 = birth_weight_lbs * 453.6, 
         weight_kg2 = birth_weight_oz * 28.35,  
         devhx_12_p = ifelse(devhx_12a_p == 0 & is.na(devhx_12_p), 0, ifelse(devhx_12a_p == 0  & devhx_12_p == 999, NA, devhx_12_p)), # if not premature, premature wk 0, if premature but premature wk unknown (999) then NA and further remove
         birthweight_gr = ifelse(!is.na(weight_kg1) & !is.na(weight_kg2), weight_kg1 + weight_kg2, ifelse(!is.na(weight_kg1) & is.na(weight_kg2), weight_kg1, NA))) %>% 
  mutate_at(vars(devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack), ~na_if(., 999)) %>% 
  mutate_at(vars(devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack), ~replace_na(., 1)) %>%
  select(src_subject_id, eventname, devhx_1_p, devhx_9_alcohol, devhx_9_tobacco, devhx_9_marijuana, devhx_9_oxycont, devhx_9_her_morph, devhx_9_coc_crack, birthweight_gr, devhx_12a_p, devhx_12_p) %>% 
  filter(devhx_1_p == 0 |                         
         devhx_9_alcohol == 1 |                   
         devhx_9_tobacco == 1 |                   
         devhx_9_marijuana == 1 |                 
         devhx_9_oxycont == 1 |                   
         devhx_9_her_morph == 1 |                 
         devhx_9_coc_crack == 1 |                 
           birthweight_gr < 1200 |
           devhx_12a_p == 999 |
           devhx_12_p > 12 |
           is.na(birthweight_gr) | 
           is.na(devhx_12_p)) %>%              
  select(src_subject_id) %>% 
  pull()

length(study.excl)        # this is the amount of individuals excluded with study-specific criteria (n = 2,619)

(n.excluded2 <- n.excluded - length(study.excl))    # this is your sample (11,636 - 2,619)

length(do.not.include <- union(study.excl, abcd.excl))    #  join lists of excluded subjects from ABCD (n = 231) and study-specific criteria (n = 2,619) to prevent from enter further steps (n = 2,850)

(n.excluded3 <- original.sample.n - length(do.not.include))   # this is how the sample looks like after applying exclusion criteria 1

#5# Now that we established whose eligible (or more like who isn't) let's load demographic data

df <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_demo.csv')) %>%
  # only retain baseline and 2y observations & remove subjects from the exclusion criteria list
  filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1'),
         !src_subject_id %in% do.not.include) %>%
  mutate(child_sex = factor(ifelse(demo_sex_v2 == 1, 'Male', ifelse(demo_sex_v2 == 2, 'Female', NA))),            # we are setting intersex participants to NA
         education_parent1 = ifelse(demo_prnt_ed_v2 == 777 | demo_prnt_ed_v2 == 999, NA, demo_prnt_ed_v2),        # set 999 (do not know/do not remember) and 777 (refuse to answer) to NA
         education_parent2 = ifelse(demo_prtnr_ed_v2 == 777 | demo_prtnr_ed_v2 == 999, NA, demo_prtnr_ed_v2),  
         # find the maximum value of education IN YEARS and assign that value to the household
         education1 = ifelse(!is.na(education_parent1) & is.na(education_parent2), education_parent1,                  # if parent 1 education is NOT NA but parent 2 is, keep parent 1 education level
                                ifelse(is.na(education_parent1) & !is.na(education_parent2), education_parent2,        # if parent 1 education is NA but parent 2 is not, keep parent 2 education level
                                       ifelse(!is.na(education_parent1) & !is.na(education_parent2), pmax(education_parent1, education_parent2), NA))),  # if both have education, keep max, else (both are missing) set to NA
         # Create an ordinary variable to be used in the analysis
         education2 = case_when(education1 < 13 ~ 1,               # less than high school '>HS'
                                between(education1, 13, 14) ~ 2,   # high school graduate 'HS/GED'
                                between(education1, 15, 17) ~ 3,   # Some college
                                education1 == 18 ~ 4,              # Bachelor degree
                                between(education1, 19, 21) ~ 5,   # Postgraduate degree
                                TRUE ~ NA),
         # Recode the ordinary variable as a 5-level factor but only for tables
         education3 = factor(case_when(education2 == 1 ~ '>HS',
                                       education2 == 2 ~ 'HS/GED',
                                       education2 == 3 ~ 'Some college',
                                       education2 == 4 ~ 'Bachelor degree',
                                       education2 == 5 ~ 'Postgrad degree',
                                       TRUE ~ NA), levels = c('>HS', 'HS/GED', 'Some college', 'Bachelor degree', 'Postgrad degree'))) %>%    
  # group by subject to map/carry observations from baseline to follow-up (year 2) -- basically repeat these time-invariant values so these are not NA in follow-up
  group_by(src_subject_id) %>% 
  fill(c(child_sex, education1, education2, education3), .direction = 'updown') %>% 
  ungroup() %>%  # important, ungroup so any transformation is no longer done at the src_subject_id level
  select(src_subject_id, eventname, child_sex, education1, education2, education3) %>%    # keep only these variables
  # append child's age at visit (months), site, and family ID and type of visit
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_y_lt.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%    # no longer needed to further exclude for subjects bc inner_join matches the ones available in this already curated data object
               # rename a few things for convenience
               rename(site = site_id_l,
                      child_age = interview_age) %>% 
               # extend rel_family_id to 2y follow-up
               group_by(src_subject_id) %>% 
               fill(rel_family_id, .direction ='updown') %>% 
               ungroup() %>% 
               # recode visit type for interpretation
               mutate(visit_type = as.factor(case_when(visit_type == 1 ~ 'In-person',
                                                       visit_type == 2 ~ 'Remote',
                                                       visit_type == 3 ~ 'Hybrid',
                                                       TRUE ~ NA))) %>% 
               select(src_subject_id, eventname, visit_type, site, rel_family_id, child_age),    # keep these variables only
             by = c('src_subject_id', 'eventname')) %>%  # join this dataset to the previous object by subject and eventname 
  # append handedness 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_ehis.csv')) %>% 
               # handedness should not change so only keep baseline + the other timepoint is 4y follow-up which is incomplete
               filter(eventname == 'baseline_year_1_arm_1') %>% 
               mutate(handedness_raw = rowMeans(select(., c('ehi1b', 'ehi2b', 'ehi3b', 'ehi4b')), na.rm = T),   # continous variable designating handedness from -100 (pure left) to 100 (pure right)
                      handedness = factor(case_when(handedness_raw > 60 ~ 'Right',                              # recode handedness factor variable for interpretation
                                                    handedness_raw < -60 ~ 'Left', 
                                                    between(handedness_raw, -60, 60) ~ 'Mixed',
                                                    TRUE ~ NA))) %>% 
              select(src_subject_id, handedness), 
             # no need to append by eventname because handedness is time invariant; joining by src_subject_id will populate both eventnames regardless
             by = 'src_subject_id') %>%    
  # append other baselien data that is important and should not be NA
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>%
               # only consider baseline data
               filter(eventname == 'baseline_year_1_arm_1') %>% 
                mutate(weight_gr1 = birth_weight_lbs * 453.6, # transform pounds to grams
                       weight_gr2 = birth_weight_oz * 28.35,  # transform ounces to grams
                       birthweight_gr = ifelse(!is.na(weight_gr1) & !is.na(weight_gr2), weight_gr1 + weight_gr2, ifelse(!is.na(weight_gr1) & is.na(weight_gr2), weight_gr1, NA)),  # if grams not missing, add. If grams 1 only, grams 1 only. Ifelse --> NA
                       premature = factor(ifelse(devhx_12a_p == 0, 'No', 'Yes')),
                       premature_wk = ifelse(premature == 'No' & is.na(devhx_12_p), 0, ifelse(premature == 'Yes' & devhx_12_p == 999, NA, devhx_12_p)),    # if premature = 'No' then premature_wk = 0, if premature = 'Yes' and premature_wk is not known = NA
                       # breastfeeding frequency (bfq) in months asked at baseline with open-ended question
                       bfq_months_y0 = ifelse(devhx_18_p > 72, NA, devhx_18_p),            # remove imposible values (182, 184 months of lactation --> likely 82 and 84 but who knows)
                       # recode open ended question as ordinary variable that is at year 3
                       bfq_duration_y0 = case_when(bfq_months_y0 == 0 ~ 0,                 # 0 months
                                                   # 1 at year 3 is reserved to 'Several days' which has no direct translation here
                                                  between(bfq_months_y0, 1, 3) ~ 2,        # 1 t o 3 months
                                                  between(bfq_months_y0, 4, 6) ~ 3,        # 4 t o 6 months
                                                  between(bfq_months_y0, 7, 9) ~ 4,        # 7 t o 9 months
                                                  between(bfq_months_y0, 10, 12) ~ 5,        # 10 t o 12 months
                                                  between(bfq_months_y0, 13, 18) ~ 6,       # 13 to 18 months
                                                  between(bfq_months_y0, 19, 24) ~ 7,        # 19 to 24 months
                                                  bfq_months_y0 > 24 ~ 8,                  # > 24 months
                                                  TRUE ~ NA)) %>% 
                select(src_subject_id, bfq_months_y0, bfq_duration_y0, birthweight_gr, premature, premature_wk),
             by = 'src_subject_id') %>% 
  # make sure to make baselien the reference value in eventname
  mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
  # drop levels of factor variables that are not populated (= 0)
  droplevels()

# Function to check how many missing per eventname in df variables are

df.split <- split(df, df$eventname)   # separate the dataframe into two based on eventname (dataframe 1 = baseline, dataframe 2 = follow-up)

# Create the function

count_missing <- function(data) {
  sapply(data, function(x) sum(is.na(x)))
}

(missing_counts_by_event <- lapply(df.split, count_missing))        # you shouldn't have missing data for past screening/exclusion (i.e., premature_wk, birthweight_gr)
 
(before.na.rm <- df %>%                                             # before removing NAs we had 9,015 at baseline and 8,344 at follow-up
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) %>% 
  select(n) %>% 
  pull())

# Drop NAs -- education and breastfeeding frequency are important variables in our study and cannot be missing

df <- df %>% 
  drop_na()

(after.na.rm <- df %>%                                             # after removing NAs we have 8,873 at baseline and 8,222 at follow-up
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) %>% 
  select(n) %>% 
  pull())

before.na.rm;after.na.rm;before.na.rm[1:2]-after.na.rm[1:2]

# Count number of subjects

df %>% 
  group_by(eventname) %>% 
  count()

# Append breastfeeding frequency assessments at year 3 for cross-reference with baseline --> At year 3, a 8-item Likert question was used instead of open

df2 <- df %>% 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_bfq.csv')) %>%
               # recode 8-likert (0-8) breastfeeding question at year 3 so it matches y0 report (bfq_duration_y0, also 0-8)
               mutate(bfq_breastfeed_time_p = ifelse(bfq_breastfeed_time_p == 999, NA, bfq_breastfeed_time_p)) %>% 
               rename(bfq_duration_y3 = bfq_breastfeed_time_p) %>% 
               select(src_subject_id, bfq_breastfeed_p, bfq_duration_y3), by = 'src_subject_id')

filtered_df <- df2 %>%
  # if breastfeeding at y0 is not NA and is NA at y3, replace with y0
  mutate(bfq_duration_y3 = ifelse(!is.na(bfq_duration_y0) & is.na(bfq_duration_y3), bfq_duration_y0, bfq_duration_y3)) %>% 
  # if breastfeeding at y0 and y3 is reported 0 then keep; if breastfeeding at y0 and y3 is different than 0 and the difference is greater than 1, lose. 
  mutate(keep = ifelse(bfq_duration_y0 == 0 & bfq_duration_y3 == 0, 'keep',
                       ifelse(bfq_duration_y0 > 0 & bfq_duration_y3 > 0 & abs(bfq_duration_y0 - bfq_duration_y3) <= 1, 'keep', 'lose')))    # 1 equals a ±3 month difference, only applicable if there's some breastfeeding reported in both

df3 <- filtered_df %>% 
  filter(keep == 'keep') %>% 
  # we are interested in the first report of breastfeeding at y0, year 3 was just to ensure accuracy
  mutate(bfq_duration = case_when(bfq_duration_y0 == 0 ~ 0,                # 0 means no breastfeeding
                                  between(bfq_duration_y0, 2, 3) ~ 1,      # 1 means 1-6 months breastfeeding (several days, 1-3 and 4-6 months according to year 3)
                                  between(bfq_duration_y0, 4, 5) ~ 2,      # 2 means 7-12 months (7-9 and 10-12 months according to year 3)
                                  bfq_duration_y0 > 5 ~ 3)) %>%            # 3 means > 12 months (13-18, 19-24 and >24 according to year 3) 
  select(-c(bfq_duration_y0, bfq_duration_y3, bfq_breastfeed_p))
  
df3 %>% 
  group_by(eventname, bfq_duration) %>% 
  summarize(n = length(src_subject_id))  

df3 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))  

#4# Load and append MRI data -- note: more exclusion criteria applied (bad scans, incidental findings)

bad.mrif.y0 <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_qc_clfind.csv')) %>%
            filter(eventname == 'baseline_year_1_arm_1',
                   !between(mrif_score, 1, 2)) %>%                                     # flag those with artifacts (0), abnormalities considering referral (3) and immediate referral (4) --> 581 at baseline
            select(src_subject_id) %>% 
            pull()

bad.mrif.y2 <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_qc_clfind.csv')) %>%
  filter(eventname == '2_year_follow_up_y_arm_1',
         !between(mrif_score, 1, 2)) %>%                                              # flag those with artifacts (0), abnormalities considering referral (3) and immediate referral (4) --> 776 excluded at follow up 
  select(src_subject_id) %>% 
  pull()

length(bad.mrif.y0)
length(bad.mrif.y2)

length(bad.mrif <- union(bad.mrif.y0, bad.mrif.y2))          # 664 participants with bad data at baseline and/or follow-up

# T1w and T2w passed ALL criteria to be included (see https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html)
bad.images.y0 <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_qc_incl.csv')) %>%
  filter(eventname == 'baseline_year_1_arm_1',
         imgincl_t2w_include == 0) %>%                                           # flag those with bad t2w --> 1,217 at baseline
  select(src_subject_id) %>% 
  pull()

bad.images.y2 <- fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_qc_incl.csv')) %>%
  filter(eventname == '2_year_follow_up_y_arm_1',
           imgincl_t2w_include == 0) %>%                                          # flag those with bad t2w --> 600 at baseline
  select(src_subject_id) %>% 
  pull()

length(bad.images.y0)
length(bad.images.y2)

length(bad.images <- union(bad.images.y0, bad.images.y2))       # 1,649 participants with bad data at baseline and/or follow-up

df4 <- df3 %>% 
  # append mri id
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_adm_info.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1'),
                      # make sure to not include subjects with bad images (bad,.images, bad.mrif)
                      !src_subject_id %in% bad.images, 
                      !src_subject_id %in% bad.mrif) %>%
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1'),
                      mri_id = factor(mri_info_deviceserialnumber),
                      mri_date = ymd(mri_info_studydate)) %>% 
               select(src_subject_id, eventname, mri_id, mri_date), by = c('src_subject_id', 'eventname')) %>% 
  # append freesurfer data (CT data, desterieux)
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_thk_dst.csv')) %>%      #mri_y_smr_area_dst.csv mri_y_smr_thk_dst.csv mri_y_smr_vol_dst.csv
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               rename_with(.fn = ~paste0("mrisdp_", seq_along(.), '_cth'), .cols = starts_with("mrisdp_")) %>%    
               select(src_subject_id, eventname, !ends_with(c('149_cth', '150_cth', '151_cth'))),
             by = c('src_subject_id', 'eventname')) %>%
  # append freesurfer data (AREA data, desterieux)
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_area_dst.csv')) %>%      #mri_y_smr_area_dst.csv mri_y_smr_thk_dst.csv mri_y_smr_vol_dst.csv
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               rename_with(.fn = ~paste0("mrisdp_", seq_along(.), '_area'), .cols = starts_with("mrisdp_")) %>%    # mrisdp or dmri
               select(src_subject_id, eventname, !ends_with(c('149_area', '150_area', '151_area'))),
             by = c('src_subject_id', 'eventname')) %>%
  # append ICV 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_vol_aseg.csv')) %>%
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1'),
                      smri_vol_scs_intracranialv = smri_vol_scs_intracranialv / 1e6) %>%        # rescale ICV to liters instead of cm3
               select(src_subject_id, eventname, smri_vol_scs_intracranialv), by = c('src_subject_id', 'eventname')) %>% 
  # append intracortical myelin data (IMC data, destrieux)
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_t1_gray_dst.csv')) %>%      # mri_y_smr_t1_contr_dst.csv --> higher t1/t2 ratio fluid (t2) < fat (t1); mri_y_smr_t1_gray_dst.csv or mri_y_smr_t1_white_dst.csv
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               rename_with(.fn = ~paste0("mrisdp_", seq_along(.), '_myelin1'), .cols = starts_with("mrisdp_")) %>%    # mrisdp or dmri
               select(src_subject_id, eventname, !ends_with(c('149_myelin1', '150_myelin1', '151_myelin1'))), 
             by = c('src_subject_id', 'eventname')) %>% 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/imaging/mri_y_smr_t2_gray_dst.csv')) %>%      # mri_y_smr_t1_contr_dst.csv --> higher t1/t2 ratio fluid (t2) < fat (t1); mri_y_smr_t1_gray_dst.csv or mri_y_smr_t1_white_dst.csv
               mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
               rename_with(.fn = ~paste0("mrisdp_", seq_along(.), '_myelin2'), .cols = starts_with("mrisdp_")) %>%    # mrisdp or dmri
               select(src_subject_id, eventname, !ends_with(c('149_myelin2', '150_myelin2', '151_myelin2'))), 
             by = c('src_subject_id', 'eventname')) %>% 
  mutate(across(ends_with("myelin1"), .names = "{.col}_ratio",
                ~ . / get(sub("myelin1$", "myelin2", cur_column())))) %>% 
  select(-c(ends_with(c('myelin1', 'myelin2')))) %>%
  rename_at(vars(ends_with('myelin1_ratio')), ~str_replace(., 'myelin1_ratio', 't1t2_ratio')) %>% 
  drop_na() 

df4 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))  

# check out-of-range values in CTH (normally goes from 1 to 4.5 mm)
(cth.features <- names(df4 %>% select(ends_with('cth'))))

check <- data.frame(matrix(NA, nrow = length(cth.features), ncol = 3))

for (i in 1:length(cth.features)) {
  var <- df4 %>% select(all_of(cth.features[i]))
  check[i,1] <- names(var)
  check[i,2] <- min(var)
  check[i,3] <- max(var)
  colnames(check) <- c('featurenames', 'min', 'max')
}

which(check$max > 4.5)
which(check$min < 1)

# check for negative values in SA and myelin (should not be negative)

(area.features <- names(df4 %>% select(ends_with(c('area', 't1t2_ratio')))))

check <- data.frame(matrix(NA, nrow = length(cth.features), ncol = 3))

for (i in 1:length(area.features)) {
  var <- df4 %>% select(all_of(area.features[i]))
  check[i,1] <- names(var)
  check[i,2] <- min(var)
  check[i,3] <- max(var)
  colnames(check) <- c('featurenames', 'min', 'max')
}

which(check$min < 0)

(names <- check$featurenames[which(check$min < 0)])   # check for negative values

df5 <- df4

# this is your sample after meeting MRI criteria (mrif_score, imgincl_t1w_include, imgincl_t2w_include)
df5 %>% 
  group_by(eventname, bfq_duration) %>% 
  summarize(n = length(src_subject_id)) 

# create a vector of participants with data at y2 but not at baseline (bc of exclusion criteria applied)
length(not.at.baseline <- setdiff(df5$src_subject_id[df5$eventname == '2_year_follow_up_y_arm_1'], 
                                  df5$src_subject_id[df5$eventname == 'baseline_year_1_arm_1'])) 

# remove subjects who had no baseline data
df6 <- df5 %>% 
  filter(!src_subject_id %in% not.at.baseline) %>% 
  drop_na()

df6 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id)) 

#5# Picking one sibling at random

set.seed(91218)       # set seed reproducibility 91218, because one subject of the same family will be randomly picked 

unique_participants <- df6 %>%
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

length(df6$src_subject_id[df6$eventname == 'baseline_year_1_arm_1']) - length(unique_participants)    # number of discarded

length(unique_participants)   # number of unique participants (of 5,101)

# remove random sibling
df7 <- df6 %>%
  filter(src_subject_id %in% unique_participants)

# confirmation
df7 %>%         
  group_by(rel_family_id) %>%
  summarize(count = n()) %>%
  ungroup() %>% 
  filter(count > 2) %>% 
  nrow()                        # this should result in 0; no participants with > 2 rel_family_id (no siblings)

df7 %>%                         # numbers will be the same regardless of seed but want to ensure we pick the same sibling if data is good on both
  group_by(eventname) %>%
  summarize(n = length(src_subject_id))

# Append the rest of demographic information such as ethnicity, race, income, pds (pubertal stage) that is not important if missing or not (just to report in tables)

df8 <- df7 %>% 
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/abcd-general/abcd_p_demo.csv')) %>%
               filter(eventname == 'baseline_year_1_arm_1') %>% 
                  mutate(mixed = rowSums(select(., c(demo_race_a_p___10:demo_race_a_p___25)), na.rm = T),     # if mixed > 0 then the children is mixed race
                         race = factor(case_when(demo_race_a_p___10 == 1 ~ 'White',
                                                 demo_race_a_p___11 == 1 ~ 'Black',
                                                 demo_race_a_p___12 == 1 ~ 'Native',
                                                 demo_race_a_p___13 == 1 ~ 'Alaska',
                                                 demo_race_a_p___14 == 1 ~ 'Hawaian',
                                                 demo_race_a_p___15 == 1 ~ 'Gaumanian',
                                                 demo_race_a_p___16 == 1 ~ 'Samoan',
                                                 demo_race_a_p___17 == 1 ~ 'Other pacific',
                                                 demo_race_a_p___18 == 1 ~ 'Asian indian',
                                                 demo_race_a_p___19 == 1 ~ 'Chinese',
                                                 demo_race_a_p___20 == 1 ~ 'Filipino',
                                                 demo_race_a_p___21 == 1 ~ 'Japanese',
                                                 demo_race_a_p___22 == 1 ~ 'Korean',
                                                 demo_race_a_p___23 == 1 ~ 'Vietnamese',
                                                 demo_race_a_p___24 == 1 ~ 'Other asian',
                                                 demo_race_a_p___25 == 1 ~ 'Other race',
                                                 demo_race_a_p___99 == 1 ~ 'Do not know',
                                                 demo_race_a_p___77 == 1 ~ 'Refuse to answer')),
                        race = factor(case_when(mixed > 1 ~ 'Mixed',
                                                  race %in% c('White', 'Black', 'Do not know', 'Refuse to answer', 'Other race') ~ race,
                                                  race %in% c('Asian indian', 'Chinese', 'Filipino', 'Japanese', 'Korean', 'Other asian') ~ 'Asian',
                                                  race %in% c('Hawaian', 'Gaumanian', 'Samoan', 'Other pacific') ~ 'NHPI',
                                                  race %in% c('Native', 'Alaska') ~ 'American indian, native american',
                                                  TRUE ~ 'Do not answer')),
                        annual_income = factor(case_when(between(demo_comb_income_v2, 1, 6) ~ '<50,000$',
                                                         between(demo_comb_income_v2, 7, 8) ~ '100,000$',
                                                         between(demo_comb_income_v2, 9, 10) ~ '>100,000$',
                                                         demo_comb_income_v2 == 999 ~ 'Do not know',
                                                         demo_comb_income_v2 == 777 ~ 'Refuse to answer')),
                      ethnicity = factor(case_when(demo_ethn_v2 == 1 ~ 'Hispanic',
                                                   demo_ethn_v2 == 2 ~ 'Non-hispanic',
                                                   demo_ethn_v2 == 999 ~ 'Do not know',
                                                   demo_ethn_v2 == 777 ~ 'Refuse to answer',
                                                   is.na(demo_ethn_v2) ~ 'Missing'))) %>% 
               select(src_subject_id, race, ethnicity, annual_income), by = 'src_subject_id') %>%    
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_dhx.csv')) %>%
                            filter(eventname == 'baseline_year_1_arm_1') %>% 
                            mutate_at(vars(devhx_10a3_p:devhx_10m3_p), ~ifelse(is.na(.), 888, .)) %>% 
                            mutate_at(vars(devhx_14a3_p:devhx_14b3_p), ~ifelse(is.na(.), 888, .)) %>% 
                            mutate(mode_of_delivery = factor(ifelse(devhx_13_3_p == 1, 'Cesarean', ifelse(devhx_13_3_p == 0, 'Vaginal', 'Do not know'))),           # if not 1 or 0, means it's 999 --> set to NA to further exclude
                                   pregnancy_complications = as.factor(ifelse(between(rowSums(select(., c(devhx_10a3_p:devhx_10l3_p))), 1, 12), 1,
                                                                              ifelse(rowSums(select(., c(devhx_10a3_p:devhx_10l3_p))) == 0, 0, 999))),  # if ANY complication during pregnancy, 1; na.rm = T bc if they have NA in a few it doesn't matter
                                   birth_complications = as.factor(ifelse(between(rowSums(select(., c(devhx_14a3_p:devhx_14g3_p))), 1, 7), 1, 
                                                                          ifelse(rowSums(select(., c(devhx_14a3_p:devhx_14g3_p))) == 0, 0, 999))),    # if ANY complication during birth, 1; na.rm = T bc if they have NA in a few it doesn't matter
                                   tobacco = as.factor(ifelse(devhx_8_tobacco == 1, 1, ifelse(devhx_8_tobacco == 0, 0, devhx_8_tobacco))),                # use of tobacco BEFORE pregnancy
                                   alcohol = as.factor(ifelse(devhx_8_alcohol == 1, 1, ifelse(devhx_8_alcohol == 0, 0, devhx_8_tobacco))),                # use of alcohol BEFORE  pregnancy
                                   marihuana = as.factor(ifelse(devhx_8_marijuana == 1, 1, ifelse(devhx_8_marijuana == 0, 0, devhx_8_tobacco))),          # use of cannabis BEFORE  pregnancy
                                   cocaine = as.factor(ifelse(devhx_8_coc_crack == 1, 1, ifelse(devhx_8_coc_crack == 0, 0, devhx_8_coc_crack))),          # use of coc/crak BEFORE  pregnancy
                                   opiates = as.factor(ifelse(devhx_8_her_morph == 1, 1, ifelse(devhx_8_her_morph == 0, 0, devhx_8_her_morph))),          # use of her/morph BEFORE  pregnancy
                                   opioids = as.factor(ifelse(devhx_8_oxycont == 1, 1, ifelse(devhx_8_oxycont == 0, 0, devhx_8_oxycont)))) %>% 
               select(src_subject_id, mode_of_delivery, pregnancy_complications, birth_complications, alcohol, tobacco, marihuana, cocaine, opiates, opioids),
             by = 'src_subject_id') %>% 
    # append PDS; take the average of parent and youth report if possible; if not only parent or youth 
    inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_p_pds.csv')) %>%
                 filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
                 left_join(fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/physical-health/ph_y_pds.csv')), 
                           by = c('src_subject_id', 'eventname')) %>%
                 mutate(across(c('pds_bdyhair_y', 'pds_2_p',      # body hair questions (1 = not yet, 2 = barely started, 3 = underway, 4 = complete, 999/777 = NA)
                                 'pds_m4_y', 'pds_m4_p',         # males only, deepening of your voice? (seem 1-4 scale as above)
                                 'pds_m5_y', 'pds_m5_p',         # males only, face hair?
                                 'pds_f4_2_y','pds_f4_p',       # females only, breasts grown?
                                 'pds_f5_y', 'pds_f5b_p'),                    # females only, menstruation? (1 = no, 4 = yes, 777/999 = NA)
                               ~replace(., . %in% c(999, 777), NA)),
                        pds_female_y = pds_bdyhair_y + pds_f4_2_y,
                        pds_female_p = pds_2_p + pds_f4_p,
                        pds_male_y = pds_bdyhair_y + pds_m4_y + pds_m5_y,
                        pds_male_p = pds_2_p + pds_m4_p + pds_m5_p,
                        pds_male_y = case_when(pds_male_y == 3 ~ 1,
                                               between(pds_male_y, 4, 5) ~ 2,
                                               between(pds_male_y, 6, 8) ~ 3,
                                               between(pds_male_y, 9, 11) ~ 4,
                                               pds_male_y > 11 ~ 5,
                                               TRUE ~ NA),
                        pds_male_p = case_when(pds_male_p == 3 ~ 1,
                                               between(pds_male_p, 4, 5) ~ 2,
                                               between(pds_male_p, 6, 8) ~ 3,
                                               between(pds_male_p, 9, 11) ~ 4,
                                               pds_male_p > 11 ~ 5,
                                               TRUE ~ NA),
                        pds_female_y = case_when(pds_female_y == 3 & pds_f5_y == 1 ~ 2,
                                                 pds_female_y == 4 & pds_f5_y == 1 ~ 3,
                                                 between(pds_female_y, 4, 7) & pds_f5_y == 4 ~ 4,
                                                 between(pds_female_y, 4, 7) & pds_f5_y == 1 ~ 4,
                                                 pds_female_y >= 8 & pds_f5_y == 4 ~ 5,
                                                 TRUE ~ NA),
                        pds_female_p = case_when(pds_female_p == 2 & pds_f5b_p == 1 ~ 1,
                                                 pds_female_p == 3 & pds_f5b_p == 1 ~ 2,
                                                 pds_female_p > 3 & pds_f5b_p == 1 ~ 3,
                                                 between(pds_female_p, 4, 7) & pds_f5b_p == 1 ~ 4,
                                                 between(pds_female_p, 4, 7) & pds_f5b_p == 4 ~ 4,
                                                 pds_female_p >= 8 & pds_f5b_p == 4 ~ 5,
                                                 TRUE ~ NA),
                        pds_male = ifelse(!is.na(pds_male_y) & !is.na(pds_male_p), I(pds_male_y + pds_male_p) / 2,
                                          ifelse(is.na(pds_male_y) & !is.na(pds_male_p), pds_male_p, 
                                                 ifelse(!is.na(pds_male_y) & is.na(pds_male_p), pds_male_y, NA))),
                        pds_female = ifelse(!is.na(pds_female_y) & !is.na(pds_female_p), I(pds_female_y + pds_female_p) / 2,
                                            ifelse(is.na(pds_female_y) & !is.na(pds_female_p), pds_female_p,
                                                   ifelse(!is.na(pds_female_y) & is.na(pds_female_p), pds_female_y, NA))),
                        pds = coalesce(pds_male, pds_female),
                        eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) %>% 
                 select(src_subject_id, eventname, pds),
               by = c('src_subject_id', 'eventname')) %>% 
  mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1'))

# check skewness and kurtosis in numeric IV

e1071::skewness(df8$bfq_duration[df8$eventname == 'baseline_year_1_arm_1'])
e1071::kurtosis(df8$bfq_duration[df8$eventname == 'baseline_year_1_arm_1'] - 3)
hist(df8$bfq_duration[df8$eventname == 'baseline_year_1_arm_1'])

e1071::skewness(df8$bfq_duration[df8$eventname != 'baseline_year_1_arm_1'])
e1071::kurtosis(df8$bfq_duration[df8$eventname != 'baseline_year_1_arm_1'] - 3)
hist(df8$bfq_duration[df8$eventname != 'baseline_year_1_arm_1'])

e1071::skewness(df8$education2[df8$eventname == 'baseline_year_1_arm_1'])
e1071::kurtosis(df8$education2[df8$eventname == 'baseline_year_1_arm_1'] - 3)

e1071::skewness(df8$education2[df8$eventname != 'baseline_year_1_arm_1'])
e1071::kurtosis(df8$education2[df8$eventname != 'baseline_year_1_arm_1'] - 3)

e1071::skewness(df8$child_age[df8$eventname == 'baseline_year_1_arm_1'])
e1071::kurtosis(df8$child_age[df8$eventname == 'baseline_year_1_arm_1'] - 3)

e1071::skewness(df8$child_age[df8$eventname != 'baseline_year_1_arm_1'])
e1071::kurtosis(df8$child_age[df8$eventname != 'baseline_year_1_arm_1'] - 3)

write.csv(df8, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_NO_COMBAT.csv'), row.names = F)    # save files BEFORE removing scanner site with combat in case someone wants the raw MRI data

#6# Remove MRI device effects from cortical thickness and t1w/t2w data (Destrieux)

library(longCombat)

(featurenames <- c(names(df8 %>% select(., starts_with('mrisdp'))), 'smri_vol_scs_intracranialv'))

length(featurenames)    # 148 CT + 148 SA + 148 myelin + 1 ICV = 445 MRI features to harmonize

covs2scale <- c('child_age', 'bfq_duration', 'education2')     # combat expects the lme to be THE SAME that is going to be run for the principal analysis, so mean-center and effects code accordingly

# generate tmp dataframe to store combat results

tmp <- df8 %>% 
  mutate_at(vars(all_of(covs2scale)), ~as.numeric(scale(., scale = F))) %>% 
  mutate(child_sex = ifelse(child_sex == 'Male', 1, -1),
         premature = ifelse(premature == 'Yes', 1, -1)) %>%
  select(src_subject_id, mri_id, child_sex, premature, handedness, education2, all_of(c(covs2scale, featurenames)))

contrasts(tmp$handedness) <- contr.sum(levels(tmp$handedness))

# apply longCombat to other features

simdata_combat <- longCombat(idvar = 'src_subject_id', 
                             timevar = 'child_age',        # time variable is not eventname but (centered) age
                             batchvar = 'mri_id',          # batch effect is MRI serial device ID, not site
                             features = featurenames,      # featurenames = MRI features to remove MRI effects from
                             formula = 'child_age * bfq_duration + premature + education2 + child_sex + handedness',   
                             ranef= '(1 | src_subject_id)',   # random intercept for subject bc repeated measurements
                             data = tmp)

# Update df8 with harmonized MRI features from tmp

df8 %>% 
  droplevels() %>% 
  select(., -all_of(featurenames)) %>%       # this will remove the original features from the dataframe 
  cbind(., simdata_combat$data_combat %>%    # this is the dataframe w/features resulting from longCombat 
          rename_at(vars(all_of(paste0(featurenames, '.combat'))), ~str_remove(., '.combat')) %>%    # we need to remove the '.combat' suffix
          select(., all_of(featurenames))) -> df.mri     # update df

# This will check that there are no negative values in the MRI features

count_neg_values <- function(data, featurenames) {
  total_negative <- 0
  for (var_name in featurenames) {
    neg_rows <- data %>% filter(!!sym(var_name) < 0) %>% select(src_subject_id, all_of(var_name))
    # Count negative and positive values for the current variable
    neg_count <- nrow(neg_rows)
    # Accumulate the counts
    total_negative <- total_negative + neg_count
  }
  # Print total counts
  cat("Total number of negative values:", total_negative, "\n")
}

count_neg_values(df8, featurenames)     # is everything ok?

write.csv(df.mri, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_original.csv'), row.names = F)     # save the harmonized combat values

# Check skewness and kurtosis in the MRI features

names(mri.check <- df.mri %>% select(starts_with(c('mrisdp'))))     

before <- as.data.frame(matrix(NA, ncol = 8, nrow = length(mri.check)))    # empty dataframe to store before summaries

names(before) <- c('ROI', 'n', 'min', 'max', 'Sk', 'Sk.se', 'Skw_ratio', 'Kurtosis')    # insert names to columns

for (i in 1:ncol(mri.check)) {
  before[,1] <- names(mri.check)       # insert feature MRI name
  column_data <- mri.check[[i]]        # store tmp data
  before[i,2] <- length(column_data)   # compute n 
  before[i,3] <- min(column_data)      # compute min
  before[i,4] <- max(column_data)      # compute max
  before[i,5] <- summarytools::descr(column_data)[[11]]   # calculate skewness
  before[i,6] <- summarytools::descr(column_data)[[12]]   # calculate skewness standard error
  before[i,7] <- before[i,5] / before[i,6]                # calculate skeness / skeness se ratio
  before[i,8] <- summarytools::descr(column_data)[[13]]   # calculate kurtosis
}

sum(abs(before$Sk) > 1)                    # amount of variables that are significantly skewed --> 6 out of 445
sum((before$Kurtosis - 3) > 3)             # amount of variables with significantly altered kurtosis --> 8 out of 445
sum(abs(before$Sk) > 1 & (before$Kurtosis - 3) > 3)  # amount of variables with both altered skewness and kurtosis --> 2 out of 445

(to.fix <- subset(before$ROI, (abs(before$Sk) > 1) | (before$Kurtosis - 3) > 3))   # these are the variables with altered skewness or kurtosis 

# remove outliers from those variables to see if that helps

df.mri2 <- df.mri %>%
  group_by(eventname) %>%
  mutate(across(all_of(to.fix), ~ifelse(scale(.) > 3 | scale(.) < -3, NA, .))) %>%      
  ungroup() 

mri.check2 <- df.mri2 %>% select(all_of(to.fix))

after <- as.data.frame(matrix(NA, ncol = 8, nrow = length(mri.check2)))                # empty dataframe to test the variables we just removed outliers from

names(after) <- c('ROI', 'n', 'min', 'max', 'Sk', 'Sk.se', 'Skw_ratio', 'Kurtosis')

for (i in 1:ncol(mri.check2)) {
  after[,1] <- names(mri.check2)    # names
  column_data <- mri.check2[[i]]     # store tmp data
  after[i,2] <- sum(!is.na(column_data))   # n 
  after[i,3] <- min(column_data, na.rm = T)   # min
  after[i,4] <- max(column_data, na.rm = T)   # max
  after[i,5] <- summarytools::descr(column_data)[[11]]   # skewness
  after[i,6] <- summarytools::descr(column_data)[[12]]   # skewness standard error
  after[i,7] <- after[i,5] / after[i,6]                        # skeness / skeness se ratio
  after[i,8] <- summarytools::descr(column_data)[[13]]    # kurtosis
}

after                                # is everything ok?
sum(abs(after$Sk) > 1)               # amount of variables that are skewed --> 0 of 12
sum((after$Kurtosis - 3) > 3)        # amount of variables with altered kurtosis --> 0 of 12

# because we drop data from non-normal variables, we now must secure that the NA, if present in baseline, is also present at follow-up; BUT if present only in follow-up, keep baseline

(participants_with_missing_baseline_but_not_followup <- df.mri2 %>%
  filter(eventname == "baseline_year_1_arm_1" & if_any(all_of(to.fix), is.na)) %>%  # NA in baseline
  select(src_subject_id, all_of(to.fix)) %>%
  inner_join(df.mri2 %>% filter(eventname == "2_year_follow_up_y_arm_1" & if_any(all_of(to.fix), ~ !is.na(.))),  # Valid data in follow-up
    by = "src_subject_id") %>%
  pull(src_subject_id))

df.mri2 %>% 
  filter(src_subject_id %in% participants_with_missing_baseline_but_not_followup) %>% 
  select(src_subject_id, eventname, all_of(to.fix)) %>% 
  head()

pre_check <- df.mri2 %>%
  group_by(src_subject_id) %>%
  filter(any(eventname == "baseline_year_1_arm_1" & if_any(all_of(to.fix), is.na)) &
           any(eventname == "2_year_follow_up_y_arm_1" & if_any(all_of(to.fix), is.na))) %>%
  select(src_subject_id, eventname, all_of(to.fix)) %>%
  ungroup()

pre_check

df.mri3 <- df.mri2 %>%
  group_by(src_subject_id) %>%
  mutate(across(all_of(to.fix), ~ ifelse(eventname == "2_year_follow_up_y_arm_1" & 
                                        any(eventname == "baseline_year_1_arm_1" & is.na(.x)), NA, .x))) %>%
  ungroup()

post_check <- df.mri3 %>%
  group_by(src_subject_id) %>%
  filter(any(eventname == "baseline_year_1_arm_1" & if_any(all_of(to.fix), is.na)) &
           any(eventname == "2_year_follow_up_y_arm_1" & if_any(all_of(to.fix), is.na))) %>%
  select(src_subject_id, eventname, all_of(to.fix)) %>%
  ungroup()

post_check

write.csv(df.mri3, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_MRI_data_no_outliers.csv'), row.names = F)     # save the final dataset, with harmonized and normally-distributed data

#~#~#~#~#~#~#~#~#~#~#
# Append cognition ~#
#~#~#~#~#~#~#~#~#~#~#

df9 <- df8 %>% 
  select(., -starts_with('mrisdp')) %>% 
  # append dataset with vision scores for fluid cognition analysis
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_svs.csv')) %>% 
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               rename(vision = snellen_va_y) %>% 
               select(src_subject_id, eventname, vision) %>% 
               drop_na(),
             by = c('src_subject_id', 'eventname')) %>% 
  # append raw crystalized and fluid uncorrected scores from NIH toolbox
  inner_join(., fread(paste0(datadir, 'ABCD_data/ABCD5.1/core/neurocognition/nc_y_nihtb.csv')) %>%
               filter(eventname %in% c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1')) %>%
               select(src_subject_id, eventname, nihtbx_flanker_uncorrected, nihtbx_pattern_uncorrected, nihtbx_picture_uncorrected) %>%
               drop_na(), by = c('src_subject_id', 'eventname')) %>% 
  mutate(eventname = relevel(factor(eventname), ref = 'baseline_year_1_arm_1')) 

df9 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))

df10 <- df9 %>%   #  get rid of remote/hybrid evaluations8
  filter(visit_type == 'In-person')

df10 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))

# check who is in follow-up but not at baseline

length(not.at.baseline <- setdiff(df10$src_subject_id[df10$eventname == '2_year_follow_up_y_arm_1'], 
                                  df10$src_subject_id[df10$eventname == 'baseline_year_1_arm_1']))

# remove those who are not at baseline

df10 <- df10 %>%
  filter(!src_subject_id %in% not.at.baseline)

df10 %>% 
  group_by(eventname) %>% 
  summarize(n = length(src_subject_id))

# harmonize cognition with longcombat

(featurenames2 <- names(df10 %>% select(starts_with('nih'))))

(covs2scale2 <- c('child_age', 'bfq_duration', 'education2', 'vision'))

tmp2 <- df10 %>% 
  mutate_at(vars(all_of(covs2scale2)), ~scale(., scale = F)) %>% 
  mutate(child_sex = ifelse(child_sex == 'Male', 1, -1),
        premature = ifelse(premature == 'Yes', 1, -1)) %>%
  select(src_subject_id, site, child_sex, premature, education2, handedness, all_of(c(covs2scale2, featurenames2)))

simdata_combat <- longCombat(idvar = 'src_subject_id', 
                             timevar = 'child_age',
                             batchvar = 'site', 
                             features = featurenames2, 
                             formula = 'child_age * bfq_duration + premature + education2 + child_sex + vision',
                             ranef= '(1| src_subject_id)',
                             data = tmp2)

df10 %>% 
  droplevels() %>% 
  select(., -matches(featurenames2)) %>%
  cbind(., simdata_combat$data_combat %>% 
          rename_at(vars(matches(featurenames2)), ~str_remove(., '.combat')) %>% 
          select(., matches(featurenames2))) -> df.cog

df.cog <- df.cog %>% 
  group_by(eventname) %>%
  mutate_at(vars(starts_with('nihtbx')), ~scales::rescale(., to = c(0, 19))) %>%
  ungroup() %>%
  mutate(nih_fluid = rowMeans(select(., starts_with('nihtbx')), na.rm = T))

# skewness and kurtosis ok on numeric variables?

e1071::skewness(df.cog$vision[df.cog$eventname == 'baseline_year_1_arm_1'])
hist(df.cog$vision[df.cog$eventname == 'baseline_year_1_arm_1'])
e1071::kurtosis(df.cog$vision[df.cog$eventname == 'baseline_year_1_arm_1'] - 3)

e1071::skewness(df.cog$vision[df.cog$eventname != 'baseline_year_1_arm_1'])
hist(df.cog$vision[df.cog$eventname != 'baseline_year_1_arm_1'])
e1071::kurtosis(df.cog$vision[df.cog$eventname != 'baseline_year_1_arm_1'] - 3)

e1071::skewness(df.cog$nih_fluid[df.cog$eventname == 'baseline_year_1_arm_1'])
hist(df.cog$nih_fluid[df.cog$eventname == 'baseline_year_1_arm_1'])
e1071::kurtosis(df.cog$nih_fluid[df.cog$eventname == 'baseline_year_1_arm_1'] - 3)

e1071::skewness(df.cog$nih_fluid[df.cog$eventname != 'baseline_year_1_arm_1'])
hist(df.cog$nih_fluid[df.cog$eventname != 'baseline_year_1_arm_1'])
e1071::kurtosis(df.cog$nih_fluid[df.cog$eventname != 'baseline_year_1_arm_1'] - 3)
  
write.csv(df.cog, paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_data/ABCD_COG_data.csv'), row.names = F)    # save harmonized cognitive data

