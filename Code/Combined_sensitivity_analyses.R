#Libraries

library(tidyverse)

#Import individual results tables

sens_analysis_additional_conditions <- read.csv(file = "path/HR_sens_analysis_additional_conditions.csv",
                                                check.names = FALSE)
sens_analysis_past_hospitalizations <- read.csv(file = "path/HR_sens_analysis_past_hospitalizations.csv",
                                                check.names = FALSE)
sens_analysis_past_hospitalizations_minus_smi <- read.csv(file = "path/HR_sens_analysis_past_hospitalizations_minus_SMI.csv",
                                                          check.names = FALSE)
sens_analysis_past_sud <- read.csv(file = "path/HR_sens_analysis_past_SUD.csv",
                                   check.names = FALSE)
sens_analysis_additional_sociodemo <- read.csv(file = "path/HR_sens_analysis_additional_sociodemo.csv",
                                               check.names = FALSE)
sens_analysis_only_natural_causes_of_death <- read.csv(file = "path/HR_sens_analysis_only_natural_causes_of_death.csv",
                                                       check.names = FALSE)

#Combine results

mget(ls(pattern = "^sens_analysis")) %>%
  reduce(left_join,
         by = "cohort") %>%
  mutate(order = case_when(cohort == "Diseases of the circulatory system" ~ 1,
                           cohort == "Hypertension" ~ 2,
                           cohort == "Ischemic heart disease" ~ 3,
                           cohort == "Atrial fibrillation" ~ 4,
                           cohort == "Heart failure" ~ 5,
                           cohort == "Peripheral artery occlusive disease" ~ 6,
                           cohort == "Stroke" ~ 7,
                           cohort == "Diseases of the endocrine system" ~ 8,
                           cohort == "Diabetes mellitus" ~ 9,
                           cohort == "Thyroid disorder" ~ 10,
                           cohort == "Chronic pulmonary diseases" ~ 11,
                           cohort == "Diseases of the gastrointestinal system" ~ 12,
                           cohort == "Ulcer or chronic gastritis" ~ 13,
                           cohort == "Chronic liver disease" ~ 14,
                           cohort == "Inflammatory bowel disease" ~ 15,
                           cohort == "Diverticular disease of intestine" ~ 16,
                           cohort == "Diseases of the urogenital system" ~ 17,
                           cohort == "Chronic kidney disease" ~ 18,
                           cohort == "Prostate disorders" ~ 19,
                           cohort == "Connective tissue disorders" ~ 20,
                           cohort == "Cancers" ~ 21,
                           cohort == "Diseases of the neurological system" ~ 22,
                           cohort == "Epilepsy" ~ 23,
                           cohort == "Parkinson's disease" ~ 24,
                           cohort == "Multiple sclerosis" ~ 25,
                           cohort == "Infectious and parasitic diseases" ~ 26,
                           cohort == "Tuberculosis" ~ 27,
                           cohort == "Chronic viral hepatitis" ~ 28)) %>%
  arrange(order) %>%
  select(cohort,
         `Additional cond`,
         `Past hosp`,
         `Past hosp minus SMI`,
         `Past SUD`,
         `Additional sociodemo`,
         `Natural causes of death`) %>%
  write.csv(file = "path/HR_combined_sens_analyses.csv",
            row.names = FALSE)
