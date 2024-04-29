#Libraries

library(data.table)
library(tidyverse)
library(lubridate)
library(broom)
library(survival)

#Variable definition

#RODCIS2 = unique personal identifier
#DATPRI = date of admission
#DATUKO = date of discharge
#VEK = age
#POHL = sex
#NBYDL = region of residence
#ZDG = primary diagnosis
#DAUMR = date of death

#Import data for the main analysis

load(file = "path/Comorbidity_SMI/Data/data_main_analysis.RData")

#Establish mortality

data_models <- map(.x = data_main_analysis,
                   ~ .x %>%
                     mutate(mortality_natural_causes = as.numeric(int_overlaps(interval(ymd(DATUKO), ymd("2017-12-31")), interval(DAUMR, DAUMR)) * str_detect(external_cause_of_death, "^V|^X|^Y", negate = TRUE)),
                            mortality_natural_causes = replace(mortality_natural_causes, is.na(mortality_natural_causes), 0),
                            mortality_external_causes = as.numeric(int_overlaps(interval(ymd(DATUKO), ymd("2017-12-31")), interval(DAUMR, DAUMR)) * str_detect(external_cause_of_death, "^V|^X|^Y")),
                            mortality_external_causes = replace(mortality_external_causes, is.na(mortality_external_causes), 0),
                            years_until_death_natural_causes_or_censoring = case_when(mortality_natural_causes == 1 ~ years_until_death,
                                                                                      mortality_external_causes == 1 ~ years_until_death,
                                                                                      TRUE ~ as.numeric(years_diff_mortality_followup))))



#Proportion of natural and unnatural causes of death

imap_dfr(data_models,
         ~ .x %>%
           group_by(exposure) %>%
           summarise(cohort = .y,
                     mortality_natural_causes_prop = sum(mortality_natural_causes)/sum(mortality) * 100,
                     mortality_external_causes_prop = sum(mortality_external_causes)/sum(mortality) * 100) %>%
           ungroup()) %>%
  pivot_wider(names_from = exposure,
              values_from = c(mortality_natural_causes_prop, mortality_external_causes_prop)) %>%
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
  transmute(cohort,
            across(starts_with("mortality"), 
                   ~ formatC(round(., 2), format = "f", digits = 2))) %>%
  select(cohort, 
         ends_with("unexposed"),
         ends_with("exposed")) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/Causes_of_death.csv",
            row.names = FALSE)

#Stratified Cox proportional hazards models
#Fully adjusted

mfull <- imap_dfr(data_models,
                  ~ tidy(coxph(Surv(years_until_death_natural_causes_or_censoring, mortality_natural_causes) ~ exposure + VEK + POHL + discharge_year + strata(strata), 
                               data = .x),
                         conf.int = TRUE,
                         exponentiate = TRUE) %>%
                    mutate(cohort = .y) %>%
                    filter(term == "exposureexposed"))

#Export table with results

mfull %>%
  arrange(cohort) %>%
  transmute(cohort,
            `Natural causes of death` = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                              paste0("(", 
                                                     formatC(round(conf.low, 2), format = "f", digits = 2),
                                                     "; ",
                                                     formatC(round(conf.high, 2), format = "f", digits = 2),
                                                     ")"))) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/HR_sens_analysis_only_natural_causes_of_death.csv",
            row.names = FALSE)