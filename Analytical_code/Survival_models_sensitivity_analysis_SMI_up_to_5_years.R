#Libraries

library(tidyverse)
library(broom)
library(survival)
library(data.table)

#Import data

load(file = "path/Comorbidity_SMI/Data/data_senstivity_analysis_SMI_up_to_5_years.RData")

#Testing the proportionality assumption using Schoenfeld residuals

map(.x = data_senstivity_analysis_SMI_up_to_5_years,
    ~ cox.zph(coxph(Surv(years_until_death_or_censoring, mortality) ~ exposure + VEK + POHL + discharge_year + strata(strata),
                    data = .x)))

#Stratified Cox proportional hazards models
#Fully adjusted

mfull <- imap_dfr(data_senstivity_analysis_SMI_up_to_5_years,
                  ~ tidy(coxph(Surv(years_until_death_or_censoring, mortality) ~ exposure + VEK + POHL + discharge_year + strata(strata), 
                               data = .x),
                         conf.int = TRUE,
                         exponentiate = TRUE) %>%
                    mutate(cohort = .y) %>%
                    filter(term == "exposureexposed"))

#Export table with results

mfull %>%
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
            `Severe mental illness 5<= years before the physical health condition` = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                                                                           paste0("(", 
                                                                                                  formatC(round(conf.low, 2), format = "f", digits = 2),
                                                                                                  "; ",
                                                                                                  formatC(round(conf.high, 2), format = "f", digits = 2),
                                                                                                  ")"))) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/HR_sens_analysis_SMI_up_to_5_years.csv",
            row.names = FALSE)
