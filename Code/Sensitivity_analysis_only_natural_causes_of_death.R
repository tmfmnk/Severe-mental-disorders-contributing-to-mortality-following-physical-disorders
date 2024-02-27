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

load(file = "path/data_main_analysis.RData")

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
  write.csv(file = "path/HR_sens_analysis_only_natural_causes_of_death.csv",
            row.names = FALSE)