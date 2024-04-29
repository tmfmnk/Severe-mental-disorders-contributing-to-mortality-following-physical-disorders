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
#STAV = marital status
#ZAM = work status
#ZDG = primary diagnosis
#DAUMR = date of death

#Import data for the main analysis

load(file = "path/Comorbidity_SMI/Data/data_main_analysis.RData")

#Keeping records with valid values on work and marital status

data_models <- map(.x = data_main_analysis,
                   ~ .x %>%
                     group_by(strata) %>%
                     filter(all(ZAM %in% as.character(c(0:9)) & STAV %in% as.character(c(0:5)))) %>%
                     ungroup())

#Sanity check
#Was anyone excluded?

all(map2_lgl(.x = data_models,
             .y = data_main_analysis,
             ~ nrow(.x) == nrow(.y)))

#Stratified Cox proportional hazards models
#Fully adjusted

mfull <- imap_dfr(data_models,
                  ~ tidy(coxph(Surv(years_until_death_or_censoring, mortality) ~ exposure + VEK + POHL + as.factor(discharge_year) + as.factor(ZAM) + as.factor(STAV) + strata(strata), 
                               data = .x),
                         conf.int = TRUE,
                         exponentiate = TRUE) %>%
                    mutate(cohort = .y) %>%
                    filter(term == "exposureexposed"))

#Export table with results

mfull %>%
  arrange(cohort) %>%
  transmute(cohort,
            `Additional sociodemo` = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                           paste0("(", 
                                                  formatC(round(conf.low, 2), format = "f", digits = 2),
                                                  "; ",
                                                  formatC(round(conf.high, 2), format = "f", digits = 2),
                                                  ")"))) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/HR_sens_analysis_additional_sociodemo.csv",
            row.names = FALSE)