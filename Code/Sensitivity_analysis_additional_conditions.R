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

#Import all hospitalizations from 1994 to 2017

hospitalizations_1994_2015 <-  list.files(path = "path/UZIS_data_raw", 
                                          full.names = TRUE) %>%
  .[str_detect(., "1[6-7].csv$", negate = TRUE)] %>%
  map_dfr(~ fread(.,
                  select = c("RODCIS2" = "character", 
                             "DATPRI" = "integer", 
                             "DATUKO" = "integer", 
                             "VEK"= "integer", 
                             "POHL" = "integer",
                             "NBYDL"= "character",
                             "ZDG" = "character"),
                  header = TRUE,
                  sep = ",",
                  dec = ".",
                  fill = TRUE,
                  encoding = "Latin-1",
                  nThread = 8))

hospitalizations_2016_2017 <- list.files(path = "path/UZIS_data_raw", 
                                         full.names = TRUE) %>%
  .[str_detect(., "1[6-7].csv$", negate = FALSE)] %>%
  map_dfr(~ fread(.,
                  select = c("RODCIS" = "character", 
                             "DATPRI" = "integer", 
                             "DATUKO" = "integer", 
                             "VEK"= "integer", 
                             "POHL" = "integer",
                             "NBYDL"= "character",
                             "ZDG" = "character"),
                  header = TRUE,
                  sep = ";",
                  dec = ".",
                  fill = TRUE,
                  encoding = "Latin-1",
                  nThread = 8)) %>%
  rename(RODCIS2 = RODCIS)

#Combine all data 

hospitalizations_1994_2017 <- hospitalizations_1994_2015 %>%
  bind_rows(hospitalizations_2016_2017)

#Remove partial data

rm(hospitalizations_1994_2015)
rm(hospitalizations_2016_2017)

#Import data on mortality

deaths_1994_2013 <- fread(file = "path/zem_1994_2013.csv")
deaths_2014 <- fread(file = "path/zem_2014.csv")
deaths_2015 <- fread(file = "path/zem_2015.csv")
deaths_2016 <- fread(file = "path/zem_2016.csv")
deaths_2017 <- fread(file = "path/zem_2017.csv")

#Unifying the format of mortality data

deaths_1994_2017 <- deaths_1994_2013 %>%
  transmute(RODCIS2 = RC,
            DAUMR = dmy(DAUMR),
            cause_of_death = trimws(DGP),
            external_cause_of_death = trimws(DGE)) %>%
  bind_rows(deaths_2014 %>%
              transmute(RODCIS2 = RC,
                        DAUMR = dmy(DAUMR),
                        cause_of_death = trimws(DGP),
                        external_cause_of_death = trimws(DGE)),
            deaths_2015 %>%
              transmute(RODCIS2 = RODCIS2,
                        DAUMR = ymd(DAUMR),
                        cause_of_death = trimws(DGP),
                        external_cause_of_death = trimws(DGE)),
            deaths_2016 %>%
              transmute(RODCIS2 = RCZEMAN2,
                        DAUMR = ymd(paste0(UMROK, str_pad(UMRMM, 2, pad = "0"), UMRDD)),
                        cause_of_death = trimws(DGUMR),
                        external_cause_of_death = trimws(DGUMR2)),
            deaths_2017 %>%
              transmute(RODCIS2 = RCZEMAN2,
                        DAUMR = ymd(paste0(UMROK, str_pad(UMRMM, 2, pad = "0"), UMRDD)),
                        cause_of_death = trimws(DGUMR),
                        external_cause_of_death = trimws(DGUMR2)))

#Remove partial data

rm(deaths_1994_2013)
rm(deaths_2014)
rm(deaths_2015)
rm(deaths_2016)
rm(deaths_2017)

#Excluding records with missing values on key variables
#Excluding records with invalid dates

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
  filter(rowSums(is.na(across(c(RODCIS2, DATPRI, DATUKO, VEK, POHL, NBYDL, ZDG)))) == 0) %>%
  filter(!is.na(ymd(DATPRI)) & !is.na(ymd(DATUKO)))

#Excluding individuals with more than one date of death
#Excluding individuals with hospitalizations after death

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
  anti_join(hospitalizations_1994_2017 %>%
              inner_join(deaths_1994_2017, 
                         by = c("RODCIS2" = "RODCIS2")) %>%
              mutate(DATUKO = ymd(DATUKO)) %>%
              group_by(RODCIS2) %>%
              filter(DAUMR < max(DATUKO) | n_distinct(DAUMR) > 1) %>%
              ungroup(),
            by = c("RODCIS2" = "RODCIS2"))

#Excluding records with invalid date overlaps (discharge date after the admission date of another record)

hospitalizations_1994_2017 <- hospitalizations_1994_2017 %>%
  anti_join(map_dfr(.x = hospitalizations_1994_2017 %>%
                      group_split(split_ID = frank(RODCIS2, ties.method = "dense") %/% 100000),
                    ~ .x %>% 
                      select(RODCIS2,
                             DATPRI_index = DATPRI,
                             DATUKO_index = DATUKO) %>%
                      inner_join(.x %>%
                                   select(RODCIS2,
                                          DATPRI_non_index = DATPRI,
                                          DATUKO_non_index = DATUKO), 
                                 by = c("RODCIS2" = "RODCIS2")) %>%
                      filter(DATPRI_non_index < DATPRI_index & DATUKO_non_index > DATPRI_index) %>%
                      pivot_longer(-RODCIS2,
                                   names_to = c(".value", "type"), 
                                   names_pattern = "([^_]+)_(.*)")),
            by = c("RODCIS2" = "RODCIS2",
                   "DATPRI" = "DATPRI",
                   "DATUKO" = "DATUKO"))

#Keeping hospitalizations that occurred between 1st January 1999 and 31st December 2017

hospitalizations_1999_2017 <- hospitalizations_1994_2017 %>%
  filter(year(ymd(DATPRI)) >= 1999 & year(ymd(DATUKO)) <= 2017) 

#Defining the codes and names of physical health conditions

dg_codes_and_names <- set_names(c("^I1[0-3]|^I15|^I2[0-5]|^I48|^I50|^I7[0-4]|^I6[0-4]|^I69",
                                  "^E1[0-4]|^E0[0-5]|^E06[1-9]|^E07",
                                  "^K221|^K2[5-8]|^K29[3-5]|^B1[6-9]|^K70|^K74|^K766|^I85|^K5[0-1]|^K57",
                                  "^N03|^N11|^N1[8-9]|^N40",
                                  "^M0[5-6]|^M0[8-9]|^M3[0-6]|^D86",
                                  "^C0|^C1|^C2|^C3|^C4[0-3]|^C4[5-9]|^C5|^C6|^C7|^C8|^C9[0-7]",
                                  "^G4[0-1]|^G2[0-2]|^G35",
                                  "^A15|^A16|^A17|^B18|^B2[0-4]",
                                  "^J4[0-7]",
                                  "^I1[0-3]|^I15",
                                  "^I2[0-5]",
                                  "^I48",
                                  "^I50",
                                  "^I7[0-4]",
                                  "^I6[0-4]|^I69",
                                  "^E1[0-4]",
                                  "^E0[0-5]|^E06[1-9]|^E07",
                                  "^K221|^K2[5-8]|^K29[3-5]",
                                  "^B1[6-9]|^K70|^K74|^K766|^I85",
                                  "^K5[0-1]",
                                  "^K57",
                                  "^N03|^N11|^N1[8-9]",
                                  "^N40",
                                  "^G4[0-1]",
                                  "^G2[0-2]",
                                  "^G35",
                                  "^A15|^A16|^A17",
                                  "^B18"),
                                c("Diseases of the circulatory system",
                                  "Diseases of the endocrine system",
                                  "Diseases of the gastrointestinal system",
                                  "Diseases of the urogenital system",
                                  "Connective tissue disorders",
                                  "Cancers",
                                  "Diseases of the neurological system",
                                  "Infectious and parasitic diseases",
                                  "Chronic pulmonary diseases",
                                  "Hypertension",
                                  "Ischemic heart disease",
                                  "Atrial fibrillation",
                                  "Heart failure",
                                  "Peripheral artery occlusive disease",
                                  "Stroke",
                                  "Diabetes mellitus",
                                  "Thyroid disorder",
                                  "Ulcer or chronic gastritis",
                                  "Chronic liver disease",
                                  "Inflammatory bowel disease",
                                  "Diverticular disease of intestine",
                                  "Chronic kidney disease",
                                  "Prostate disorders",
                                  "Epilepsy",
                                  "Parkinson's disease",
                                  "Multiple sclerosis",
                                  "Tuberculosis",
                                  "Chronic viral hepatitis"))

#Import data for the main analysis

load(file = "path/data_main_analysis.RData")

#Other physical health conditions
#Establish occurrences in the past 5 years                           

past_history_all_cond <- map(.x = data_main_analysis,
                             ~ .x %>%
                               inner_join(hospitalizations_1994_2017 %>%
                                            select(RODCIS2,
                                                   DATPRI_historic = DATPRI,
                                                   DATUKO_historic = DATUKO,
                                                   ZDG_historic = ZDG),
                                          by = c("RODCIS2" = "RODCIS2")) %>%
                               filter(int_overlaps(interval(ymd(DATPRI) %m-% years(5), ymd(DATPRI)), 
                                                   interval(ymd(DATPRI_historic), ymd(DATUKO_historic)))))

#Create subsets
#Diagnostic groups

past_history_dg_groups <- past_history_all_cond[which(names(past_history_all_cond) %in% c("Diseases of the circulatory system",
                                                                                          "Diseases of the endocrine system",
                                                                                          "Diseases of the gastrointestinal system",
                                                                                          "Diseases of the urogenital system",
                                                                                          "Connective tissue disorders",
                                                                                          "Cancers",
                                                                                          "Diseases of the neurological system",
                                                                                          "Infectious and parasitic diseases",
                                                                                          "Chronic pulmonary diseases"))]


#Specific diseases

past_history_spec_dis <- past_history_all_cond[which(names(past_history_all_cond) %in% c("Hypertension",
                                                                                         "Ischemic heart disease",
                                                                                         "Atrial fibrillation",
                                                                                         "Heart failure",
                                                                                         "Peripheral artery occlusive disease",
                                                                                         "Stroke",
                                                                                         "Diabetes mellitus",
                                                                                         "Thyroid disorder",
                                                                                         "Ulcer or chronic gastritis",
                                                                                         "Chronic liver disease",
                                                                                         "Inflammatory bowel disease",
                                                                                         "Diverticular disease of intestine",
                                                                                         "Chronic kidney disease",
                                                                                         "Prostate disorders",
                                                                                         "Epilepsy",
                                                                                         "Parkinson's disease",
                                                                                         "Multiple sclerosis",
                                                                                         "Tuberculosis",
                                                                                         "Chronic viral hepatitis"))]


#Exclude the target physical health condition and add the rest as individual variables
#Diagnostic groups

data_models_dg_groups <- imap(past_history_dg_groups,
                              function(data, cohort_names) {
                                data %>%
                                  mutate(map_dfr(ZDG_historic,
                                                 function(ICD_code_historic) {
                                                   set_names(str_detect(ICD_code_historic, setdiff(dg_codes_and_names, dg_codes_and_names[match(cohort_names, names(dg_codes_and_names))])),
                                                             setdiff(dg_codes_and_names, dg_codes_and_names[match(cohort_names, names(dg_codes_and_names))]))
                                                 })) %>%
                                  group_by(across(strata:age_death_or_censoring)) %>%
                                  summarise(across(any_of(dg_codes_and_names[match(names(past_history_dg_groups), names(dg_codes_and_names))]), any)) %>%
                                  ungroup()
                              })

#Specific diseases

data_models_spec_dis <- imap(past_history_spec_dis,
                             function(data, cohort_names) {
                               data %>%
                                 mutate(map_dfr(ZDG_historic,
                                                function(ICD_code_historic) {
                                                  set_names(str_detect(ICD_code_historic, setdiff(dg_codes_and_names, dg_codes_and_names[match(cohort_names, names(dg_codes_and_names))])),
                                                            setdiff(dg_codes_and_names, dg_codes_and_names[match(cohort_names, names(dg_codes_and_names))]))
                                                })) %>%
                                 group_by(across(strata:age_death_or_censoring)) %>%
                                 summarise(across(any_of(dg_codes_and_names[match(names(past_history_spec_dis), names(dg_codes_and_names))]), any)) %>%
                                 ungroup()
                             })


#Models
#Diagnostic groups

mfull_dg_groups <- imap_dfr(.x = data_models_dg_groups,
                            ~ tidy(coxph(as.formula(paste("Surv(years_until_death_or_censoring, mortality) ~ exposure + VEK + POHL + discharge_year",
                                                          paste(paste0("`", names(.x)[(which(names(.x) == "age_death_or_censoring")+1):length(names(.x))], "`"), collapse = "+"),
                                                          "strata(strata)",
                                                          sep = "+")),
                                         data = .x),
                                   conf.int = TRUE,
                                   exponentiate = TRUE)%>%
                              mutate(cohort = .y) %>%
                              filter(term == "exposureexposed"))

#Specific diseases

mfull_spec_dis <- imap_dfr(.x = data_models_spec_dis,
                           ~ tidy(coxph(as.formula(paste("Surv(years_until_death_or_censoring, mortality) ~ exposure + VEK + POHL + discharge_year",
                                                         paste(paste0("`", names(.x)[(which(names(.x) == "age_death_or_censoring")+1):length(names(.x))], "`"), collapse = "+"),
                                                         "strata(strata)",
                                                         sep = "+")),
                                        data = .x),
                                  conf.int = TRUE,
                                  exponentiate = TRUE) %>%
                             mutate(cohort = .y)  %>%
                             filter(term == "exposureexposed"))

#Combine results

mfull <- mfull_dg_groups %>%
  bind_rows(mfull_spec_dis)

#Export table with results

mfull %>%
  arrange(cohort) %>%
  transmute(cohort,
            `Additional cond` = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                      paste0("(", 
                                             formatC(round(conf.low, 2), format = "f", digits = 2),
                                             "; ",
                                             formatC(round(conf.high, 2), format = "f", digits = 2),
                                             ")"))) %>%
  write.csv(file = "path/HR_sens_analysis_additional_conditions.csv",
            row.names = FALSE)