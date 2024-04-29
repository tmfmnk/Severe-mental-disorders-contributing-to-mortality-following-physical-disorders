#Libraries

library(data.table)
library(tidyverse)
library(lubridate)
library(collapse)

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
                             "STAV" = "character",
                             "ZAM" = "character",
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
                             "STAV" = "character",
                             "ZAM" = "character",
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

deaths_1994_2013 <- fread(file = "path/UZIS_mortality_raw/zem_1994_2013.csv")
deaths_2014 <- fread(file = "path/UZIS_mortality_raw/zem_2014.csv")
deaths_2015 <- fread(file = "path/UZIS_mortality_raw/zem_2015.csv")
deaths_2016 <- fread(file = "path/UZIS_mortality_raw/zem_2016.csv")
deaths_2017 <- fread(file = "path/UZIS_mortality_raw/zem_2017.csv")

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

#Establishing the cohorts of people with physical health conditions
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

#Retrieving all people with the target physical health condition
#Identifying first occurrence per individual

set.seed(123)
cohorts <- map(.x = dg_codes_and_names,
               ~ hospitalizations_1999_2017 %>%
                 filter(grepl(.x, ZDG)) %>%
                 group_by(RODCIS2, DATPRI, DATUKO) %>%
                 filter(row_number() == sample(1:n(), 1)) %>%
                 group_by(RODCIS2) %>%
                 filter(DATPRI == min(DATPRI)) %>%
                 filter(DATUKO == min(DATUKO)) %>%
                 ungroup())

#Establishing the past history of the target physical health condition

set.seed(123)
past_history <- map(.x = dg_codes_and_names,
                    ~ hospitalizations_1994_2017 %>%
                      filter(grepl(.x, ZDG)) %>%
                      group_by(RODCIS2, DATPRI, DATUKO) %>%
                      filter(row_number() == sample(1:n(), 1)) %>%
                      group_by(RODCIS2) %>%
                      filter(DATPRI == min(DATPRI)) %>%
                      filter(DATUKO == min(DATUKO)) %>%
                      ungroup() %>%
                      select(RODCIS2,
                             DATUKO_historic = DATUKO))

#Excluding incident cases 

cohorts <- map2(.x = cohorts,
                .y = past_history,
                ~ .x %>%
                  anti_join(.x %>%
                              inner_join(.y,
                                         by = "RODCIS2") %>%
                              filter(DATUKO_historic < DATPRI),
                            by = "RODCIS2"))

#Excluding individuals with residence outside of Czechia

cohorts <- map(.x = cohorts,
               ~ .x %>%
                 filter(!grepl("^99", NBYDL))) 

#Establishing the history of severe mental illness
#Identifying first occurrence per individual

set.seed(123)
smi_history <- hospitalizations_1999_2017 %>%
  filter(grepl("^F2|^F31|^F32[2,3]|^F33[2,3]", ZDG)) %>%
  group_by(RODCIS2, DATPRI, DATUKO) %>%
  filter(row_number() == sample(1:n(), 1)) %>%
  group_by(RODCIS2) %>%
  filter(DATPRI == min(DATPRI)) %>%
  filter(DATUKO == min(DATUKO)) %>%
  ungroup() %>%
  select(RODCIS2,
         DATUKO_smi = DATUKO)

#Establishing the exposure status

cohorts <- map(.x = cohorts,
               ~ .x %>%
                 mutate(smi_up_to_5_years_start = ymd(DATPRI) %m-% years(5),
                        smi_up_to_5_years_end =  ymd(DATPRI) - 1) %>%
                 left_join(smi_history,
                           by = "RODCIS2") %>%
                 mutate(exposure = case_when(int_overlaps(interval(smi_up_to_5_years_start, smi_up_to_5_years_end), interval(ymd(DATUKO_smi), ymd(DATUKO_smi))) ~ "exposed",
                                             is.na(DATUKO_smi) ~ "unexposed",
                                             TRUE ~ "unexposed")) %>%
                 select(-DATUKO_smi))

#Matching

all_matched_pairs <- map2(.x = cohorts,
                          .y = cohorts,
                          ~ .x %>%
                            filter(exposure == "exposed") %>%
                            transmute(strata = RODCIS2,
                                      VEK_exposed = VEK,
                                      POHL,
                                      discharge_year = year(ymd(DATUKO))) %>%
                            inner_join(.y %>%
                                         filter(exposure == "unexposed") %>%
                                         transmute(RODCIS2,
                                                   POHL,
                                                   VEK,
                                                   DATPRI,
                                                   DATUKO,
                                                   discharge_year = year(ymd(DATUKO)),
                                                   ZDG,
                                                   NBYDL,
                                                   STAV,
                                                   ZAM,
                                                   exposure),
                                       by = c("POHL", 
                                              "discharge_year")) %>%
                            filter(data.table::between(VEK, VEK_exposed - 3, VEK_exposed + 3)))

#Recursive random sampling
#Define function

recursive_sampling <- function(data, n) {
  setDT(data)
  u <- 1:nrow(data)
  if (anyDuplicated(data, by = c("strata", "RODCIS2"))) {
    u <- u[
      data[,i := .I][
        sample(.N), -i[duplicated(.SD)], .SDcols = c("strata", "RODCIS2")
      ]
    ]
    data[,i := NULL]
  } 
  sampled <- vector(mode(data$RODCIS2), length(u))
  k <- 0L
  data[
    u, {
      i <- which(RODCIS2 %!in% sampled[seq_len(k)])
      i <- i[sample.int(length(i), min(length(i), n))]
      sampled[seq.int(k + 1L, along.with = i)] <- RODCIS2[i]
      k <- k + length(i)
      .SD[i]
    }, strata
  ]
}

#Perform sampling

set.seed(123)
all_matched_sampled_pairs <- map(.x = all_matched_pairs,
                                 ~ .x %>%
                                   group_by(strata_factor = factor(strata, levels = sample(unique(strata)))) %>%
                                   mutate(random_id = cur_group_id()) %>%
                                   ungroup() %>%
                                   arrange(random_id) %>%
                                   recursive_sampling(., 5L) %>%
                                   select(-VEK_exposed,
                                          -random_id,
                                          -strata_factor))

#Sanity check
#Are all unexposed individuals present only once?

all(map_lgl(.x = all_matched_sampled_pairs,
            ~ .x %>%
              group_by(RODCIS2) %>%
              summarise(cond = n() == 1) %>%
              ungroup() %>%
              summarise(cond = all(cond)) %>%
              pull(cond)))

#Combine datasets
#Leaving out unmatched exposed individuals

cohorts_baseline <- map2(.x = all_matched_sampled_pairs,
                         .y = cohorts,
                         ~ .x %>%
                           bind_rows(.y %>% 
                                       filter(exposure == "exposed") %>%
                                       filter(RODCIS2 %in% .x$strata) %>%
                                       mutate(strata = RODCIS2,
                                              discharge_year = year(ymd(DATUKO)))))

#Sanity check
#Are only matched individuals present?

all(map2_lgl(.x = cohorts_baseline,
             .y = all_matched_pairs,
             ~ all(.x$strata %in% .y$strata)))

#Number of unmatched individuals per cohorts

imap_dfr(map2(.x = cohorts,
              .y = cohorts_baseline,
              ~ .x %>%
                filter(exposure == "exposed") %>%
                summarise(n_unmatched = sum(!RODCIS2 %in% .y$strata),
                          prop_unmatched = n_unmatched/n() * 100)),
         ~ .x %>%
           mutate(cohort = .y) %>%
           select(cohort, everything())) %>%
  arrange(cohort) 

#Establish mortality

data_senstivity_analysis_SMI_up_to_5_years <- map(.x = cohorts_baseline,
                                                  ~ .x %>%
                                                    left_join(deaths_1994_2017,
                                                              by = c("RODCIS2" = "RODCIS2")) %>%
                                                    mutate(mortality = as.numeric(int_overlaps(interval(ymd(DATUKO), ymd("2017-12-31")), interval(DAUMR, DAUMR))),
                                                           mortality = replace(mortality, is.na(mortality), 0),
                                                           years_diff_mortality_followup = as.duration(ymd(DATUKO) %--% ymd("2017-12-31"))/dyears(1),
                                                           years_until_death = as.duration(ymd(DATUKO) %--% DAUMR)/dyears(1),
                                                           years_until_death_or_censoring = ifelse(mortality == 1, years_until_death, years_diff_mortality_followup),
                                                           age_death_or_censoring = ifelse(mortality == 0, 
                                                                                           VEK + years_diff_mortality_followup + 1/365.25,
                                                                                           VEK + years_until_death + 1/365.25)) %>%
                                                    mutate(exposure = factor(exposure, 
                                                                             levels = c("unexposed", "exposed"))))


#Description of cohorts

imap_dfr(.x = data_senstivity_analysis_SMI_up_to_5_years,
         ~ .x %>%
           group_by(exposure) %>%
           summarise(cohort = .y,
                     overall_n = formatC(n(), big.mark = " "),
                     age =  paste(formatC(round(mean(VEK), 2), format = "f", digits = 2),
                                  paste0("(",
                                         formatC(round(sd(VEK), 2), format = "f", digits = 2),
                                         ")")),
                     females = paste(formatC(sum(POHL == 2), big.mark = " "),
                                     paste0("(",
                                            formatC(round(sum(POHL == 2)/n() * 100, 2), format = "f", digits = 2),
                                            ")")),
                     discharge_year = paste(median(discharge_year),
                                            paste0("(",
                                                   paste0(quantile(discharge_year, 0.25, na.rm = TRUE),
                                                          "-",
                                                          quantile(discharge_year, 0.75, na.rm = TRUE)),
                                                   ")"))) %>%
           ungroup()) %>%
  pivot_wider(names_from = "exposure",
              values_from = c("overall_n", "age", "females", "discharge_year")) %>%
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
  select(-order) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/Descriptives_senstivity_analysis_SMI_up_to_5_years.csv",
            row.names = FALSE)

#Export data

save(data_senstivity_analysis_SMI_up_to_5_years, 
     file = "path/Comorbidity_SMI/Data/data_senstivity_analysis_SMI_up_to_5_years.RData")
