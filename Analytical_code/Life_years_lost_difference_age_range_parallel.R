#Libraries

library(data.table)
library(tidyverse)
library(lillies)
library(patchwork)
library(furrr)

#Import data with cohorts

load(file = "path/Comorbidity_SMI/Data/data_main_analysis.RData")

#Prepare the data for the procedure

data_lyl_no_smi <- map(.x = data_main_analysis,
                ~ .x %>%
                  filter(exposure == "unexposed") %>%
                  select(mortality,
                         age_death_or_censoring,
                         age_first_hosp_comorbid_condition = VEK) %>%
                  as.data.frame())

data_lyl_smi <- map(.x = data_main_analysis,
                    ~ .x %>%
                      filter(exposure == "exposed") %>%
                      select(mortality,
                             age_death_or_censoring,
                             age_first_hosp_comorbid_condition = VEK) %>%
                      as.data.frame())

#Set up parallelization

future::plan(multisession, workers = 24)

#Estimate LYL 

lyl_smi <- future_map(.x = data_lyl_smi,
                      ~ lyl_ci(lyl_range(data = .x,
                                         t0 = age_first_hosp_comorbid_condition,
                                         t = age_death_or_censoring,
                                         status = mortality,
                                         age_begin = 0,
                                         age_end = pmin(80, max(.x$age_first_hosp_comorbid_condition)),
                                         tau = 81),
                               niter = 10000),
                      .options = furrr_options(seed = 123))

lyl_no_smi <- future_map2(.x = data_lyl_no_smi,
                          .y = data_lyl_smi,
                          ~ lyl_ci(lyl_range(data = .x,
                                             t0 = age_first_hosp_comorbid_condition,
                                             t = age_death_or_censoring,
                                             status = mortality,
                                             age_begin = 0,
                                             age_end = pmin(80, max(.y$age_first_hosp_comorbid_condition)),
                                             tau = 81),
                                   niter = 10000),
                          .options = furrr_options(seed = 123))

#Save data

save(lyl_smi, file = "path/Comorbidity_SMI/Data/lyl_smi.RData") 
save(lyl_no_smi, file = "path/Comorbidity_SMI/Data/lyl_no_smi.RData") 

#Load data

load(file = "path/Comorbidity_SMI/Data/lyl_smi.RData")
load(file = "path/Comorbidity_SMI/Data/lyl_no_smi.RData") 

#Differences in life-years lost

lyl_diff_res <- imap_dfr(pmap(list(lyl_smi,
                                   lyl_no_smi,
                                   map(.x = data_lyl_smi,
                                       ~ .x %>%
                                         select(age_first_hosp_comorbid_condition) %>%
                                         pull())),
                              function(x, y, z) {
                                lyl_diff(x, 
                                         y,
                                         weights = z) %>%
                                  as.data.frame()}),
                         ~ .x %>%
                           mutate(cohort = .y))

#Graph

lyl_diff_res %>%
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
  mutate(cohort = fct_rev(factor(cohort, levels = unique(cohort))),
         cohort = fct_relabel(cohort, 
                              ~ ifelse(.x %in% c("Diseases of the circulatory system",
                                                 "Diseases of the endocrine system",
                                                 "Diseases of the gastrointestinal system",
                                                 "Diseases of the urogenital system",
                                                 "Connective tissue disorders",
                                                 "Cancers",
                                                 "Diseases of the neurological system",
                                                 "Infectious and parasitic diseases",
                                                 "Chronic pulmonary diseases"),
                                       .x,
                                       paste0("       ", .x)))) %>%
  ggplot(aes(x = cohort, 
             y = lyl_estimate.TotalLYL)) +
  geom_col(position = position_dodge(width = 0.8), 
           color = "#999999", 
           fill = "#999999",
           width = 0.7) +
  geom_errorbar(aes(ymin = lyl_ci_left.TotalLYL, ymax = lyl_ci_right.TotalLYL),
                position = position_dodge(width = 0.8),
                width = 0.35) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_y_continuous(name = "Difference in life-years lost (95% CI)",
                     limits = c(-10, 15),
                     labels = c(-10, -5, -1, 0, 1, 2, 3, 5, 7, 10, 15),
                     breaks = c(-10, -5, -1, 0, 1, 2, 3, 5, 7, 10, 15)) +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text =element_text(size = 18,
                                 vjust = 0.5,
                                 hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 0),
        axis.text.x = element_text(size = 12))

ggsave(file = "path/Comorbidity_SMI/Results/Diff_life_years_lost_age_range.eps",
       device = "eps",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 300)

#Table

lyl_diff_res %>%
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
            lyl_diff_estimate_with_ci = paste(formatC(round(lyl_estimate.TotalLYL, 2), format = "f", digits = 2),
                                              paste0("(", 
                                                     formatC(round(lyl_ci_left.TotalLYL, 2), format = "f", digits = 2),
                                                     "; ",
                                                     formatC(round(lyl_ci_right.TotalLYL, 2), format = "f", digits = 2),
                                                     ")"))) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/Diff_life_years_lost.csv",
            row.names = FALSE)
