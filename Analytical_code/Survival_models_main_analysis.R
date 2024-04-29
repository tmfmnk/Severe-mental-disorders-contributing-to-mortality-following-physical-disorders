#Libraries

library(tidyverse)
library(broom)
library(survival)
library(data.table)
library(EValue)
library(patchwork)
library(survminer)
library(ggpubr)
library(scales)

#Import data

load(file = "path/Comorbidity_SMI/Data/data_main_analysis.RData")

#Testing the proportionality assumption using Schoenfeld residuals

map(.x = data_main_analysis,
    ~ cox.zph(coxph(Surv(years_until_death_or_censoring, mortality) ~ exposure + VEK + POHL + discharge_year + strata(strata),
                    data = .x)))

#Stratified Cox proportional hazards models
#Fully adjusted

mfull <- imap_dfr(data_main_analysis,
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
            estimate = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                             paste0("(", 
                                    formatC(round(conf.low, 2), format = "f", digits = 2),
                                    "; ",
                                    formatC(round(conf.high, 2), format = "f", digits = 2),
                                    ")"))) %>%
  write.csv(file = "path/Comorbidity_SMI/Results/HR_main_analysis.csv",
            row.names = FALSE)

#Plotting the results from stratified Cox proportional hazards models
#Subplot with HRs

reg_plot <- mfull %>%
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
                                       paste0("       ", .x))),
         color = rep(c("white", "gray95"), length.out = n())) %>%
  {ggplot(data = ., aes(x = estimate, y = cohort, xmin = conf.low, xmax = conf.high)) +
      geom_hline(aes(yintercept = cohort, color = color), size = 7) + 
      geom_pointrange(shape = 22, fill = "black") +
      geom_vline(xintercept = 1, linetype = 3) +
      theme_classic() +
      scale_colour_identity() +
      scale_x_log10(limits = c(0.25, 7), 
                    breaks = c(0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.5, 2, 3, 5, 7)) +
      xlab("aHR (95% CI)") +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank())}

#Subplot with counts
#Prepare data with counts

count_tab <- imap_dfr(.x = data_main_analysis,
                      ~ .x %>%
                        group_by(exposure, cohort = .y) %>%
                        summarise(n = paste0(sum(mortality), "/", n())) %>%
                        ungroup() %>%
                        pivot_wider(names_from = exposure,
                                    values_from = n,
                                    names_prefix = "n_")) %>%
  left_join(mfull %>%
              transmute(cohort,
                        estimate = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                                         paste0("(", 
                                                formatC(round(conf.low, 2), format = "f", digits = 2),
                                                "; ",
                                                formatC(round(conf.high, 2), format = "f", digits = 2),
                                                ")"))),
            by = c("cohort")) %>%
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
                                       paste0("       ", .x))),
         color = rep(c("white", "gray95"), length.out = n())) %>%
  select(-order)

#Prepare subplot 

count_plot <- count_tab %>%
  pivot_longer(-c(cohort, color),
               values_transform = list(value = as.character),
               names_to = "variables",
               values_to = "values") %>%
  mutate(variables = factor(variables, 
                            levels = c("n_unexposed", "n_exposed", "estimate"))) %>%
  ggplot(aes(x = variables, y = cohort, label = values)) +
  geom_hline(aes(yintercept = cohort, color = color), size = 7) +
  geom_text(size = 3) +
  scale_x_discrete(position = "top", 
                   labels = c("Unexposed \nevents/total", "Exposed \nevents/total", "aHR \n(95% CI)")) +
  scale_colour_identity() +
  labs(y = NULL, x = NULL) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(hjust = 0),
        axis.ticks = element_blank(),
        axis.title = element_text(face = "bold"))

#Combine subplots

count_plot + reg_plot + plot_layout(widths = c(5, 10))

ggsave(file = "path/Comorbidity_SMI/Results/HR_plot_main_analysis.eps",
       device = "eps",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 300)

#Survival probability plots
#Define custom plotting function

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 10),
          axis.title.y = element_text(face = "bold", size = 10, angle = 90),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(face = "bold", size = 10),
          axis.ticks.x = element_blank(),    
          axis.ticks.y = element_blank(),    
          legend.text = element_text(face = "bold", size = 10),
          legend.title = element_blank())
}

#Programmatic ploting

imap(data_main_analysis,
     function(data, cohort_names) {
       
       p <- ggsurvplot(surv_fit(Surv(years_until_death_or_censoring, mortality) ~ exposure, 
                                data = data), 
                       conf.int = TRUE, 
                       conf.int.alpha = 0.3,
                       risk.table = TRUE,
                       cumcensor = TRUE,
                       cumevents = TRUE,
                       xlim = c(0, 20),
                       break.x.by = 5,
                       legend.labs = c("People without severe mental illness", "People with severe mental illness"),
                       title = paste0("Survival following ", ifelse(cohort_names == "Parkinson's disease", cohort_names, tolower(cohort_names)), " by severe mental illness status"),
                       xlab = "Time (years) since the health condition",
                       ylab = "Survival probability (95% CI)",
                       palette = c("#0072B2", "#D55E00"),
                       legend.title = "",
                       ggtheme = custom_theme())
       
       p1 = p$plot
       p2 = p$table
       p3 = p$ncensor.plot
       p4 <- p$cumevents
       plots = cowplot::plot_grid(p1, p2, p3, p4, align = "v", ncol = 1, rel_heights = c(4, 1, 1, 1))
       
       ggsave(plot = plots,
              filename = paste0("Survival_probability_plot_", cohort_names, ".eps"),
              path = "path/Comorbidity_SMI/Results/Survival_plots",
              device = cairo_ps,
              width = 10, 
              height = 7, 
              dpi = 300)
     }
)

#Absolute risks

imap_dfr(data_main_analysis,
         ~ .x %>%
           group_by(exposure) %>%
           summarise(cohort = .y,
                     n_deaths = sum(mortality == 1),
                     prop_deaths = sum(mortality)/n() * 100) %>%
           ungroup()) %>%
  transmute(cohort,
            exposure,
            n_prop_deaths =  paste(formatC(n_deaths, big.mark = " "),
                                   paste0("(", 
                                          formatC(round(prop_deaths, 2), format = "f", digits = 2),
                                          ")"))) %>%
  pivot_wider(names_from = "exposure",
              values_from = "n_prop_deaths") %>%
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
  write.csv(file = "path/Comorbidity_SMI/Results/Absolute_risks.csv",
            row.names = FALSE)

#E-values

mfull %>%
  transmute(cohort,
            E_value = pmap_dbl(across(c(estimate, conf.low, conf.high)), 
                               ~ as_tibble(evalues.HR(..1, ..2, ..3, rare = FALSE), 
                                           rownames = "name") %>%
                                 filter(name == "E-values") %>%
                                 pluck("point")),
            E_value = ifelse(data.table::between(1, conf.low, conf.high), NA, E_value),
            E_value = formatC(round(E_value, 2), format = "f", digits = 2)) %>%
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
  write.csv(file = "path/Comorbidity_SMI/Results/E_values.csv",
            row.names = FALSE)