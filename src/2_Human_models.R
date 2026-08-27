# manuscript: Johnson et al. 2026. Large-scale One Health surveillance at 
#             wildlife-human interfaces reveals virus  spillover risk in 
#             the wildlife trade. Nature Microbiology.
# analysis: Fig. 1d Model-adjusted odds ratios for virus detection in people (GLMM)

# Packages ---------------------------------------------------------------------

library(tidyverse)

# Load Data --------------------------------------------------------------------

model_data = read.csv("data2_viruses_humans.csv", stringsAsFactors = F) %>% 
    mutate(age = factor(age, levels = c("18-40y", "2-8y", "9-17y","41y+")))

model_list = list()

# Demographics ------------------------------------------------------------
# Supplementary Table 3.3

demographics = model_data %>% 
    mutate(region = ifelse(region_asia==1, "Asia", "Africa")) %>% 
    select(human_id, region, age, gender) %>% 
    distinct() %>% 
    group_by(region, age, gender) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = c(region, gender), values_from = n)

# Coronaviruses ----------------------------------------------------------------

model_list$cov_mm = lme4::glmer(data = model_data, family = binomial,
                       coronaviruses ~ (1|lab_name) +
                           acute_fever +
                           oral_nasal_swab)

# Influenza viruses -------------------------------------------------------------

model_list$flu_mm = lme4::glmer(data = model_data, family = binomial,
                       influenzas ~ (1|lab_name) +
                           acute_fever +
                           age +
                           oral_nasal_swab)

# Flaviviruses -----------------------------------------------------------------

model_list$flavi_mm = lme4::glmer(data = model_data, family = binomial,
                         flaviviruses ~ (1|lab_name) + 
                             region_asia + 
                             acute_fever +
                             age +
                             blood_serum)

# Paramyxoviruses --------------------------------------------------------------

model_data_pmx = model_data %>% filter(no_livelihood_data==0)

model_list$pmx_mm = lme4::glmer(data = model_data_pmx, family = binomial,
                       paramyxoviruses ~ (1|lab_name) +
                           acute_fever +
                           age +
                           oral_nasal_swab + 
                           healthcare_worker)

# Model Output -----------------------------------------------------------------

summary(model_list$cov_mm)
summary(model_list$flu_mm)
summary(model_list$flavi_mm)
summary(model_list$pmx_mm)

# Odds Ratios ------------------------------------------------------------------
# Fig 1d and Supplementary Table 3.2

ORtable = function(GLMM, model_name) {
    
    round_cols = c("estimate", "conf.low", "conf.high")
    
    t = broom.mixed::tidy(GLMM, conf.int=TRUE, 
                          exponentiate=TRUE, effects="fixed") %>% 
        mutate_at(round_cols, round, 3) %>% 
        mutate(model = model_name) %>% 
        select(model, term, estimate, conf.low, conf.high, p.value)
    
    t
}

model_list$OR = bind_rows(
    ORtable(model_list$cov_mm, "Coronaviruses"),
    ORtable(model_list$flu_mm, "Influenza Viruses"),
    ORtable(model_list$pmx_mm, "Paramyxoviruses"),
    ORtable(model_list$flavi_mm, "Flaviviruses")) %>% 
    filter(term!="(Intercept)")

write.csv(model_list$OR, "human_models_OR.csv", row.names = F)

