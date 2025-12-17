# clean_up()
library(PoSTutorial)
library(clintrialx)
library(dplyr)
library(httr)
library(jsonlite)
library(purrr)
library(tidyr)

# Fetch data from clinicaltrials.gov----
mmds <- get_studies(add_data = TRUE)

# Using multiple cores, however there is sometimes an issue with mclapply when fetching json data...
# mmds <- get_studies(add_data = FALSE)
# js = parallel::mclapply(mmds$nct, get_json, mc.cores = parallel::detectCores()-4) %>% silence
# all((lapply(js, length) %>% unlist)>1)
#
# all_bl = parallel::mclapply(js, get_baseline,
#                             mc.cores = parallel::detectCores()-4) %>% bind_rows %>% silence
# all_res = parallel::mclapply(js, get_results,
#                              mc.cores = parallel::detectCores()-4) %>% bind_rows %>% silence
# mmds = list(results = all_res, baseline=all_bl)



# Restrict to studies with MMD endpoint at 12 weeks ----
mmds$results = mmds$results %>%
  mutate(include =
           grepl("monthly", tolower(endpoint)) &
           grepl("migraine", tolower(endpoint)) &
           grepl("days", tolower(endpoint)) &
           grepl("day", tolower(unit)) &
           grepl("12", tolower(time_frame)) &
           grepl("difference", tolower(param_type)))      # reported results should be change from baseline
ncts = mmds$results %>% filter(include) %>% distinct(nct) # relevant studies

## Add inclusion logical to baseline table ----
mmds$baseline = mmds$baseline %>%
  left_join(ncts %>% mutate(include=TRUE)) %>%
  mutate(include=ifelse(is.na(include), FALSE, include)) %>%
  silence

# Get non-missing baseline values ----
mmds$baseline %>%
  mutate(param = case_when(grepl("age", param, ignore.case = TRUE)~"age",
                            grepl("height", param, ignore.case = TRUE)~"height",
                            grepl("weight", param, ignore.case = TRUE)~"weight",
                            grepl("body mass index", param, ignore.case = TRUE)~"bmi",
                            grepl("migraine days", param, ignore.case = TRUE)~"md",
                            TRUE ~ param)) %>%
  filter(param %in% c("age","height","weight","bmi","md"),
         grepl("mean", param_type, ignore.case = TRUE),
         include) %>%
  select(-param_type, -limit_type, -lowerLimit, -upperLimit) -> tmp

## Overview of reported baseline characteristics ----
tmp %>%
  distinct(nct, param) %>%
  mutate(value=1) %>%
  pivot_wider(names_from = param,
              values_from = value)
  # Conlusion: only age and mean monthly migraine days is included on clinicaltrials.gov for almost all studies

## Manually lookup values for missing MMD ----
  # NCT01513291 MMD
  # source: https://journals.sagepub.com/doi/10.1177/0333102414544979
  # mean (sd)
  # Placebo   Active
  # 120       115
  # 6.2 (1.9)	6.2 (1.9)

  # NCT02025556 MMD
  # source: https://www.sciencedirect.com/science/article/pii/S1474442215002495
  # mean (sd)
  # Placebo     TEV-48125 225 mg  TEV-48125 675 mg
  # 104	        96                97
  # 11.5 (2.24)	11.5 (1.9)	  11.3 (2.2)

  # NCT03927144 MMD
  # source: https://jamanetwork.com/journals/jamaneurology/fullarticle/2816800
  # mean (sd)
  # Erenumab    OMPMs       Total
  # 413         208         621
  # 9.5 (2.7)   9.1 (2.9)   9.4 (2.8)

  # NCT03855137 MMD
  # source: https://journals.sagepub.com/doi/10.1177/0333102414544979
  # mean (range)
  # Placebo	          Atogepant 30 mg BID	  Atogepant 60 mg QD  Total
  # 246               253                   256                 755
  # 18.00 (8.0, 28.0) 18.00 (8.0, 28.0)   19.00 (9.0, 28.0)     18.2 (8.0, 28.0)



miss = tibble(nct = c("NCT01513291",
                      "NCT02025556",
                      "NCT03927144",
                      "NCT03855137"),
              n = c(120+115,
                    104+96+97,
                    621,
                    755),
              param = "md",
              value = c(6.2,
                        (104*11.5+96*11.5+97*11.3)/(104+96+97), # weighted average means
                        9.4,
                        18.2),
              spread = c(1.9,
                         sqrt((104*2.24^2+96*1.9^2+97*2.2^2)/(104+96+97)), # weighted average of sd's
                         2.8,
                         range2sd(8,28)))

# Add missing mmds to table of baseline characteristics
tmp <- bind_rows(tmp %>%
                   filter(grepl("age|md",param),
                          grepl("total",title, ignore.case = TRUE)) %>%
                   select(nct, n, param, value, spread) %>%
                   mutate(n = as.numeric(n), value=as.numeric(value), spread = as.numeric(spread)),
                 miss)
tmp = tmp  %>%
  group_by(nct, param) %>%
  summarize(n = mean(n),
            mean = mean(value),
            sd = sqrt(mean(spread^2))) %>%
  pivot_wider(names_from = param,
              values_from = c(n, mean,sd)) %>%
  select(-n_age) %>%
  ungroup %>% silence

# Make results table including baseline values ----
mmds$results %>% filter(include) %>%
  left_join(tmp) -> tmp

tmp %>%
  rename(ntotal = n_md,
         group = group_name,
         age_mean = mean_age,
         age_sd = sd_age,
         md_mean = mean_md,
         md_sd = sd_md) %>%
  select(nct, endpoint, group, ntotal, n, starts_with("age"), starts_with("md"), value, lower, upper, p_value) %>%
  filter(!grepl("acute",endpoint, ignore.case = TRUE), # disregard endpoints on acute migraine drugs
         !grepl("mitt",endpoint, ignore.case = TRUE),  # disregard special population
         !grepl("othe",endpoint, ignore.case = TRUE),  # disregard special population
          grepl("placebo",group, ignore.case = TRUE)   # restrict to placebo comparisons
         ) %>%
  mutate(dose = str_extract_all(group, "(?<!\\d)(\\d+)(?=\\s*mg\\b)") %>% as.numeric) %>%
  relocate(dose, .after=group) %>%
  rename(est = value) -> mmd

save(mmd, file = "data/mmd.Rda")
