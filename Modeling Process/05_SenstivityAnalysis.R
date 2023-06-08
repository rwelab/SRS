# Overview:
# We performed sensitivity analysis to check the robustness of results from `04_SeavueComparison.R` by reproducing results after removing each of the following:
# * two trials (PRECISE 1, ENACT) associated with the greatest degree of outcome missing data 
# * participant data associated with missing outcome measurements (week 8 CDAI)
# * participant data from any placebo recipient who participated in a trial of ustekinumab from the dataset used to train our placebo model. (data leakage)

# Outputs:
# * Table 3: Sensitivity Analysis


library(tidyverse)  # For general data processing
library(lme4)       # Mixed effect modeling
library(lmerTest)   # Mixed effect model p-values

################################################################################

## Remission Rate Function

RemissionRateSummaryTable <- function(crohnsData, covariate_list){
  
  ## Placebo Model
  placebo_df <- crohnsData %>% filter(TRTGRP=='Placebo')
  f2 <- paste0(paste('CDAI_REDUCTION', paste(covariate_list, collapse='+'), sep='~'), '+(1|TRIAL)')
  mm <- lmerTest::lmer(f2, data=placebo_df)
  
  ## Placebo Imputation
  crohnsData <- crohnsData %>% 
    mutate(PLACEBO_ATTRIBUTABLE = predict(mm, newdata=crohnsData, re.form=~0)) %>% 
    mutate(DRUG_REDUCTION       = CDAI_REDUCTION - PLACEBO_ATTRIBUTABLE)
  
  ## ADA Model
  covariate_list <- covariate_list[-1] # remove YEAR_NORM
  ADA <- crohnsData %>%  filter(DRUG == 'ADA')
  f2 <- paste0(paste('DRUG_REDUCTION', paste(covariate_list, collapse='+'), sep='~'), '+(1|TRIAL)')
  m2_ada <- lmerTest::lmer(f2, data=ADA)
  
  ## UST Naive Active
  UST_active <- crohnsData %>% filter(DRUG == 'UST') %>% 
    filter(TRTGRP=='Active') %>% filter(HxOfTNFi == 0)
  
  UST_summary <- UST_active %>% 
    mutate(REMISSION = if_else(CDAIL_WEEK8 < 150, 1, 0)) %>% 
    summarize(
      DRUG = 'UST Active TNF-Naive',
      N = n(),
      RemCount = sum(REMISSION),
      RemRate  = mean(REMISSION),
      CI_low   = RemRate - qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N),
      CI_high  = RemRate + qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N)
    )
  
  ## Simulated ADA Naive Active (UST Naive Placebo)
  UST_placebo <- crohnsData %>% filter(DRUG == 'UST') %>% 
    filter(TRTGRP=='Placebo') %>% filter(HxOfTNFi == 0)
  
  N_ada = nrow(UST_placebo) # N of simulated ada data 
  ADA_summary <- UST_placebo %>%
    mutate(ADA_WEEK8 = CDAIL_WEEK8 - predict(m2_ada, newdata=UST_placebo, re.form=~0)) %>%
    mutate(REMISSION = if_else(ADA_WEEK8 < 150, 1, 0))%>% 
    summarize(
      DRUG = 'ADA Sim. Active TNF-Naive',
      N = n(),
      RemCount = sum(REMISSION),
      RemRate  = mean(REMISSION),
      CI_low   = RemRate - qt(0.975, df=N_ada)*sqrt((RemRate)*(1-RemRate)/N_ada),
      CI_high  = RemRate + qt(0.975, df=N_ada)*sqrt((RemRate)*(1-RemRate)/N_ada)
    )
  
  ## Output Remission Rate
  rem_summary <- rbind(UST_summary, ADA_summary)
  return(rem_summary)
}

################################################################################

## Missingness Heatmap

crohnsData_wk8 <- read_csv(file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8.csv')

# re-order covariates and trials
levels.x <- c('AGE','SEX','BMI','CDAI_BASELINE',
              'CRP..mgL','HxOfTNFi','STEROID','IMMUNOMOD', 'LOC.ILEAL',
              'CDAI_WEEK8', 'CDAI_REDUCTION')

levels.y <- c('CERTIFI','UNITI1','UNITI2',
              'ENACT','ENCORE',
              'PRECISE1','NCT02499783','CLASSIC','EXTEND')

miss_data <- miss_data %>% 
  mutate(variable = factor(variable, levels=levels.x)) %>%
  mutate(TRIAL = factor(TRIAL, levels=rev(levels.y)))

ggplot(miss_data, aes(x = TRIAL, y = variable, fill = pct_miss)) + 
  geom_tile() + 
  viridis::scale_fill_viridis(name = 'NA%') + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  coord_flip()

################################################################################

## Sensitivity Analysis - Test


crohnsData <- read_csv(file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_imputed.csv')
covariate_list <- c('YEAR_CENT','CDAI_BASELINE_CENT','AGE_CENT','BMI_CENT','CRP..mgL_CENT',
                    'SEX.MALE','HxOfTNFi','STEROID','IMMUNOMOD','LOC.ILEAL')

RemissionRateSummaryTable(crohnsData, covariate_list)

################################################################################

## Sensitivity Analysis - Remove PRECISE1, ENACT


sense1 <- crohnsData %>% 
  filter(!TRIAL %in% c('PRECISE1','ENACT'))

# find n missing
crohnsData %>% filter(TRIAL %in% c('PRECISE1','ENACT')) %>% group_by(TRIAL) %>% count()
nrow(sense1) # removed 1482 participants

RemissionRateSummaryTable(sense1, covariate_list)

# Fisher Exact Test
rem_rates <- data.frame(
  #                        ADA              UST
  'remission'     = c(    60,     67),
  'not-remission' = c(135-60, 149-67)
)
rem_rates

fisher.test(rem_rates)

################################################################################

## Sensitivity Analysis - Removing Missing Week8 CDAI Participants


sense2 <- crohnsData %>% 
  filter(!is.na(CDAI_WEEK8))

# find n missing
crohnsData %>% filter(is.na(CDAI_WEEK8)) %>% group_by(TRIAL) %>% count()
nrow(sense2) # removed 360 participants

RemissionRateSummaryTable(sense2, covariate_list)



# Fisher Exact Test
rem_rates <- data.frame(
  #                        ADA              UST
  'remission'     = c(    65,     66),
  'not-remission' = c(129-65, 148-66)
)
rem_rates

fisher.test(rem_rates)

################################################################################

## Information Leakage - Removing UST trials from placebo attributable model

sense3 <- crohns_data # imported dataset
covariate_list <- c('YEAR_CENT','CDAI_BASELINE_CENT','AGE_CENT','BMI_CENT','CRP..mgL_CENT',
                    'SEX.MALE','HxOfTNFi','STEROID','IMMUNOMOD','LOC.ILEAL')

## Placebo Model - remove UST
placebo_df <- sense3 %>% filter(TRTGRP=='Placebo') %>% filter(DRUG != 'UST')
f2 <- paste0(paste('CDAI_REDUCTION', paste(covariate_list, collapse='+'), sep='~'), '+(1|TRIAL)')
mm <- lmerTest::lmer(f2, data=placebo_df)
summary(mm)

## Placebo Imputation
sense3 <- sense3 %>% 
  mutate(PLACEBO_ATTRIBUTABLE = predict(mm, newdata=sense3, re.form=~0)) %>% 
  mutate(DRUG_REDUCTION       = CDAI_REDUCTION - PLACEBO_ATTRIBUTABLE)

## ADA Model
covariate_list2 <- covariate_list[-1] # remove YEAR_NORM
ADA <- sense3 %>%  filter(DRUG == 'ADA')
f2 <- paste0(paste('DRUG_REDUCTION', paste(covariate_list2, collapse='+'), sep='~'), '+(1|TRIAL)')
m2_ada <- lmerTest::lmer(f2, data=ADA)
summary(m2_ada)

## UST Naive Active
UST_active <- sense3 %>% filter(DRUG == 'UST') %>% 
  filter(TRTGRP=='Active') %>% filter(HxOfTNFi == 0)

UST_summary <- UST_active %>% 
  mutate(REMISSION = if_else(CDAIL_WEEK8 < 150, 1, 0)) %>% 
  summarize(
    DRUG = 'UST Active TNF-Naive',
    N = n(),
    RemCount = sum(REMISSION),
    RemRate  = mean(REMISSION),
    CI_low   = RemRate - qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N),
    CI_high  = RemRate + qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N)
  )

## Simulated ADA Naive Active (UST Naive Placebo)
UST_placebo <- sense3 %>% filter(DRUG == 'UST') %>% 
  filter(TRTGRP=='Placebo') %>% filter(HxOfTNFi == 0)

N_ada = nrow(ADA) # N of ada data used to train ada model
ADA_summary <- UST_placebo %>%
  mutate(ADA_WEEK8 = CDAIL_WEEK8 - predict(m2_ada, newdata=UST_placebo, re.form=~0)) %>%
  mutate(REMISSION = if_else(ADA_WEEK8 < 150, 1, 0))%>% 
  summarize(
    DRUG = 'ADA Sim. Active TNF-Naive',
    N = n(),
    RemCount = sum(REMISSION),
    RemRate  = mean(REMISSION),
    CI_low   = RemRate - qt(0.975, df=N_ada)*sqrt((RemRate)*(1-RemRate)/N_ada),
    CI_high  = RemRate + qt(0.975, df=N_ada)*sqrt((RemRate)*(1-RemRate)/N_ada)
  )

## Output Remission Rate
rem_summary <- rbind(UST_summary, ADA_summary)
rem_summary

################################################################################

## Negative Control

# ADA Negative Control
ADA %>% 
  mutate(REMISSION = if_else(CDAIL_WEEK8 < 150, 1, 0)) %>% 
  summarize(
    DRUG = 'ADA Negative Control',
    N = n(),
    RemCount = sum(REMISSION),
    RemRate  = mean(REMISSION),
    CI_low   = RemRate - qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N),
    CI_high  = RemRate + qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N)
  )

# UST Negative Control
UST %>% 
  filter(HxOfTNFi == 0) %>% 
  filter(TRTGRP == 'Active') %>%
  mutate(REMISSION = if_else(CDAIL_WEEK8 < 150, 1, 0)) %>% 
  summarize(
    DRUG = 'UST Negative Control',
    N = n(),
    RemCount = sum(REMISSION),
    RemRate  = mean(REMISSION),
    CI_low   = RemRate - qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N),
    CI_high  = RemRate + qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N)
  )

# Fisher Exact Test
rem_rates <- data.frame(
  #                       ADA     UST
  'remission'     = c(    119,     67),
  'not-remission' = c(239-119, 149-67)
)
rem_rates

fisher.test(rem_rates)
