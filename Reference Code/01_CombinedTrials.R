# Overview:
# Steps on how individual trials were merged together(IndividualTrials folder). Includes which variables were selected, how they were refactored, missingness per variable, and how they were centered.  
# 
# Outputs: 
# * crohnsData_MASTER
# * crohnsData_wk8
# * Missingness Heatmap
# * Table 1 Baseline Characteristics
# * Tables of N and pct missing by trial

library(tidyverse)
library(ggplot2)
# install.packages('naniar')
library(naniar)
library(sjPlot)
library(table1)
library(data.table)

#################################################################################################

# ## Individual files -> Combined file (MASTER)
# Re-factor columns to either <chr> or <dbl>

# get a list of all csv files from dir (individual trials)
dir <- 'G:/Shan/Week 8 Identification/IndividualTrials'
files <- list.files(path=dir, pattern='.csv', full.names=TRUE)
length(files)

# merge all csvs to crohnsData_MASTER
crohnsData_MASTER <- data.frame()
for (file in files) {
  trial <-read_csv(file=file) %>%
    # as.character
    mutate_at(.vars=c('ID','TRIAL','DRUG','SEX','HxOfTNFi','STEROID','IMMUNOMOD','LOC.ILEAL',
                      'RACE','ETHNIC','LOC.COLON','SMOKING','CURR.FISTULA','PERIANAL','CURR.PRIOR.STRICTURE'), 
              ~ as.character(.x)) %>% 
    # as.double
    mutate_at(.vars=c('AGE','BMI','CRP..mgL','HEIGHT..cm','WEIGHT..kg','DURATION','ALBUMIN..gL'), ~ as.double(.x)) %>% 
    mutate_at(vars(contains('CDAI')), ~ as.double(.x))
  
  crohnsData_MASTER <- rbind(crohnsData_MASTER, trial)
  
  rm(trial)
}

#################################################################################################

# Impute NA -> No for categorical vars (IMMUNOMOD, STEROID, HxOfTNFi, LOC.ILEAL)
crohnsData_MASTER <- crohnsData_MASTER %>%
  # categorical variables - change NA to No
  mutate(across(c(IMMUNOMOD, STEROID, HxOfTNFi, LOC.ILEAL, LOC.COLON, PERIANAL, CURR.FISTULA, CURR.PRIOR.STRICTURE), ~ replace_na(.x, 'No'))) %>% 
  # categorical varialbes -> binary
  mutate(IMMUNOMOD    = if_else(IMMUNOMOD=='Yes',1,0)) %>% 
  mutate(STEROID      = if_else(STEROID=='Yes',1,0)) %>%
  mutate(HxOfTNFi     = if_else(HxOfTNFi=='Yes',1,0)) %>%
  mutate(LOC.ILEAL    = if_else(LOC.ILEAL=='Yes',1,0)) %>%
  mutate(LOC.COLON    = if_else(LOC.COLON=='Yes',1,0)) %>%
  mutate(SMOKING      = if_else(SMOKING=='Yes',1,0)) %>% 
  mutate(PERIANAL     = if_else(PERIANAL=='Yes',1,0)) %>%
  mutate(CURR.FISTULA = if_else(CURR.FISTULA=='Yes',1,0)) %>%
  mutate(CURR.PRIOR.STRICTURE = if_else(CURR.PRIOR.STRICTURE=='Yes',1,0)) %>%
  # standardize ethnicity: UNKNOWN, NOT REPORTED, NA -> NA
  mutate(ETHNIC = if_else(tolower(ETHNIC) %like% 'not hispanic', 'NOT HISPANIC OR LATINO', 
                  if_else(tolower(ETHNIC) =='hispanic or latino', 'HISPANIC OR LATINO', NA_character_))) %>% 
  mutate(Sex_Male     = if_else(SEX=='M',1,0)) %>% 
  mutate(Race_White   = if_else(RACE=='WHITE',1,0)) %>% 
  mutate(Ethnicity_Hispanic = if_else(ETHNIC=='HISPANIC OR LATINO',1,0))

write.csv(crohnsData_MASTER,
          file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_MASTER.csv',
          row.names=FALSE)

#################################################################################################

crohnsData_MASTER %>% 
  group_by(TRIAL) %>% 
  dplyr::summarize(N = n())

#################################################################################################

# Missingness for MASTER
crohnsData_MASTER %>% 
  # find those who don't have any cdai measurements btw weeks 1-8
  mutate(WK1_WK8_NA = if_else(is.na(CDAI_WEEK1) & 
                              is.na(CDAI_WEEK2) & 
                              is.na(CDAI_WEEK3) & 
                              is.na(CDAI_WEEK4) &
                              is.na(CDAI_WEEK6) &
                              is.na(CDAI_WEEK8), NA, TRUE)) %>% 
  select(TRIAL, CDAI_BASELINE, WK1_WK8_NA, CDAI_WEEK8, BMI, CRP..mgL, IMMUNOMOD, STEROID, HxOfTNFi, LOC.ILEAL, 
         Race_White, Ethnicity_Hispanic, LOC.COLON, SMOKING, PERIANAL, CURR.FISTULA, CURR.PRIOR.STRICTURE,
         ALBUMIN..gL) %>% 
  group_by(TRIAL) %>%
  summarise(across(everything(), ~sum(is.na(.))))

#################################################################################################

## MASTER -> Week 8 Data

# Isolate patients with no cdai data btw weeks 1-8
id_na <- crohnsData_MASTER %>%
  select(ID, CDAI_WEEK1:CDAI_WEEK8) %>%
  filter_at(vars(-ID), ~is.na(.)) %>%
  select(ID)

crohnsData_wk8 <- crohnsData_MASTER %>%
  # remove patients with no cdai baseline
  filter(!is.na(CDAI_BASELINE)) %>%
  # remove patients with no cdai data btw weeks 1-8
  anti_join(., id_na) %>% 
  # select important columns
  mutate(CDAI_REDUCTION = CDAI_BASELINE - CDAIL_WEEK8, na.rm=T) %>% 
  select(ID, TRIAL, YEAR, DRUG, TRTGRP,  # identifiers
         AGE, SEX, Sex_Male, BMI, CRP..mgL, HxOfTNFi, STEROID, IMMUNOMOD, LOC.ILEAL,  # main covariates
         CDAI_BASELINE:CDAI_WEEK8, # original data 
         CDAIL_WEEK8,              # week 8 LVCF
         CDAI_REDUCTION,           # week 8 cdai reduction
         RACE, ETHNIC, Race_White, Ethnicity_Hispanic, HEIGHT..cm, WEIGHT..kg, LOC.COLON, DURATION, 
         SMOKING, ALBUMIN..gL, CURR.FISTULA, PERIANAL, CURR.PRIOR.STRICTURE) # additional covariates

# N removed = 143 (113 no week1:week8 data, 30 no baseline)
crohnsData_MASTER %>% nrow()
crohnsData_wk8 %>% nrow()

write.csv(crohnsData_wk8,
          file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8.csv',
          row.names=FALSE)

#################################################################################################

# Missingness for wk8
print.data.frame( 
  crohnsData_wk8 %>% 
  select(TRIAL, TRTGRP, BMI, CRP..mgL, CDAI_WEEK8) %>% 
  group_by(TRTGRP, TRIAL) %>% 
  summarise(across(everything(), ~sum(is.na(.)))) 
  )

#################################################################################################

# N Active
dim(crohnsData_wk8 %>% filter(TRTGRP == 'Active'))

# N Placebo
dim(crohnsData_wk8 %>% filter(TRTGRP == 'Placebo'))

crohnsData_wk8 %>% 
  # find those who don't have any cdai measurements btw weeks 1-8
  mutate(WK1_WK8_NA = if_else(is.na(CDAI_WEEK1) & 
                              is.na(CDAI_WEEK2) & 
                              is.na(CDAI_WEEK3) & 
                              is.na(CDAI_WEEK4) &
                              is.na(CDAI_WEEK6) &
                              is.na(CDAI_WEEK8), NA, TRUE)) %>% 
  select(TRTGRP, CDAI_BASELINE, WK1_WK8_NA, CDAI_WEEK8, BMI, CRP..mgL, IMMUNOMOD, STEROID, HxOfTNFi, LOC.ILEAL, 
         Race_White, Ethnicity_Hispanic, LOC.COLON, SMOKING, PERIANAL, CURR.FISTULA, CURR.PRIOR.STRICTURE,
         ALBUMIN..gL) %>% 
  group_by(TRTGRP) %>%
  summarise(across(everything(), ~sum(is.na(.))))

#################################################################################################

## Summary of Categorical Variables

summaryTable <- crohnsData_wk8 %>% 
  select(TRIAL, TRTGRP, IMMUNOMOD, STEROID, HxOfTNFi, LOC.ILEAL) %>% 
  group_by(TRIAL, TRTGRP) %>% 
  summarize_all(mean)

write.csv(summaryTable,
          file='G:/Shan/Week 8 Identification/Tables/CovariateSummaryByStudy.csv',
          row.names=FALSE)

print.data.frame(summaryTable)

#################################################################################################

## Outcomes Summary (Response, Remission, Endpoint)

summaryTable <- crohnsData_wk8 %>%
  mutate(REMISSION = if_else(CDAIL_WEEK8    <= 150,    1, 0)) %>%
  mutate(RESPONSE  = if_else(CDAI_REDUCTION >= 100,    1, 0)) %>%
  mutate(ENDPOINT  = if_else(RESPONSE==1|REMISSION==1, 1, 0)) %>% 
  group_by(TRIAL, TRTGRP, DRUG) %>%
  summarise(N = n(), 
            Response = round(sum(RESPONSE)/n(),3), 
            Remission = round(sum(REMISSION)/n(),3),
            Endpoint = round(sum(ENDPOINT)/n(),3))

write.csv(summaryTable,
          file='G:/Shan/Week 8 Identification/Tables/OutcomeSummaryByStudy.csv',
          row.names=FALSE)

print.data.frame(summaryTable)

#################################################################################################

# Find n, pct missing by trial

# re-order covariates and trials
levels.x <- c('Age','Sex','BMI','Baseline CDAI', 'CRP (mg/L)','HxOfTNFi','Steroid Use','Immunomod Use', 'Ileal Disease', 'Week 8 CDAI')
levels.y <- c('CERTIFI','UNITI1','UNITI2', 'ENACT','ENCORE', 'PRECISE1','NCT02499783','CLASSIC','EXTEND')

miss_data <- crohnsData_wk8 %>%
  select(TRIAL, AGE, SEX, BMI, CDAI_BASELINE, CDAI_WEEK8, CRP..mgL, HxOfTNFi, STEROID, IMMUNOMOD, LOC.ILEAL) %>% 
  rename(Age=AGE) %>%
  rename(Sex=SEX) %>%
  rename('Baseline CDAI'=CDAI_BASELINE) %>% 
  rename('CRP (mg/L)'=CRP..mgL) %>%
  rename('Steroid Use'=STEROID) %>% 
  rename('Immunomod Use'=IMMUNOMOD) %>%
  rename('Ileal Disease'=LOC.ILEAL) %>% 
  rename('Week 8 CDAI'=CDAI_WEEK8) %>% 
  group_by(TRIAL) %>% naniar::miss_var_summary() %>% ungroup() %>% 
  mutate(variable = factor(variable, levels=levels.x)) %>%
  mutate(TRIAL = factor(TRIAL, levels=rev(levels.y)))

p1 <- ggplot(miss_data, aes(x = TRIAL, y = variable, fill = pct_miss)) + 
  geom_tile() + 
  viridis::scale_fill_viridis(limits=c(0,50), name = 'NA%') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  coord_flip()

ggsave(filename = "G:/Shan/Week 8 Identification/Figures/MissingnessSummaryByStudy.pdf",
       plot = p1, height = 5, width = 7)
p1

#################################################################################################

## Missingness n, pct tables


pctMissing <- miss_data %>% 
  select(-n_miss) %>% 
  pivot_wider(names_from='variable', values_from='pct_miss')

write_csv(pctMissing,
          file = "G:/Shan/Week 8 Identification/Tables/PctMissingSummaryByStudy.csv")

nMissing <- miss_data %>% 
  select(-pct_miss) %>% 
  pivot_wider(names_from='variable', values_from='n_miss')

write_csv(nMissing,
          file = "G:/Shan/Week 8 Identification/Tables/nMissingSummaryByStudy.csv")

#################################################################################################

## LVCF Evaluation

crohnsData_wk8 %>% 
  select(ID, TRIAL, TRTGRP, CDAI_WEEK1:CDAI_WEEK8, CDAIL_WEEK8) %>% 
  filter(is.na(CDAI_WEEK8)) %>%  # 360 missing week 8
  filter(!is.na(CDAI_WEEK6)) %>% # 119 LVCF from week 6
  nrow()

crohnsData_wk8 %>% 
  select(ID, TRIAL, TRTGRP, CDAI_WEEK1:CDAI_WEEK8, CDAIL_WEEK8) %>% 
  filter(is.na(CDAI_WEEK8)) %>%  # 360 missing week 8
  filter(is.na(CDAI_WEEK6)) %>%  # 119 LVCF from week 6
  filter(!is.na(CDAI_WEEK4)) %>% # 136 LVCF from week 4
  nrow()

crohnsData_wk8 %>% 
  select(ID, TRIAL, TRTGRP, CDAI_WEEK1:CDAI_WEEK8, CDAIL_WEEK8) %>% 
  filter(is.na(CDAI_WEEK8)) %>%  # 360 missing week 8
  filter(is.na(CDAI_WEEK6)) %>%  # 119 LVCF from week 6
  filter(is.na(CDAI_WEEK4)) %>%  # 136 LVCF from week 4
  filter(!is.na(CDAI_WEEK3)) %>% # 15  LVCF from week 3
  nrow()

crohnsData_wk8 %>% 
  select(ID, TRIAL, TRTGRP, CDAI_WEEK1:CDAI_WEEK8, CDAIL_WEEK8) %>% 
  filter(is.na(CDAI_WEEK8)) %>%  # 360 missing week 8
  filter(is.na(CDAI_WEEK6)) %>%  # 119 LVCF from week 6
  filter(is.na(CDAI_WEEK4)) %>%  # 136 LVCF from week 4
  filter(is.na(CDAI_WEEK3)) %>%  # 15  LVCF from week 3
  filter(!is.na(CDAI_WEEK2)) %>% # 90  LVCF from week 2
  nrow()

#################################################################################################

# ## Missing value imputation
# * 360  CDAI_WEEK8  = LVCF imputation 
# *  80  CRP..mgL    = median imputation by trial
# *   5  BMI         = median imputation by trial

crohnsData_wk8_imputed <- crohnsData_wk8 %>%
  # median imputation by Trial
  group_by(TRIAL) %>% 
  mutate(BMI = if_else(is.na(BMI), median(BMI, na.rm=T), BMI)) %>%
  mutate(CRP..mgL = if_else(is.na(CRP..mgL), median(CRP..mgL, na.rm=T), CRP..mgL)) %>% 
  ungroup() %>% 
  # convert Yes/No to 1/0
  mutate(across(c(where(is.character), -c(ID, TRIAL, TRTGRP, DRUG, SEX)), ~ recode(., 'Yes'=1, 'No'=0))) %>% 
  # convert M/F to 1/0
  mutate(SEX.MALE = if_else(SEX=='M', 1, 0)) %>% 
  # center continuous variables
  mutate(CRP..mgL_CENT = CRP..mgL-10) %>%
  mutate(BMI_CENT      = BMI-20) %>%
  mutate(AGE_CENT      = AGE-35) %>% 
  mutate(CDAI_BASELINE_CENT = CDAI_BASELINE-300) %>%
  mutate(YEAR_CENT     = YEAR-2000) %>% 
  # re-order variables
  select(ID, TRIAL, DRUG, TRTGRP,                                          # identifiers
         YEAR_CENT, AGE_CENT, BMI_CENT, CRP..mgL_CENT, CDAI_BASELINE_CENT, # centered covariates
         YEAR, AGE, BMI, CRP..mgL, CDAI_BASELINE,                          # original covariates
         SEX, SEX.MALE, HxOfTNFi, STEROID, IMMUNOMOD, LOC.ILEAL,           # categorical covariates
         CDAI_REDUCTION, CDAIL_WEEK8,                                      # Y, week 8 reference
         CDAI_WEEK1:CDAI_WEEK8)                                            # original cdai data

write.csv(crohnsData_wk8_imputed,
          file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_imputed.csv',
          row.names=FALSE)

#################################################################################################

## Characteristics Table 1

tbl1 <- crohnsData_wk8_imputed %>% 
  # re-factor categorical variables
  mutate(SEX.FEMALE = abs(SEX.MALE - 1)) %>% 
  mutate(SEX.FEMALE = factor(SEX.FEMALE, 
                             levels=c(1,0), 
                             labels=c('Female','Male'))) %>% 
  mutate(HxOfTNFi   = factor(HxOfTNFi,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(STEROID    = factor(STEROID,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(IMMUNOMOD  = factor(IMMUNOMOD,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(LOC.ILEAL  = factor(LOC.ILEAL,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  filter(DRUG != 'ADA')

label(tbl1$TRTGRP)         <- 'Treatment Group'
label(tbl1$CRP..mgL)       <- 'CRP'
units(tbl1$CRP..mgL)       <- 'mg/L'
label(tbl1$AGE)            <- 'Age'
label(tbl1$BMI)            <- 'BMI'

label(tbl1$CDAI_BASELINE)  <- 'Baseline CDAI'
label(tbl1$CDAI_REDUCTION) <- 'CDAI Reduction'

label(tbl1$SEX.FEMALE) <- 'Sex:Female'
label(tbl1$HxOfTNFi)   <- 'History of TNF Intolerance'
label(tbl1$STEROID)    <- 'Steroid Use'
label(tbl1$IMMUNOMOD)  <- 'Immunomodulator Use'
label(tbl1$LOC.ILEAL)  <- 'Location:Ileal'

my.render.cont <- function(x){
  with(stats.apply.rounding(stats.default(x), digits=2), 
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

my.render <- function(x,...) {
  y <- render.default(x,...)
  if(is.factor(x)) y[2] else y
}

table1(~ TRTGRP + AGE + SEX.FEMALE + BMI + CDAI_BASELINE + CRP..mgL + HxOfTNFi + 
         STEROID + IMMUNOMOD + LOC.ILEAL | DRUG*TRIAL, 
       data=tbl1,
       overall=FALSE,
       # render = my.render,
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)

#################################################################################################
