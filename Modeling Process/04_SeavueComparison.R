# Overview:
# We trained ADA-attributable model using all ADA participants
# We calculated the remission rate of UST active participants who are TNFi naive
# Using UST placebo participants who are TNFi naive as simulated ADA cohort, we generated predictions of ADA-attributable effect of those participants using ADA-attributable model


# Outputs:
# * remission rates table
# * Table 1-ADA Trial Baseline Characteristics
# * Table 1- UST active vs UST placebo (ADA simulated) Baseline Characteristics 


################################################################################

library(tidyverse)  # For general data processing
library(lme4)       # Mixed effect modeling
library(lmerTest)   # Mixed effect model p-values
library(ggplot2)    # Plotting
library(gridExtra)  # Subplots
library(sjPlot)     # Publication ready tables
library(ggpubr)     # Publication ready figures
library(table1)     # table 1

################################################################################

crohnsData <- read_csv(file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_placebo_imp.csv')

# ADA participants: N = 239 (all active)
ADA <- crohnsData %>%  filter(DRUG == 'ADA')
dim(ADA)

# UST participants that's TNFi naive: N = 284 
# (Placebo: 135; Active: 149)
UST <- crohnsData %>% filter(DRUG == 'UST') %>% 
  # TNF naive
  filter(HxOfTNFi == 0)

dim(UST %>% filter(TRTGRP=='Placebo'))
dim(UST %>% filter(TRTGRP=='Active'))

################################################################################

## ADA Attributable Model

# build formula objects
dep_var <- 'DRUG_REDUCTION'    #DRUG_REDUCTION = Week_8_CDAI_Reductoin - Predicted_Placebo_Effect #see '01_PlaceboAttributableModel.R'
covariate_list <- c('CDAI_BASELINE_CENT','AGE_CENT','BMI_CENT','CRP..mgL_CENT',
                    'SEX.MALE','HxOfTNFi','STEROID','IMMUNOMOD','LOC.ILEAL')    

# lmer DRUG_REDUCTION ~ covariates + (1|TRIAL) 
f1 <- paste(dep_var, paste(covariate_list, collapse = '+'), sep='~')
f2 <- paste0(f1, '+(1|TRIAL)')

# ada model
m2_ada <- lmerTest::lmer(f2, data=ADA)

summary(m2_ada)
ranef(m2_ada)

# save model
PATH = 'G:/Shan/Week 8 Identification/Export Jun23/Paper 1 Week 8/Objects/'
sink(paste(PATH,'m2_ada_readable.txt',sep='')); print(summary(m2_ada)); sink()

# ICC, R2 (sjPlot)
tab_model(m2_ada,
          show.ci=F,
          show.se=T,
          pred.labels = c('Intercept','Baseline CDAI (Centered)','Age (Centered)','BMI (Centered)',
                          'CRP (mg/L) (Centered)','Sex: Male','HxOfTNFi','Steroid Use','Immunomodulator Use','Ileal Disease'),
          dv.labels = c('ADA Attributable (Mixed)'),
          file='G:/Shan/Week 8 Identification/Tables/AdaAttributable.html')

################################################################################

## UST Active TNF-Naive: calculate  Remission Rate

UST_active <- UST %>% 
  filter(TRTGRP=='Active') %>%
  mutate(REMISSION = if_else(CDAIL_WEEK8 < 150, 1, 0))

UST_summary <- UST_active %>% 
  dplyr::summarize(
    DRUG = 'UST Active TNF-Naive',
    N = n(),
    RemCount = sum(REMISSION),
    RemRate  = mean(REMISSION),
    CI_low   = RemRate - qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N),
    CI_high  = RemRate + qt(0.975, df=N)*sqrt((RemRate)*(1-RemRate)/N)
  )


## UST Placebo TNF-Naive (==simulated ADA)

UST_placebo <- UST %>% 
  filter(TRTGRP=='Placebo')

# calculate week 8 cdai if UST placebo tnf-naive patients took ada instead: 
# * ADA-drug-effect-caused-CDAI-reduction = predict (m2_ada, newdata=UST_placebo)
# * simulated ADA cohort week 8 CDAI= CDAI_week8 - ADA-drug-effect-caused-CDAI-reduction

ADA_active <- UST_placebo %>%
  mutate(ADA_WEEK8 = CDAIL_WEEK8 - predict(m2_ada, newdata=UST_placebo, re.form=~0)) %>%
  mutate(REMISSION = if_else(ADA_WEEK8 < 150, 1, 0))

N_ada = 135 # N of simulated ada data 
ADA_summary <- ADA_active %>% 
  dplyr::summarize(
    DRUG = 'ADA Sim. Active TNF-Naive',
    N = n(),
    RemCount = sum(REMISSION),
    RemRate  = mean(REMISSION),
    CI_low   = RemRate - qt(0.975, df=N_ada)*sqrt((RemRate)*(1-RemRate)/N_ada),
    CI_high  = RemRate + qt(0.975, df=N_ada)*sqrt((RemRate)*(1-RemRate)/N_ada)
  )

################################################################################

## Remission Rates

rem_summary <- rbind(UST_summary, ADA_summary)
rem_summary

# Fisher Exact Test
rem_rates <- data.frame(
  #                        ADA              UST
  'remission'     = c(    62,     67),
  'not-remission' = c(135-62, 149-67)
)
rem_rates

fisher.test(rem_rates)

################################################################################

## Characteristic Tables

## ADA Trials
ADA_tbl <- ADA %>% 
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
                             labels=c('Yes','No')))

label(ADA_tbl$CRP..mgL)       <- 'CRP'
units(ADA_tbl$CRP..mgL)       <- 'mg/L'
label(ADA_tbl$AGE)            <- 'Age'
label(ADA_tbl$BMI)            <- 'BMI'

label(ADA_tbl$CDAI_BASELINE)  <- 'Baseline CDAI'
label(ADA_tbl$CDAI_REDUCTION) <- 'CDAI Reduction'

label(ADA_tbl$SEX.FEMALE) <- 'Sex: Female'
label(ADA_tbl$HxOfTNFi)   <- 'HxOfTNFi'
label(ADA_tbl$STEROID)    <- 'Steroid Use'
label(ADA_tbl$IMMUNOMOD)  <- 'Immunomodulator Use'
label(ADA_tbl$LOC.ILEAL)  <- 'Location: Ileal'

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

table1(~ AGE + SEX.FEMALE + BMI + CDAI_BASELINE + CRP..mgL + HxOfTNFi + 
         STEROID + IMMUNOMOD + LOC.ILEAL + CDAI_REDUCTION | TRIAL, 
       data=ADA_tbl,
       overall=FALSE,
       render = my.render,
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)

################################################################################

## UST Active vs UST Placebo (ADA Simulated) -- TNF-naive

# UST Active = UST Active 
# ADA Simulated = UST Placebo
UST_tbl <- UST %>% 
  mutate(POPULATION = if_else(TRTGRP=='Active', 'Observed UST', 'Simulated ADA')) %>% 
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
  # add ada cdai reduction
  mutate(ADA_ATTRIBUTABLE = predict(m2_ada, newdata=UST, re.form=~0)) %>% 
  mutate(CDAI_REDUCTION = if_else(TRTGRP=='Placebo',
                                  CDAI_REDUCTION + ADA_ATTRIBUTABLE,
                                  CDAI_REDUCTION))

label(UST_tbl$CRP..mgL)       <- 'CRP'
units(UST_tbl$CRP..mgL)       <- 'mg/L'
label(UST_tbl$AGE)            <- 'Age'
label(UST_tbl$BMI)            <- 'BMI'

label(UST_tbl$CDAI_BASELINE)  <- 'Baseline CDAI'
label(UST_tbl$CDAI_REDUCTION) <- 'CDAI Reduction'

label(UST_tbl$SEX.FEMALE) <- 'Sex:Female'
label(UST_tbl$HxOfTNFi)   <- 'HxOfTNFi'
label(UST_tbl$STEROID)    <- 'Steroid Use'
label(UST_tbl$IMMUNOMOD)  <- 'Immunomodulator Use'
label(UST_tbl$LOC.ILEAL)  <- 'Location:Ileal'

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

table1(~ AGE + SEX.FEMALE + BMI + CDAI_BASELINE + CRP..mgL + HxOfTNFi + 
         STEROID + IMMUNOMOD + LOC.ILEAL + CDAI_REDUCTION | POPULATION, 
       data=UST_tbl,
       overall=FALSE,
       render = my.render,
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)

################################################################################
