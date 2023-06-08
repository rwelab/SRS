# Overview:
# Training placebo model and imputing placebo attributable effect for all participants. 
# Leave-one-out analysis done to check for residual normality. 

# Outputs:
# * LeaveOneTrialOut plots (QQ, residual)
# * crohnsData_imputed

library(tidyverse)  # For general data processing
library(lme4)       # Mixed effect modeling
# install.packages('lmerTest')
library(lmerTest)   # Mixed effect model p-values
library(caret)
library(ggplot2)    # Plotting
library(gridExtra)  # Subplots
library(sjPlot)     # Publication ready tables
library(ggpubr)     # Publication ready figures
library(sjstats)    # ICC, R2

#####################################################################################################

crohnsData <- read_csv(file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_imputed.csv')

crohnsData %>% colnames()

placebo_df <- crohnsData %>% 
  filter(TRTGRP=='Placebo') %>% 
  # sensitivity analysis
  #filter(DRUG != 'UST') %>% 
  select(TRIAL, YEAR_CENT, CDAI_BASELINE_CENT, AGE_CENT, BMI_CENT, CRP..mgL_CENT, 
         SEX.MALE, HxOfTNFi, STEROID, IMMUNOMOD, LOC.ILEAL, CDAI_REDUCTION) 
  
#####################################################################################################

## Fixed vs Mixed (Baseline/Intercept)


# baseline fixed, mixed models
bfm <- lm(CDAI_REDUCTION ~ 1, data=placebo_df)
bmm <- lmerTest::lmer(CDAI_REDUCTION ~ (1|TRIAL), data=placebo_df)

summary(bfm)
summary(bmm)
ranef(bmm)

#####################################################################################################

## Fixed vs Mixed (With Covariates)


# build formula objects - makes predict(lmer) work
dep_var <- 'CDAI_REDUCTION'
covariate_list <- c('YEAR_CENT','CDAI_BASELINE_CENT','AGE_CENT','BMI_CENT','CRP..mgL_CENT',
                    'SEX.MALE','HxOfTNFi','STEROID','IMMUNOMOD','LOC.ILEAL')

# formula 1: lm CDIA_REDUCTION ~ covariates 
f1 <- paste(dep_var, paste(covariate_list, collapse = '+'), sep='~')
# formula 2: lmer CDAI_REDUCTION ~ covariates + (1|TRIAL)
f2 <- paste0(f1, '+(1|TRIAL)')

# fixed, mixed models
fm <- lm(f1, data=placebo_df)
mm <- lmerTest::lmer(f2, data=placebo_df)

# goodness of fit
anova(mm, fm)

#####################################################################################################

# R-outputs
summary(fm)
summary(mm)
ranef(mm)

# Save models - human readable
PATH = 'G:/Shan/Week 8 Identification/Export Jun23/Paper 1 Week 8/Objects/'
sink(paste(PATH,'m1_plac_fix_readable.txt',sep='')); print(summary(fm)); sink()
sink(paste(PATH,'m1_plac_mix_readable.txt',sep='')); print(summary(mm)); sink()

#####################################################################################################

# ICC, R2 (sjPlot)
tab_model(fm, mm,
          show.ci=F,
          show.se=T,
          pred.labels = c('Intercept','Year (Centered)','Baseline CDAI (Centered)','Age (Centered)','BMI (Centered)',
                          'CRP (mg/L) (Centered)','Sex: Male','HxOfTNFi','Steroid Use','Immunomodulator Use','Ileal Disease'),
          dv.labels = c('Fixed','Mixed'),
          file='G:/Shan/Week 8 Identification/Tables/PlaceboFixedMixed.html')

#####################################################################################################

## Leave One Out (LOO) Analysis

set.seed(8)

# create index for TRIALs -- stratify by TRIAL
placebo_df$ind <- as.numeric(as.factor(placebo_df$TRIAL))

# leave-one-out
trial_pred_info <- NULL
for(i in 1:6) {
  # train, test
  train <- placebo_df[placebo_df$ind != i,]
  test  <- placebo_df[placebo_df$ind == i,]
  
  # fit, predict on trial left out
  fit  <- lmer(f2, data=train)
  pred <- predict(fit, newdata=test, re.form=~0)
  
  # prediction results (y - y_hat)
  TRIAL <- test$TRIAL[1]
  pred.cdai_red.mean <- mean(pred)
  pred.cdai_red.se   <- sd(pred)
  rmse <- caret::RMSE(pred, test$CDAI_REDUCTION)
  
  # concatenate results
  trial_pred_info <- rbind(trial_pred_info, 
                           data.frame(TRIAL, pred.cdai_red.mean, pred.cdai_red.se, rmse))
}

# true vs predicted
loo_results <- placebo_df %>% 
  group_by(TRIAL) %>%
  dplyr::summarize(n = n(), 
            true.cdai_red.mean = mean(CDAI_REDUCTION),
            true.cdai_red.se   = sd(CDAI_REDUCTION)) %>% 
  full_join(trial_pred_info, by='TRIAL')
loo_results

# write.csv(loo_results,
#           file='G:/Shan/Week 8 Identification/Tables/LeaveOneOutTruePredResults.csv',
#           row.names=FALSE)

#####################################################################################################

## Plots
levels.y <- c('CERTIFI','UNITI1','UNITI2','ENACT','ENCORE','PRECISE1')

tidy_results <- loo_results %>% 
  select(-c(n, true.cdai_red.se, pred.cdai_red.se, rmse)) %>% 
  pivot_longer(!TRIAL, names_to='type', values_to='mean') %>% 
  mutate(type = replace(type, type == 'true.cdai_red.mean', 'true')) %>%
  mutate(type = replace(type, type == 'pred.cdai_red.mean', 'pred')) %>% 
  mutate(TRIAL = factor(TRIAL, levels=levels.y))

p1 <- ggplot(tidy_results, aes(x=TRIAL, y=mean)) +
  geom_point(aes(color=type), size=3) + geom_line() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
  labs(title='Actual vs Predicted Mean Placebo CDAI Reduction by Trial') + 
  ylab('CDAI Reduction')

p2 <- ggqqplot(loo_results$pred.cdai_red.mean - loo_results$true.cdai_red.mean,
               title = 'Mean Residual by Trial')

pFinal <- grid.arrange(p1, p2, nrow=2)
ggsave(filename = "G:/Shan/Week 8 Identification/Figures/LeaveOneTrialOutPlots.png", 
       plot = pFinal, height = 9, width = 9)

#####################################################################################################

# Shapiro-Wilks
shapiro.test(loo_results$pred.cdai_red.mean - loo_results$true.cdai_red.mean)

#####################################################################################################

## Placebo Imputation (all trials)

# CDAI_BASELINE - CDAI_WEEK8 = CDAI_REDUCTION
# CDAI_REDUCTION = PLACEBO_REDUCTION + DRUG_REDUCTION

# Reduction = predicted against (lmer(*_REDUCTION ~ ))
# Attributable = predicted results (*_ATTRIBUTABLE = predict())
# Both mean same thing (amount cdai reduced), but using this as naming convention

# CDAI_REDUCTION (TOTAL) = DRUG_REDUCTION + PLACEBO_REDUCTION (predicted; mm)

crohnsData_imputed <- crohnsData %>% 
  mutate(PLACEBO_ATTRIBUTABLE = predict(mm, newdata=crohnsData, re.form=~0)) %>% 
  mutate(DRUG_REDUCTION       = CDAI_REDUCTION - PLACEBO_ATTRIBUTABLE)

write.csv(crohnsData_imputed,
          file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_placebo_imp.csv',
          row.names=FALSE)

#####################################################################################################
