# Overview:
# * We checked the model assumptions for the Placebo Attributable Model: 
# ** Independence: Residual vs Index
# ** Equal Variance: Residuals vs Fitted Values
# ** Normality: qq plot 
# ** Linearity: Residuals vs Covariates

# Outputs:
# *Extended Data Figure 5: Model Assumptions

################################################################################

library(tidyverse)
library(lme4)
library(ggplot2)
library(ggpubr)
library(patchwork)

################################################################################

# N = 1310
placebo_df <- read.csv(file = "G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_imputed.csv") %>% 
  filter(TRTGRP == 'Placebo')
 
rm <- readRDS(file='G:/All Analysis/Vignesh/Export/Paper 1 Week 8/Objects/placebo_mixed.rds')

# linearity assumption check
placebo_df$resid <- resid(rm)
placebo_df$fitted <- fitted(rm)

################################################################################

norm.check <- data.frame(index = 1:1310,
                         resid = resid(rm),
                         fitted = fitted(rm))

resid.sd <- sd(resid(rm))

p1 <- ggplot(norm.check, aes(x=index, y=resid)) +
  geom_point() + 
  geom_hline(yintercept=0, linetype='dashed') + 
  geom_hline(yintercept=2*resid.sd, linetype='dashed') + 
  geom_hline(yintercept=-2*resid.sd, linetype='dashed') + 
  ggtitle('Check for Independence') + 
  xlab('Index') + 
  ylab('Residual')

p2 <- ggplot(norm.check, aes(x=fitted, y=resid)) +
  geom_point() + 
  geom_hline(yintercept=0, linetype='dashed') + 
  geom_hline(yintercept=2*resid.sd, linetype='dashed') + 
  geom_hline(yintercept=-2*resid.sd, linetype='dashed') +
  ggtitle('Check for Equal Variance') + 
  xlab('Fitted') + 
  ylab('Residual')

p3 <- ggqqplot(norm.check$resid) + 
  ggtitle('Check for Normality') + 
  xlab('Theoretical') + 
  ylab('Sample')

p4 <- placebo_df %>%
  mutate(CDAI_baseline = CDAI_BASELINE_CENT + 300,
         Age = AGE_CENT + 35,
         BMI = BMI_CENT + 20,
         CRP = CRP..mgL_CENT + 10) %>%
  select(c(CDAI_baseline, Age, BMI, CRP, resid)) %>%
  gather(-resid, key = 'some_var_name', value='some_value_name') %>%
  ggplot(aes(x=some_value_name, y=resid)) +
    geom_point() + 
    facet_wrap(~ some_var_name, scales='free') + 
    ggtitle('Check for Linearity') + 
    xlab('Covariates') + 
    ylab('Residual')

pFinal <- (p1 + p2) / (p3 + p4)
pFinal

# Save the joint plot to file
ggsave(filename = "G:/Shan/Week 8 Identification/Figures/NormalityAssumptions.png", 
       plot = pFinal, height = 6, width = 9)

ggsave(filename = "G:/Shan/Week 8 Identification/Figures/NormalityAssumptions_Independence.png", 
       plot = p1, height = 6, width = 9)

ggsave(filename = "G:/Shan/Week 8 Identification/Figures/NormalityAssumptions_Homoscedasticity.png", 
       plot = p2, height = 6, width = 9)

ggsave(filename = "G:/Shan/Week 8 Identification/Figures/NormalityAssumptions_Normality.png", 
       plot = p3, height = 6, width = 9)

ggsave(filename = "G:/Shan/Week 8 Identification/Figures/NormalityAssumptions_Linearity.png", 
       plot = p4, height = 6, width = 9)


#####################################################################################################
