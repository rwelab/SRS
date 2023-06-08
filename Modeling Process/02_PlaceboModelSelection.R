# Overview:
# * We compared the predictive performance of the mixed effect model for placebo effect against other 
#   statistical and machine learning models
# * We compared fixed effect linear model, mixed effect linear model. Lasso, Random Forest and XGBoost
# * We used 5-fold cross-validation to calculate average RMSE as the performance metric 

# Outputs:
# * Extended Table 2:	Model Selection Results for Placebo Effects

################################################################################

library(tidyverse)  # For general data processing
library(lme4)       # Mixed effect modeling
library(lmerTest)   # Mixed effect model p-values
library(caret)      # RMSE calculation
library(glmnet)     # LASSO regression
library(randomForest)
library(xgboost)

################################################################################

## Placebo Data


placebo_df <- read_csv(file='G:/Shan/Week 8 Identification/CombinedTrials/crohnsData_wk8_imputed.csv') %>% 
  filter(TRTGRP=='Placebo') %>% 
  select(TRIAL, YEAR_CENT, CDAI_BASELINE_CENT, AGE_CENT, BMI_CENT, CRP..mgL_CENT, 
         SEX.MALE, HxOfTNFi, STEROID, IMMUNOMOD, LOC.ILEAL, CDAI_REDUCTION)

################################################################################

## Cross Validation Function

# stratify dataset
set.seed(123)
folds <- createFolds(as.factor(placebo_df$CDAI_REDUCTION), k=5)

# Folds are stratified
# placebo_df[folds$Fold1,] %>% summarise_all(list(mean))
# placebo_df[folds$Fold2,] %>% summarise_all(list(mean))
# placebo_df[folds$Fold3,] %>% summarise_all(list(mean))
# placebo_df[folds$Fold4,] %>% summarise_all(list(mean))
# placebo_df[folds$Fold5,] %>% summarise_all(list(mean))


covariate_list <- placebo_df %>% select(-c('TRIAL','CDAI_REDUCTION')) %>% colnames()
f1 <- paste('CDAI_REDUCTION', paste(covariate_list, collapse = '+'), sep='~') # fixed effect formula
f2 <- paste0(f1, '+(1|TRIAL)') # random effect formula

myCV <- function(data, folds, k, method, param){
  rmse_out <- c()
  for (i in 1:k){
    train <- data[-folds[[i]], ]
    test  <- data[ folds[[i]], ]
    
    if (method=='lm'){
      cv.fit  <- lm(f1, data=train)
      cv.pred <- predict(cv.fit, newdata=test)  
    }
    else if (method=='lmer'){
      cv.fit  <- lmer(f2, data=train)
      cv.pred <- predict(cv.fit, newdata=test, re.form=~0) 
    }
    else if (method=='lasso'){
      cv.fit  <- glmnet(x=data.matrix(train[,2:11]),
                        y=data.matrix(train[,  12]),
                        alpha=1,
                        lambda=param$lambda)
      cv.pred <- predict(cv.fit, s=param$lambda, newx=data.matrix(test[,2:11])) 
    }
    else if (method=='rf'){
      # ntree = 500
      cv.fit  <- randomForest(CDAI_REDUCTION~., data=train[,2:12], mtry=2)
      cv.pred <- predict(cv.fit, newdata=test[,2:11])  
    }
    else if (method=='xgboost'){
      # convert train to DMatrix for xgboost
      dtrain <- xgb.DMatrix(label = train$CDAI_REDUCTION, data = as.matrix(train[,2:11]))
      dtest  <- xgb.DMatrix(label = test$CDAI_REDUCTION,  data = as.matrix(test[,2:11]))

      # using known best parameters
      cv.fit <- xgboost::xgb.train(data = dtrain,
                                   objective = 'reg:squarederror',
                                   eval_metric = 'rmse',
                                   nrounds   = param$nrounds,
                                   max_depth = param$max_depth,
                                   eta       = param$eta,
                                   # rate_drop = param$rate_drop,
                                   # skip_drop = param$skip_drop,
                                   min_child_weight = param$min_child_weight,
                                   gamma     = param$gamma,
                                   colsample_bytree = param$colsample_bytree)
      cv.pred <- predict(cv.fit, dtest) 
    }
    else{
      print('Wrong method input')
      return(-1)
    }

    rmse_out <- c(rmse_out, caret::RMSE(test$CDAI_REDUCTION, cv.pred))
  }
  return(rmse_out)
}

################################################################################

## Linear Regression

set.seed(8)

lm.results = myCV(data = placebo_df, folds=folds, k=5, method='lm')
lm.results
mean(lm.results)

################################################################################

## Mixed Effect Linear Regression

set.seed(8)

lmer.results = myCV(data = placebo_df, folds=folds, k=5, method='lmer')
lmer.results
mean(lmer.results)

################################################################################

## Linear LASSO Regression

set.seed(8)

bstGlm <- cv.glmnet(x=data.matrix(placebo_df[,2:11]), y=data.matrix(placebo_df[,12]), alpha=1)
param = expand.grid(lambda=bstGlm$lambda.min)

lasso.results = myCV(data = placebo_df, folds=folds, k=5, method='lasso', param=param)
lasso.results
mean(lasso.results)

################################################################################

## Random Forest

set.seed(8)

rf.results = myCV(data = placebo_df, folds=folds, k=5, method='rf')
rf.results
mean(rf.results)

################################################################################

## XGBoost (DART)

set.seed(8)

# Best hyperparameters from caret package (see below)
# nrounds = 150, max_depth = 2, eta = 0.3, gamma = 0, subsample = 0.5, colsample_bytree = 0.8, 
# rate_drop = 0.5, skip_drop = 0.05 and min_child_weight = 1

param = expand.grid(nrounds = 150, max_depth = 2, eta = 0.3, rate_drop = 0.5, skip_drop = 0.05,
                    min_child_weight = 1, gamma = 0, subsample = 0.5, colsample_bytree = 0.8)

xgb.results = myCV(data = placebo_df, folds=folds, k=5, method='xgboost', param=param)
xgb.results
mean(xgb.results)

################################################################################

