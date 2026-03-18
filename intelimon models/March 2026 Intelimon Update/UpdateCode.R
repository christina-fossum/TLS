
# 3/18/2026


# Look at pre/post/ yr-01 burn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# tree count
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(leaps)
library(performance)
library(ggplot2)

pred <- read.csv("RMNP Intelimon Update/2. ROMO_East_Trees/a. ROMO_East_TotalTreeCount/ROMO_EAST_TotalTreeCount_pred.csv")
vars <- read.csv("RMNP Intelimon Update/2. ROMO_East_Trees/a. ROMO_East_TotalTreeCount/ROMO_EAST_TotalTreeCount_vars.csv")
vars <- vars %>% mutate(across(where(is.character), as.factor)) %>% select(where(~!is.factor(.) || nlevels(.) > 1)) %>% select(where(~sum(is.na(.)) == 0))



# Create training/test data
set.seed(123)
sample_size <- floor(0.2* nrow(vars))
test_indices <- sample(seq_len(nrow(vars)), size = sample_size)
train_vars <- vars[-test_indices, ]
test_vars <- vars[test_indices, ]
train_pred <- pred[-test_indices, ]
test_pred <- pred[test_indices, ]

regfit = regsubsets(train_pred$TreeCount~., train_vars,  nvmax = 5 , method = "seqrep", really.big = T) 

# Visual checks to see how many variables you should use
regsum<- summary(regfit)
plot(regsum$rsq, type = "l") 
plot(regsum$rss, type = "l") 
plot(regsum$bic, type = "l") 
which.min(regsum$bic) 
plot(regfit, scale = "r2") 

# Change # to best model
vcov(regfit,2) 
coef(regfit,2)

# 0.1 for test data
# 2: h_MS_cnt + fuel0_3l1_kurt  
# 3:  h_MS_cnt + h_MS_per + fuel0_3l1_kurt 
# 4: h_MS_cnt + h_zq20 + vox_l1_cnt + fuel0_3l1_kurt 
# 5: h_l5_std + h_MS_cnt + h_zq20 + fDim + fuel0_3l1_kurt 



attach(train_vars)

lm1 <- lm(formula = train_pred$TotalAll ~   h_ground_cnt + s_lng_zero_per + shrubs_l1_skew)



summary(lm1)
r2(lm1) 
rmse(lm1) 


check_model(lm1) 

lm_out <- predict.lm(lm1, vars)

data<-data.frame(x=lm_out, y=pred$TotalAll)

ggplot(data,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='TotalFuel Load') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 

test_out<- predict.lm(lm1,test_vars)

errors<- test_out - test_pred$TotalAll
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
rmse 

data2<- data.frame(x=test_out, y=test_pred$TotalAll)

ggplot(data2,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='test') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 

saveRDS(lm1, "RMNP Intelimon Update/1. ROMO_East_Fuels/b. ROMO_East_Fuels_POST_TotalAll/ROMO_East_Fuels_POST_TotalAll.rda")

