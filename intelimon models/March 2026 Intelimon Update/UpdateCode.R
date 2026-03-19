
# 3/18/2026


# Look at pre/post/ yr-01 burn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Substrate Burn Severity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(leaps)
library(performance)
library(ggplot2)

pred <- read.csv("RMNP Intelimon Update/3. ROMO_East_RX/a. ROMO_East_RX_SubSev/ROMO_East_RX_SubSev_pred.csv")
vars <- read.csv("RMNP Intelimon Update/3. ROMO_East_RX/a. ROMO_East_RX_SubSev/ROMO_East_RX_SubSev_vars.csv")
vars <- vars %>% mutate(across(where(is.character), as.factor)) %>% select(where(~!is.factor(.) || nlevels(.) > 1)) %>% select(where(~sum(is.na(.)) == 0))



# Create training/test data
set.seed(123)
sample_size <- floor(0.3* nrow(vars))
test_indices <- sample(seq_len(nrow(vars)), size = sample_size)
train_vars <- vars[-test_indices, ]
test_vars <- vars[test_indices, ]
train_pred <- pred[-test_indices, ]
test_pred <- pred[test_indices, ]

regfit = regsubsets(train_pred$SubSev~., train_vars,  nvmax = 5 , method = "seqrep", really.big = T) 

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
# 2: s_l2_zero_per + hr0_10_l1_median
# 3:   h_l4_per     +      MeanSA + hr0_10_l1_median
# 4: h_l4_per   +     h_MS_kurt  +  s_l2_zero_per+ hr0_10_l1_median 
# 5: h_l4_per + h_MS_kurt + s_l2_zero_per + LAI + hr0_10_l1_median

# 0.2 for test data
# 2: s_l2_zero_per + hr0_10_l1_median
# 3:  s_l2_zero_per +hr0_10_l1_median     +      LF_CBD
# 4:  h_US_tgi    +    h_MS_kurt   +  s_l1_prop_ku +hr0_10_l1_median 
# 5: h_US_tgi    +    h_MS_mean    +    h_MS_kurt  +   s_l1_prop_ku +hr0_10_l1_median 


# 0.3 for test data
# 2: h_l4_skew + s_l4_prop_ku
# 3:   h_l4_skew + s_l4_prop_ku + fine_l1_tgi
# 4: h_l1_kurt  +  h_l2_std + h_US_median + s_l4_na_per 
# 5:   h_l1_kurt   +   h_l2_std +  h_US_median    +  h_US_tgi +  s_l4_na_per


attach(train_vars)

lm1 <- lm(formula = train_pred$SubSev ~    hr0_10_l1_median + s_l2_zero_per + h_MS_kurt )



summary(lm1)
r2(lm1) 
rmse(lm1) 


check_model(lm1) 

lm_out <- predict.lm(lm1, vars)

data<-data.frame(x=lm_out, y=pred$SubSev)

ggplot(data,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='?') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 

test_out<- predict.lm(lm1,test_vars)

errors<- test_out - test_pred$SubSev
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
rmse 

data2<- data.frame(x=test_out, y=test_pred$SubSev)

ggplot(data2,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='test') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 

saveRDS(lm1, "RMNP Intelimon Update/3. ROMO_East_RX/a. ROMO_East_RX_SubSev/ROMO_East_RX_SubSev.rda")

