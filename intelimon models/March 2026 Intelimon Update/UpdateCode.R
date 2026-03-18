
# 3/18/2026


# Look at pre/post/ yr-01 burn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pre-Burn
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(leaps)
library(performance)
library(ggplot2)

pred <- read.csv("RMNP Intelimon Update/1. ROMO_East_Fuels/a. ROMO_East_Fuels_PRE_TotalAll/ROMO_East_Fuels_PRE_TotalAll_pred.csv")
vars <- read.csv("RMNP Intelimon Update/1. ROMO_East_Fuels/a. ROMO_East_Fuels_PRE_TotalAll/ROMO_East_Fuels_PRE_TotalAll_vars.csv")
vars <- vars %>% mutate(across(where(is.character), as.factor)) %>% select(where(~!is.factor(.) || nlevels(.) > 1)) %>% select(where(~sum(is.na(.)) == 0))



# Create training/test data
set.seed(123)
sample_size <- floor(0.2* nrow(vars))
test_indices <- sample(seq_len(nrow(vars)), size = sample_size)
train_vars <- vars[-test_indices, ]
test_vars <- vars[test_indices, ]
train_pred <- pred[-test_indices, ]
test_pred <- pred[test_indices, ]

regfit = regsubsets(train_pred$TotalAll~., train_vars,  nvmax = 5 , method = "seqrep", really.big = T) 

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
# 2: MaxSD + hr100_1000_l1_median
# 3: h_l2_median + MaxSD + hr100_1000_l1_median 
# 4: h_l2_median + h_US_cnt+MaxSD + hr100_1000_l1_median 
# 5: h_l2_median + h_US_cnt + MaxSD + hr100_1000_l1_median + hr100_1000_l1_tgi

# 0.05
# 2: hr100_1000_l1_median + hr100_1000_l1_tgi 
# 3: hr0_10_l1_vari + hr100_1000_l1_median + hr100_1000_l1_tgi 
# 4: s_lng_prop_sk + hr0_10_l1_vari + hr100_1000_l1_median + hr100_1000_l1_tgi 
# 5: h_l2_median + s_lng_prop_sk + hr0_10_l1_vari + hr100_1000_l1_median + hr100_1000_l1_tgi


attach(train_vars)

lm1 <- lm(formula = train_pred$TotalAll ~hr100_1000_l1_median + hr100_1000_l1_tgi )



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

saveRDS(lm1, "RMNP Intelimon Update/1. ROMO_East_Fuels/a. ROMO_East_Fuels_PRE_TotalAll/ROMO_East_Fuels_PRE_TotalAll.rda")

