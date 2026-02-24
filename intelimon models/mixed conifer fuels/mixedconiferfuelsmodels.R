# Edited by: Christina Fossum 2/09/2026
# Title: Intelimon Modeling - DEMO

# The purpose of this script is to walk you through the steps of modeling 'pred' metrics (traditional field-collected data) using 'vars' intelimon data metrics
# Ultimately, the goal is to produce good predictive models that we can send into intelimon so that they can provide us with more useful outputs from our TLS scans
# Basic workflow: 
# (1) we create reasonable models (i.e. minimal variables (2-3), good enough r^2, good enough rmse, and that these both translate from training data to test data)
# (2) Save models, send to intelimon staff
# (3) Whenever they do their updates, they will review these models and if they pass review process, they will integrate into our available intelimon TLS metrics


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1
# Load required r packages

library(tidyverse)
library(leaps)
library(performance)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2
# Load data 

pred <- read.csv("intelimon models/mixed conifer fuels/mixedconifer_fuels_pred.csv") %>% select(-1) %>% slice(-15)
vars <- read.csv("intelimon models/mixed conifer fuels/mixedconifer_fuels_vars.csv") %>% select(-1) %>% slice(-15)


# Run the following chunk so that the vars data is formatted correctly for modeling
vars <- vars %>% mutate(across(where(is.character), as.factor)) %>% select(where(~!is.factor(.) || nlevels(.) > 1)) %>% select(where(~sum(is.na(.)) == 0))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3
## Split data into training and test datasets (change sample_size if needed so that 'test_vars' and 'test_pred' contain atleast 3-4 datapoints)

set.seed(123)
sample_size <- floor(0.2 * nrow(vars))
test_indices <- sample(seq_len(nrow(vars)), size = sample_size)
train_vars <- vars[-test_indices, ]
test_vars <- vars[test_indices, ]
train_pred <- pred[-test_indices, ]
test_pred <- pred[test_indices, ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4
## Model Selection

# I want to try to model substrate burn severity, so I will set variable to 'SubSev'. nvmax is max number of variables (set to 2 or 3)
# method should be kept as 'seqrep' and really.big as 'T' although these can both be changed 
regfit = regsubsets(train_pred$TotalAll~., train_vars,  nvmax = 3 , method = "seqrep", really.big = T) 

# Visual checks to see how many variables you should use
regsum<- summary(regfit)
plot(regsum$rsq, type = "l") 
plot(regsum$rss, type = "l") 
plot(regsum$bic, type = "l") 
which.min(regsum$bic) 
plot(regfit, scale = "r2") 

# Change # to best model
vcov(regfit,2) 
coef(regfit,3)

# Whichever variables are listed for vcov and coef, these are your variables selected for models 
# Here I got: h_l2_std, s_l4_na_per, and s_l4_prop_ku

attach(train_vars)
# put the selected predictor variables in the model after the "~" seperated by "+" and run the model
lm1 <- lm(formula = train_pred$TotalAll ~ h_l5_median+ s_l2_na_per   +   LF_CBD)
lm1 <- lm(formula = train_pred$TotalAll ~ h_OS_per   +    LF_CBD)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 5
# Check model 


#look at the summary r2 values and RMSE
summary(lm1)
r2(lm1) 
rmse(lm1) 

# Visual check of model assumptions. Make sure Linearity, homogeneity of variance, influential observations, colinearity, normality of risiduals all check out
check_model(lm1) #if you have colinearity, consider using model with fewer variables

#use the model to predict out new values from variables alone
lm_out <- predict.lm(lm1, vars)

#format the predictions and observations for plotting
data<-data.frame(x=lm_out, y=pred$TotalAll)

# Plot observed values vs. predicted
ggplot(data,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='TotalFuel Load') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 
# data points should be distributed fairly evenly along lm line

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 6
# Check model's ability to predict test data

# use the model to predict test data
test_out<- predict.lm(lm1,test_vars)

# Calculate RMSE of test data
errors<- test_out - test_pred$TotalAll
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
rmse 

# format test predicted and observed
data2<- data.frame(x=test_out, y=test_pred$TotalAll)

# Plot test observed values vs. test predicted
ggplot(data2,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='test') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 

# this plot really is not great: points are not distributed evenly along lm

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If the model looks good, save to disk using naming convention you will remember meaning of
saveRDS(lm1,  "intelimon models/DEMO_ROMO_SubSev.rda")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot your model to see where you went wrong/show efficacy. y = actual data, x = predicted data 


BS <- cbind(pred, data)



BS <- BS %>% pivot_longer(cols = c(x,y), names_to = "type", values_to = "value")


ggplot(BS) + geom_col(aes(x = Macroplot, y = value, fill = type), position = position_dodge()) 


sum <- BS %>% group_by( burn, type ) %>% summarise(max = max(value), min = min(value), median = median(value), mean = mean(value), sd = sd(value)) %>% ungroup()

ggplot(sum) + geom_col(aes(x = burn, y = mean, fill = type), position = position_dodge())+ geom_errorbar(aes(x = burn, ymin = mean - sd, ymax = mean+ sd, fill = type), 
                                                                                                         position = position_dodge())













