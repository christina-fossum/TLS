

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1
# Load required r packages

library(tidyverse)
library(leaps)
library(performance)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2
# Load data (using burn severity data for this demo)

pred <- read.csv("intelimon models/pred_vars_pairs/pred_dino_bighornbasin.csv") %>% select(-1)
vars <- read.csv("intelimon models/pred_vars_pairs/vars_dino_bighornbasin.csv") %>% select(-1)

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
regfit = regsubsets(train_pred$Cover_Shrub~., train_vars,  nvmax = 3 , method = "seqrep", really.big = T) 

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


attach(train_vars)
# put the selected predictor variables in the model after the "~" seperated by "+" and run the model
lm1 <- lm(formula = train_pred$Cover_Shrub ~ vox_l1_tgi + fine_l1_vari)

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
data<-data.frame(x=lm_out, y=pred$Cover_Shrub)

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
errors<- test_out - test_pred$Cover_Shrub
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
rmse 

# format test predicted and observed
data2<- data.frame(x=test_out, y=test_pred$Cover_Shrub)

# Plot test observed values vs. test predicted
ggplot(data2,aes(x,y)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE,) +
  theme_light() +
  labs(x='Predicted Values', y='Observed Values', title='test') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold') 
  ) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If the model looks good, save to disk using naming convention you will remember meaning of
saveRDS(lm1,  "intelimon models/DEMO_ROMO_SubSev.rda")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot your model to see where you went wrong/show efficacy. y = actual data, x = predicted data 
library(stringr)

BS <- cbind(pred, data)
BS <- BS %>% select(c(x, y, MacroPlot.Name, scan_position))

BS <- BS %>% pivot_longer(cols = c(x,y), names_to = "type", values_to = "value")

BS <- BS %>% mutate(Treatment = case_when(str_detect(MacroPlot.Name, "Control") ~ "Control",
                                          str_detect(MacroPlot.Name, "Experimental") ~ "Experimental"))

ggplot(BS) + geom_boxplot(aes(x = Treatment, y = value, fill = type)) 


sum <- BS %>% group_by( Treatment, type ) %>% summarise(max = max(value), min = min(value), median = median(value), mean = mean(value), sd = sd(value)) %>% ungroup()

ggplot(sum) + geom_col(aes(x = Treatment, y = mean, fill = type), position = position_dodge())+ geom_errorbar(aes(x = Treatment, ymin = mean - sd, ymax = mean+ sd, fill = type), 
                                                                                                         position = position_dodge())




