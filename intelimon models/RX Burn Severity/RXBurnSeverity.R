


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

pred <- read.csv("intelimon models/pred_vars_pairs/pred_burnseverity.csv") %>% select(-1)
vars <- read.csv("intelimon models/pred_vars_pairs/vars_burnseverity.csv") %>% select(-1)
temp_data <- cbind(pred, vars)

temp_data <- temp_data %>% select(-c(4:12)) %>% filter(MonStatus != "YR01") 
temp_data <- temp_data %>% slice(-c(19:24)) %>% filter(MonStatus != "PR01") 
temp_data <- temp_data %>% slice(-c(1, 5, 7, 11, 17, 19:21, 25, 34:36))

# Now I have my data filtered the way I want it, separate pred and vars back out
pred <- temp_data %>% select(1:15)
vars <- temp_data %>% select(16:268)

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
regfit = regsubsets(train_pred$Sub5Pct~., train_vars,  nvmax = 3 , method = "seqrep", really.big = T) 

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
lm1 <- lm(formula = train_pred$Sub5Pct ~ h_US_vari + s_l3_prop_mn      +     CBH)

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
data<-data.frame(x=lm_out, y=pred$Sub5Pct)

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
errors<- test_out - test_pred$Sub5Pct
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
rmse 

# format test predicted and observed
data2<- data.frame(x=test_out, y=test_pred$Sub5Pct)

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
BS <- BS %>% select(c(x, y, MacroPlot.Name, MonStatus)) %>% 
  mutate(burn = case_when(MacroPlot.Name %in% c("BME_RAP001", "BME_RAP003", "BME_RAP004", "BME_RAP006", "BME_RAP009", "BME_RAP010", "BME_RAP002") ~ "Fall2024",
                          MacroPlot.Name %in% c("BME_RAP013", "BME_RAP015", "BME_RAP016", "BME_RAP014", "FPIPS1T02:08", "FPIPS1T02:16") ~ "Spring2025",
                          MacroPlot.Name %in% c("FPIPS1T02:14", "FPIPS1T02:01", "FPIPS1T02:05", "FPIPS1T02:02") ~ "Fall2025"))

BS <- BS %>% pivot_longer(cols = c(x,y), names_to = "type", values_to = "value")
BS <- BS %>% mutate(type = case_when(type == "x" ~ "Predicted", type == "y" ~ "Observed"))

sum <- BS %>% group_by( burn, type ) %>% summarise(max = max(value), min = min(value), median = median(value), mean = mean(value), sd = sd(value)) %>% 
  ungroup()


ggplot(BS) + geom_boxplot(aes(x = burn, y = value, fill = type)) 




ggplot(sum) + geom_col(aes(x = burn, y = mean, fill = type), position = position_dodge())+ 
  geom_errorbar(aes(x = burn, ymin = mean - sd, ymax = mean+ sd, fill = type), position = position_dodge())




