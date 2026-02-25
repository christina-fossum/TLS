
# Edited by: Christina Fossum 2/24/2026


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1
# Load required r packages

library(tidyverse)
library(leaps)
library(performance)
library(see)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: fuels

pred <- read.csv("DATA/pred_clean/fuels.csv")%>% select(-1)
vars <- read.csv("DATA/vars_raw/romo_intelimon2026/merged_metrics.csv") %>% select(-1) %>% rename(scan_name = h_filename)
temp <- inner_join(pred, vars)

#East side pre-burn fuels
pre1 <- temp %>% filter(MonStatus %in% c("00PRE", "01PRE", "02PR01", "02PRE", "01PR01", "00Pre", "01YR15b", "02YR15b", "01YR21b", "01YR16b")) %>% 
  slice(-c(1:24))
pre2 <- temp %>% filter(!MonStatus %in% c("01POST", "02POST", "3POST")) %>% slice(-c(1:24)) #include yr01
pre3 <- temp %>% slice(-c(1:24)) #include postburn
#East Side Post-Burn Fuels
pbf <- temp %>% filter(MonStatus %in% c("01POST", "02POST", "3POST")) %>% filter(!Macroplot %in% c("AP_02", "AP_03", "AP_04", "AP_05", "AP_06"))
#East Side yr01 fuels
yr1 <- temp %>% filter(MonStatus == "01YR01")
# West Side Fuels
wsf <- temp %>% slice(c(1:14, 15,17,19,21,23))


#format pred/vars (change this for each above)
pred <- pre2 %>% select(names(pred))
vars <- pre2 %>% select(names(vars))

vars <- vars %>% select(-c(1:3))
vars <- vars %>% mutate(across(where(is.character), as.factor)) %>% select(where(~!is.factor(.) || nlevels(.) > 1)) %>% select(where(~sum(is.na(.)) == 0))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3
## Split data into training and test datasets (change sample_size if needed so that 'test_vars' and 'test_pred' contain atleast 3-4 datapoints)

set.seed(123)
sample_size <- floor(0.1 * nrow(vars))
test_indices <- sample(seq_len(nrow(vars)), size = sample_size)
train_vars <- vars[-test_indices, ]
test_vars <- vars[test_indices, ]
train_pred <- pred[-test_indices, ]
test_pred <- pred[test_indices, ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4
## Model Selection


regfit = regsubsets(train_pred$TotalAll~., train_vars,  nvmax = 3 , method = "seqrep", really.big = T) 

# Visual checks to see how many variables you should use
regsum<- summary(regfit)
plot(regsum$rsq, type = "l") 
plot(regsum$rss, type = "l") 
plot(regsum$bic, type = "l") 
which.min(regsum$bic) 
plot(regfit, scale = "r2") 

# Change # to best model
vcov(regfit,3) 
coef(regfit,2)

# Whichever variables are listed for vcov and coef, these are your variables selected for models 


attach(train_vars)
# put the selected predictor variables in the model after the "~" seperated by "+" and run the model
lm1 <- lm(formula = train_pred$TotalAll ~  h_l2_cnt   +   h_l3_cnt+ hr100_1000_l1_cnt)

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
  labs(x='Predicted Values', y='Observed Values', title='total fuel load') +
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
sqrt(mse)
rmse #rmse = 0.62. You can expect rmse for test data to not be as good as that for training data, but should still be good enough to predict. 0.62 is not great 

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If the model looks good, save to disk using naming convention you will remember meaning of
saveRDS(lm1, "TLS/intelimon models/intelimon_update/pipo_nonburn_totalfuelload/pipo_totalfuels_not_postburn.rda")
write.csv(pred, "TLS/intelimon models/intelimon_update/pipo_nonburn_totalfuelload/pipo_totalfuels_not_postburn_pred.csv")
write.csv(vars, "TLS/intelimon models/intelimon_update/pipo_nonburn_totalfuelload/pipo_totalfuels_not_postburn_vars.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot your model to see where you went wrong/show efficacy. y = actual data, x = predicted data 


BS <- cbind(pred, data)
BS <- BS %>% select(c(Macroplot, MonStatus, scan_name, x, y))

BS <- BS %>% mutate(group = case_when(Macroplot %in% c("BME_RAP003", "BME_RAP004", "BME_RAP006", "BME_RAP009", "BME_RAP010") ~ "Fall24",
                                      Macroplot %in% c("FPIPO1T09:27", "FPIPS1T02:01", "FPIPS1T02:02", "FPIPS1T02:05", "FPIPS1T02:14", "FPIPO1T09:09") ~ "Fall25",
                                      Macroplot %in% c("BME_RAP013", "BME_RAP014", "BME_RAP015", "BME_RAP016", "FPIPS1T02:08", "FPIPS1T02:16") ~ "Spring25",
                                      Macroplot %in% c("EC_01", "EC_02", "EC_03", "EC_04", "EC_05") ~ "EagleCliff", .default = "other"
                                      ))


BS <- BS %>% pivot_longer(cols = c(x,y), names_to = "data", values_to = "value")

ggplot(subset(BS, group == "other"), aes(x = MonStatus, y = value, fill = data)) + geom_boxplot()




sum <- BS %>% group_by(MonStatus, group, data) %>% summarise(mean = mean(value), sd = sd(value)) %>% ungroup()




ggplot(BS, aes(x = MonStatus, y = value, fill = data)) + geom_boxplot()


ggplot(sum, aes(x = MonStatus, y = mean, fill = data)) + geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean+ sd), position = position_dodge())







