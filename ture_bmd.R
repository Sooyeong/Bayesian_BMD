# Calculate true BMD value & Compare BBMD & BMDS from data


# True BMD is calculated based on a grid search approach

### Generate Data for Simulation - BMDS, BMDL ### 
###
### Written by SL ###

# N_i is fixed as 50 
# Four dose levels - 0, 25, 50, 75


library(ggplot2)
library(dplyr)
#rm(list=ls())
getwd()
setwd('/Users/sooyeonglim/Desktop/r_temp')        

# Set Seed
set.seed(1234)

# d was normalized in shao's fitting example

BMR<-0.1
d_origin<-seq(0,75,1)

# Posterior weight is calculated based on BBMD fitting with an example data from Shao's webpage
post_weight<-c(0.207,0.346,0.115,0.0468,0.00224,0.0953,0.125,0.0625)

# Find a ture BMD from
model_average<-function(d_input, post_weight, d_origin){
  # Since the estimated parameter values were normalized.
  d_norm=(d_input-min(d_origin))/(max(d_origin)-min(d_origin))
  
  d<-d_norm
  # model_1 -Logistic
  a_1<- -4.08
  b_1<- 4.37
  f_1<-1/(1+exp(-a_1-b_1*d))
  
  # model_2 -Probit
  a_2<--2.29
  b_2<-2.44
  
  f_2<-pnorm(a_2+b_2*d, mean=0, sd=1, lower.tail=T, log.p=F)
  

  
  # model_3 -LogLogistic
  a_3=0.03
  b_3=3.49
  c_3=0.2
  
  f_3=a_3+(1-a_3)/(1+exp(-c_3-b_3*log(d)))
  
  # model_4- Dichotomous Hill
  a_4=0.54
  b_4=9.33
  c_4=4.48
  g_4=0.07
  
  f_4=a_4*g_4+(a_4-a_4*g_4)/(1+exp(-c_4-b_4*log(d)))
  

  
  # model_5- Quantal linear fit 
  
  a_5=0.03
  b_5=0.53
  
  f_5=a_5+(1-a_5)*(1-exp(-b_5*d))
  
  # model_6 Weibull fit
  
  a_6=0.03
  b_6=2.96
  c_6=0.83
  
  f_6=a_6+(1-a_6)*(1-exp(-c_6*d^b_6))

  # model_7 Log probit
  
  a_7=0.03
  b_7=2.06
  c_7=0.12
  
  f_7=a_7+(1-a_7)*pnorm((c_7+b_7*log(d)), mean=0, sd=1, lower.tail=T, log.p=F)
  
  
  # model_8 Multistage 2nd order 
  
  a_8=0.02
  b_8=0.12
  c_8=0.58
  
  f_8=a_8+(1-a_8)*(1-exp(-b_8*d-c_8*(d^2)))
  
  
  # 8 different models are considered here
  models<-rbind(f_1,f_2,f_3,f_4,f_5,f_6,f_7,f_8)
  models<-data.frame(models)
  return((t(models))%*%post_weight)
  
}


# Model average value is considered as true BMD
# Which one will be used as true? 

# The values of dose level 0 should be saved there


ma_0_to_75<-model_average(d_origin,post_weight,d_origin)




ma_0_to_75

plot(seq(0,75,by=1),ma_0_to_75)



plot_ma=data.frame(cbind(seq(0,75,by=1),ma_0_to_75))




# BMR =0.1
ggplot()+
  geom_point(aes(x=X1,y=X2),data=plot_ma)+geom_hline(aes(yintercept=0.02),color="red")+
  geom_hline(aes(yintercept=0.02+0.1),color="blue")+
  scale_x_continuous(limits=c(0,75),breaks=c(0,25,35.3758,50,75))+theme_bw()+labs(x="Dose level", y="Pi_hat_MA")+
  geom_vline(xintercept=35.3758,color="grey")


# first element-- non-treatment
pi_hat_ma_0<-ma_0_to_75[1]

# Return nearest two value of return 
risk_0_bmr<-pi_hat_ma_0+BMR


risk_0_bmr
# Nearest value 
ind_min_range<-which(abs(ma_0_to_75-risk_0_bmr)==min(abs(ma_0_to_75-risk_0_bmr)))

# Ma values nearest risk_0_bmr value
ma_0_to_75[ind_min_range]
ma_0_to_75[ind_min_range+1]

# Second round

d_origin[ind_min_range]

d_35_to_36<-seq(d_origin[ind_min_range]
  ,d_origin[ind_min_range]+1,0.01)



# Second grid for the value

ma_35_to_36<-model_average(d_35_to_36,post_weight,d_origin)

ind_min_range<-which(abs(ma_35_to_36-risk_0_bmr)==min(abs(ma_35_to_36-risk_0_bmr)))


#value of this


risk_0_bmr

#This is slight larger than origianl value



ma_35_to_36[ind_min_range]
ma_35_to_36[ind_min_range-1]


d_35_to_36[ind_min_range-1]
d_35_to_36[ind_min_range]


d_3537_3538<-seq(35.37,36.38,0.0001)

ma_3537_3538<-model_average(d_3537_3538,post_weight,d_origin)

ind_min_range<-which(abs(ma_3537_3538-risk_0_bmr)==min(abs(ma_3537_3538-risk_0_bmr)))

ind_min_range


# True BMD -- 35.3785 okay... BMD ture

# risk_0+0.1%
# Approximated around there true BMD

# The difference between the two values are almost close to 0
risk_0_bmr

ma_3537_3538[ind_min_range]


true_BMD<-d_3537_3538[ind_min_range]
true_BMD

# Let's say is true BMD from grid search is 35.3785

# Model fitting result analysis case

# true_BMD= 35.3785

library(readxl)
########################### decide true BMD from grid search#################################3

new_batch_sep_03_2019_0945_pm_bmrs <- read_excel("/Users/sooyeonglim/Downloads/bmd_fitting_result/bbmd_fit_result/new-batch-sep-03-2019-0945-pm-bmrs.xlsx")



bbmd_fit_result<-new_batch_sep_03_2019_0945_pm_bmrs
bbmd_fit_result


bbmd_fit_result_MA<-new_batch_sep_03_2019_0945_pm_bmrs

# Extract result of model average
bbmd_fit_result_MA<-bbmd_fit_result %>% 
  group_by(study_index) %>% 
  filter(model=='Model average')


bbmd_fit_result_MA_extra<-bbmd_fit_result_MA %>% filter(bmr %in% 'Extra') %>% mutate(true_BMD=true_BMD,
                                                                                               cover_sign=true_BMD-bmdl,
                                                                                               bias_sign=bmd-true_BMD,
                                                                                               sq_bias_sign=bias_sign^2)



bbmd_fit_result_MA_Added<-bbmd_fit_result_MA %>% filter(!bmr %in% 'Extra') %>% mutate(true_BMD=true_BMD,
                                                                                                cover_sign=true_BMD-bmdl,
                                                                                                bias_sign=bmd-true_BMD,
                                                                                                sq_bias_sign=bias_sign^2)




View(bbmd_fit_result_MA_extra)
View(bbmd_fit_result_MA_Added)

# a. Coverage : Pr(BMD>BMDL_hat)
# Coverage for each case


nrow(bbmd_fit_result_MA_extra %>% filter(cover_sign>0))
nrow(bbmd_fit_result_MA_Added %>% filter(cover_sign>0))


# coverage.. why
# BBMDS- 64% for extra
# BBMDS- 61% for added

# b. Half-Width CI: E[BMD-BMDL_hat]
# Confidence interval can should be... 

mean(bbmd_fit_result_MA_extra$cover_sign)
mean(bbmd_fit_result_MA_Added$cover_sign)

# c. Bias: B=E[BMD_hat-BMD]

mean(bbmd_fit_result_MA_Added$bias_sign)
mean(bbmd_fit_result_MA_extra$bias_sign)

# d. MSE: MSE=E[(BMD_hat-BMD)^2]=Var(BMD_hat)+B^2
mean(bbmd_fit_result_MA_Added$sq_bias_sign)
mean(bbmd_fit_result_MA_extra$sq_bias_sign)


# Caculated the value from BMDS

# Scketch of the summary statistics

# Assign & read multiple outputs all together


files <- list.files(path = "/Users/sooyeonglim/Downloads/09_29", pattern = "*.xlsx", full.names = T)
files



tbL_extra<- sapply(files[1:100], read_excel,sheet = "Summary", range = "F26:I27", 
              col_names = FALSE, simplify=FALSE) %>% bind_rows(.id = "id")

tbL_add <- sapply(files[1:100], read_excel,sheet = "Summary", range = "G41:I41", 
                  col_names = FALSE, simplify=FALSE) %>% bind_rows(.id = "id")



tbL_extra<-tbL_extra %>% filter(!is.na(...1))


tbL_add<-tbL_add %>% filter(!is.na(...1))

bmds_fit_extra_result<-tbL_extra %>% mutate(true_BMD=true_BMD,
                                      cover_sign=true_BMD-...3,
                                      bias_sign=...2-true_BMD,
                                      sq_bias_sign=bias_sign^2)



bmds_fit_add_result<-tbL_add %>% mutate(true_BMD=true_BMD,
                                            cover_sign=true_BMD-...2,
                                            bias_sign=...1-true_BMD,
                                            sq_bias_sign=bias_sign^2)



# a. Coverage : Pr(BMD>BMDL_hat)
# Make a table for this
# Coverage for each case



nrow(bmds_fit_extra_result %>% filter(cover_sign>0))
nrow(bmds_fit_add_result %>% filter(cover_sign>0))

# 100% for BMDS-Extra, 99% for BMDS-Added


# b. Half-Width CI: E[BMD-BMDL_hat]
# Confidence interval can should be... 

mean(bmds_fit_extra_result$cover_sign)
mean(bmds_fit_add_result$cover_sign)

# c. Bias: B=E[BMD_hat-BMD]

mean(bmds_fit_extra_result$bias_sign)
mean(bmds_fit_add_result$bias_sign)


# d. MSE: MSE=E[(BMD_hat-BMD)^2]=Var(BMD_hat)+B^2
mean(bmds_fit_extra_result$sq_bias_sign)
mean(bmds_fit_add_result$sq_bias_sign)


true_BMD
#Matt's package use-- debugging,
#table 100 is missing here
#BMDS should not have to have a gamma value here