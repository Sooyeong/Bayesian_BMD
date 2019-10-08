### Generate Data for Simulation - BMDS, BMDL ### 
### Date - 06/20/19 ###
### Written by SL ###


# N_i is fixed as 50 
# Four Dose level - 0, 25, 50, 75

# Logistic model ~ Parameter estimates with smallest AIC 
rm(list=ls())
getwd()
setwd('/Users/sooyeonglim/Desktop/r_temp')        
# Set Seed
set.seed(1234)


d_origin


# d was normalized in shao's fitting example
d_origin<-c(0, 25, 50, 75)

d_norm=(d_origin-min(d_origin))/(max(d_origin)-min(d_origin))
d<-d_norm

# model_1 -Logistic
a_1<- -4.08
b_1<- 4.37
f_1<-1/(1+exp(-a_1-b_1*d))


# model_2 -Probit
a_2<--2.29
b_2<-2.44

f_2<-pnorm(a_2+b_2*d, mean=0, sd=1, lower.tail=T, log.p=F)

f_2

# model_3 -LogLogistic
a_3=0.03
b_3=3.49
c_3=0.2

f_3=a_3+(1-a_3)/(1+exp(-c_3-b_3*log(d)))

f_3
# model_4- Dichotomous Hill
a_4=0.54
b_4=9.33
c_4=4.48
g_4=0.07

f_4=a_4*g_4+(a_4-a_4*g_4)/(1+exp(-c_4-b_4*log(d)))

f_4


# model_5- Quantal linear fit 

a_5=0.03
b_5=0.53

f_5=a_5+(1-a_5)*(1-exp(-b_5*d))

f_5
# model_6 Weibull fit

a_6=0.03
b_6=2.96
c_6=0.83

f_6=a_6+(1-a_6)*(1-exp(-c_6*d^b_6))

f_6

# Log probit

a_7=0.03
b_7=2.06
c_7=0.12

f_7=a_7+(1-a_7)*pnorm((c_7+b_7*log(d)), mean=0, sd=1, lower.tail=T, log.p=F)

f_7



# Multistage 2nd order 

a_8=0.02
b_8=0.12
c_8=0.58

f_8=a_8+(1-a_8)*(1-exp(-b_8*d-c_8*(d^2)))

f_8

post_weight<-c(0.207,0.346,0.115,0.0468,0.00224,0.0953,0.125,0.0625)



# 8 different models are considered here

models
models<-rbind(f_1,f_2,f_3,f_4,f_5,f_6,f_7,f_8)
models<-data.frame(models)


MA_dose_0<-sum(models[,1]*post_weight)
MA_dose_25<-sum(models[,2]*post_weight)
MA_dose_50<-sum(models[,3]*post_weight)
MA_dose_75<-sum(models[,4]*post_weight)

MA_dose_0
MA_dose_25
MA_dose_50
MA_dose_75
# Propensity for each dose level - f_0, f_25, f_50, f_75 as propensity score
# Data generation for binomial distribution

# Create bench model
# N-- number of simulation
N<-100
sample_size<-50
number_dose=4
sim_data_06_25<-matrix(,N*number_dose,4)

for (i in 1:N){
        
        ## Generate new incidences for each simulation i 
        y_0<-rbinom(1,50,MA_dose_0)
        y_25<-rbinom(1,50,MA_dose_25)
        y_50<-rbinom(1,50,MA_dose_50)
        y_75<-rbinom(1,50,MA_dose_75)
        
        incidence<-c(y_0,y_25,y_50,y_75)
        
        for (j in 1:number_dose){
                
                sim_data_06_25[(i-1)*4+j,1]<-paste0("ID",i)
                sim_data_06_25[(i-1)*4+j,2]<-d_origin[j]
                sim_data_06_25[(i-1)*4+j,3]<-sample_size
                sim_data_06_25[(i-1)*4+j,4]<-incidence[j]
                
        }
}

col_names<-c('Dataset Index', 'Dose', 'N', 'Effect')
sim_data_06_25<-data.frame(sim_data_06_25)
colnames(sim_data_06_25)<-col_names

View(sim_data_06_25)
write.csv(sim_data_06_25,'sim_data_06_25.csv',row.names=F)
