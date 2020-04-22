#############################################################################

# Date: 04/21/2020
# Program name: qard_4to6.R
# Programmer: Sooyeong Lim
# Purpose of the program: Extract QRAD database with dose levels with 4-6
# Note: ID need to be matched later to specify chemical name and other informations (e.g., source of data)

#############################################################################

# Load libraries
library(dplyr)
library(ggplot2)

# Read QRAD database
load(url("http://www.users.miamioh.edu/baileraj/research/database.RData"))

# The data is saved at final.data
# The data base has 11 columns, including data source and normalized dose level;

names(final.data)
ncol(final.data)

dose_groups<-final.data %>% 
  group_by(ID) %>%
  summarize(nDoseGrps=n())

# There are 383 dataset which has dose levels between 4~6
ID_ov4to6<-dose_groups %>% 
  filter(nDoseGrps>=4, nDoseGrps<=6)

Dose_4to6<-final.data %>% 
  filter(ID %in% ID_ov4to6$ID)




# Column names need to be changed
# ID->Dataset Index, dose->Dose, n->N, obs->Effect

selected<-Dose_4to6 %>%
  select(ID,dose,n,obs)%>%
  rename("Dataset Index"=ID, "Dose"=dose, "N"=n, "Effect"=obs)



selected[,1]<-paste0("ID",selected[,1],sep="")

write.csv(selected,"dose_4to6.csv",row.names=FALSE,fileEncoding = "UTF-8")
