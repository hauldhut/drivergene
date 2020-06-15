#library
if(!require(devtools)) install.packages("devtools")
devtools::install_github("huynguyen250896/computeC")
library(computeC) #compute the correlation

#create the necessary df
cor=t(exp_dri_norm)
cor <- as.data.frame(cor) #should be as data frame

feature = clinical_exp[,c(1,2,4)]
colnames(feature)[1:3] <- c("lymph", "npi", "stage") #change the columns’ name
feature <- as.data.frame(feature) #should be as data frame

set.seed(25896)
#####correlation between driver genes and lymph
computeC(cor,feature,"lymph")
#####correlation between driver genes and npi
computeC(cor,feature,"npi")
#####correlation between driver genes and stage
computeC(cor,feature,"stage")

