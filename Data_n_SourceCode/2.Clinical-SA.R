#library
if(!require(survival)) install.packages('survival')
if(!require(Hmisc)) install.packages('Hmisc')
library('survival')
library('Hmisc')

# create event vector for RNASeq data
#>median is up-regulated genes and <median is down regulated genes
event_rna <- apply(exp_dri_norm,1, function(x) ifelse(x > median(x),"up","down"))
# check how many altered samples we have
table(event_rna)
# event_rna
# down    up 
# 33327 33313 
event_rna <- as.data.frame(event_rna) #should be as data frame

#create new column ‘status’ with survival event as binary
clinical_exp$status <- clinical_exp$OS_STATUS
clinical_exp$status <-ifelse(clinical_exp$status == "LIVING",0,1) #set event: die = 1, survive = 0
#add time and event columns of clinical_exp to event_rna
event_rna <- cbind(event_rna,clinical_exp[,3]) #time
event_rna <- cbind(event_rna,clinical_exp[,6]) #event
colnames(event_rna)[36:37] <- c("time", "event")
event_rna[1:5,1:5] 

#check what driver genes are significantly correlated with outcome
set.seed(420)
lapply(driver,
       
       function(x) {
         
         formula <- as.formula(paste('Surv(time,event)~',as.factor(x)))
         coxFit <- coxph(formula, data = event_rna)
         
         summary(coxFit)
       })