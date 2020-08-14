#library
if(!require(devtools)) install.packages("devtools")
devtools::install_github("huynguyen250896/geneSA")
library(geneSA)
if(!require(rlist)) install.packages("https://cran.r-project.org/src/contrib/Archive/rlist/rlist_0.4.tar.gz", repos = NULL); library(rlist)

# create event vector for RNASeq data
#>median is up-regulated genes and <median is down regulated genes
event_rna <- apply(exp_dri_norm,1, function(x) ifelse(x > median(x),"up","down"))
# check how many altered samples we have
table(event_rna)
# event_rna
# down    up 
# 29519 29505 
event_rna <- as.data.frame(event_rna) #should be as data frame

#create new column ‘status’ with survival event as binary
clinical_exp$status <- clinical_exp$OS_STATUS
clinical_exp$status <-ifelse(clinical_exp$status == "LIVING",0,1) #set event: die = 1, survive = 0
#add time and event columns of clinical_exp to event_rna
event_rna <- cbind(event_rna,clinical_exp[,3]) #time
event_rna <- cbind(event_rna,clinical_exp[,6]) #event
colnames(event_rna)[32:33] <- c("time", "event")

#check what driver genes are significantly correlated with outcome
geneSA(genename = driver, event = event_rna)