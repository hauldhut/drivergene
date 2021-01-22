#library
if(!require(devtools)) install.packages("devtools")
devtools::install_github("huynguyen250896/geneSA", force=TRUE)
library(geneSA)

# create event vector for RNASeq data
#>median is up-regulated genes and <median is down regulated genes
event_rna <- apply(t(exp_dri_norm), 2, function(x) ifelse(x > median(x),"up","down")) %>% as.data.frame()
# check how many altered samples we have
table(as.matrix(event_rna))
# down    up 
# 29519 29505 

#create new column ‘status’ with survival event as binary
clinical_exp <- clinical_exp %>%
  mutate(status = ifelse(clinical_exp$OS_STATUS == "LIVING",0,1)) #set event: death = 1, alive = 0

#Make sure samples that in event_rna are also included in rows of clinical_exp and in exactly the same order
all(rownames(event_rna) == rownames(clinical_exp))
#[1] TRUE

#check which driver gene is significantly correlated with patient outcome
geneSA(data = event_rna, time = clinical_exp$OS_MONTHS, status = clinical_exp$status)

