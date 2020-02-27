rm(list=ls())
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
#library
library(dplyr)
library(tidyr)

#import data
exp_dri_norm<-read.table('data_mRNA_median_Zscores.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL)
clinical<-read.table('data_clinical_patient.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL,fill=TRUE)
#check dimension
dim(clinical) # 2509   22
dim(exp_dri_norm) # 18534  1906
#####prepare data
driver=c("MAP2K4","ARID1A", "PIK3CA", "TBX3", "MAP3K1", "CBFB", "TP53", "KMT2C",
         "AKT1", "GATA3", "RUNX1", "PTEN", "CDH1", "NF1", "PIK3R1", "RB1", "CDKN1B",
         "NCOR1", "CDKN2A", "ERBB2", "ERBB3", "FOXO3", "SMAD4", "KRAS", "BRCA2",
         "BAP1", "GPS2", "AGTR2", "ZFP36L1", "MEN1", "CHEK2", "SF3B1", "AHNAK2",
         "SYNE1", "MUC16")
length(driver) #35
table(driver %in% exp_dri_norm$Hugo_Symbol) #check, TRUE = 35
table(exp_dri_norm$Hugo_Symbol %in% driver) #check
#exp with row = driver genes
exp_dri_norm=exp_dri_norm[,-2]#remove Entrez_Gene_Id column
exp_dri_norm = exp_dri_norm %>% filter (exp_dri_norm$Hugo_Symbol %in% driver)
rownames(exp_dri_norm) = exp_dri_norm$Hugo_Symbol
exp_dri_norm=exp_dri_norm[,-1] #remove Hugo_Symbol column

exp_dri_norm = as.matrix(exp_dri_norm)
#intersect of samples between RNA_seq and clinical_exp
s_exp_cli=intersect(colnames(exp_dri_norm), clinical$PATIENT_ID)
length(s_exp_cli) #1904
c_idx_exp1<-match(s_exp_cli,colnames(exp_dri_norm))
exp_dri_norm = exp_dri_norm[,c_idx_exp1]

c_idx_cli7<-match(s_exp_cli,clinical$PATIENT_ID)
clinical_exp = clinical[c_idx_cli7,]

# keep columns having clinical features of interest
rownames(clinical_exp) = clinical_exp$PATIENT_ID
clinical_exp=data.frame(clinical_exp$LYMPH_NODES_EXAMINED_POSITIVE, clinical_exp$NPI,
                     clinical_exp$OS_MONTHS, clinical_exp$stage, clinical_exp$OS_STATUS)
colnames(clinical_exp) = c("lymph", "npi", "OS_MONTHS", "stage", "OS_STATUS")
rownames(clinical_exp) = colnames(exp_dri_norm)

dim(clinical_exp) #1904 5
dim(exp_dri_norm)# 35 1904
max(exp_dri_norm)# 13.1189
min(exp_dri_norm)# -7.5848
table(is.finite(as.matrix(exp_dri_norm))) #check missing value

exp_dri_norm[1:5,1:5]
clinical_exp[1:5,1:5]

