#import data
cna<-read.table('data_CNA.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = 1)

#check dimension and adjust slightly the form of CNV, clinical data.
cna=cna[,-1]
dim(cna) #22544  2173
#check missing value
cna=as.matrix(cna)
table(is.finite(cna))
# FALSE     TRUE 
# 2356 48985756
library(CancerSubtypes)
cna=data.imputation(cna, fun="mean")
#extract the same sample names CNV and clinical share
s_cna_cli=intersect(colnames(cna), clinical$PATIENT_ID)
length(s_cna_cli) #2173
c_idx_cna<-match(s_cna_cli,colnames(cna))
c_cna<-cna[,c_idx_cna]
c_idx_clinical<-match(s_cna_cli,clinical$PATIENT_ID)
c_clinical<-clinical[c_idx_clinical,]

dim(c_cna)
# 22544  2173
dim(c_clinical)
# 2173   22
