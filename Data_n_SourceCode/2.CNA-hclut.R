#library
library(ComplexHeatmap) 
install.packages("RColorBrewer")
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(dendextend) # for comparing two dendrograms
library(RColorBrewer)
#######heatmap.2
# the most 10 frequently amplified and 10 deleted genes of CNV

gene<-c("ERBB2", "GATA3", "PIK3CA", "RUNX1", "FOXO3", "KRAS", "NF1", "AKT1", "KMT2C", "CDKN1B",
        "BRCA2", "CBFB", "CDH1", "CDKN2A", "CHEK2", "ERBB3")
length(gene) #16
gene %in% rownames(c_cna) #check
cna_freq=c_cna
cna_freq = cna_freq[gene,]
dim(cna_freq) # 16 2173

#hierichical clustering
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
set.seed(258)
ac <- function(x) {
  agnes(t(cna_freq), method = x)$ac
}

map_dbl(m, ac)
# average    single  complete      ward 
# 0.7966950 0.7097392 0.8694155 0.9779233 
# Dissimilarity matrix
d <- dist(t(cna_freq), method = "euclidean")
# Ward's method
hc <- hclust(d, method = "ward.D2" )

#find the number of cluster
#library
pkgs <- c("factoextra",  "NbClust")
install.packages(pkgs)
library(factoextra)
library(NbClust)

#Silhouette
fviz_nbclust(t(cna_freq), hcut, method = "silhouette", k.max = 15)+
  labs(title = "Optimal number of groups", subtitle = "Silhouette method")+
  xlab("Number of groups k")#2 groups

# Cut tree into 2 groups
sub_grp <- cutree(hc, k = 2)

# Number of members in each cluster
table(sub_grp)
# sub_grp - maximum
# 1    2 
# 1733  440 

## make a named vector from the vector
info =as.data.frame(sub_grp)
colnames(info) <- c('groups')
info$groups = as.character(info$groups)
## Heatmap annotation
ha <- columnAnnotation(df = info)

Heatmap(cna_freq, name = "CNA scale", 
        show_row_names = TRUE, show_column_names = FALSE, 
        row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        top_annotation = ha
)