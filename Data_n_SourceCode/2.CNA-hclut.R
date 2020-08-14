#library
library(ComplexHeatmap) 
library(tidyverse)  # data manipulation
#######driver genes with CNV profile
cna_dri = c_cna
cna_dri = cna_dri[driver,]
dim(cna_dri) # 31 2173

#hierichical clustering
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
set.seed(258)
ac <- function(x) {
  agnes(t(cna_dri), method = x)$ac
}

map_dbl(m, ac)
# average    single  complete      ward 
# 0.7040206 0.6608054 0.7966347 0.9710952 

# Dissimilarity matrix
d <- dist(t(cna_dri), method = "euclidean")
# Ward's method
hc <- hclust(d, method = "ward.D2" )

#find the number of cluster
#library
dunn <- c("fpc","clValid","RankAggreg","kohonen","plyr","clv","cluster","stats","caret","party","partykit")
install.packages(dunn)
library("clValid")
#Dunn's index
transposed_cna=t(cna_dri); transposed_cna=as.data.frame(transposed_cna)
v <- clValid::clValid(transposed_cna, 2:15, clMethods="hierarchical",
                  validation="internal", metric = "euclidean", method = "ward")
#result
optimalScores(v)
# Score       Method Clusters
# Connectivity 263.1932540 hierarchical        2
# Dunn           0.1342312 hierarchical        2
# Silhouette     0.1382382 hierarchical        3
sizeGrWindow(4, 4)
plot(v)

# agnes() with cutree() cuts the dendrogram into 2 groups
hc_a <- agnes(t(cna_dri), method = "ward")
sub_grp= cutree(as.hclust(hc_a), k = 2)

# Number of members in each cluster
table(sub_grp)
# sub_grp
# 1    2 
# 993 1180 

## make a named vector from the vector
info =as.data.frame(sub_grp)
info$patient = rownames(info)
info <-info[order(info$sub_grp),]
info = dplyr::select(info,-patient)
colnames(info) <- c('groups')
info$groups = as.character(info$groups)
cna_dri = cna_dri[,rownames(info)] #change the order of column/patients of cna data following the variable 'info'

## Heatmap annotation
library(circlize)
ha <- columnAnnotation(df = info, col = list(groups = c("1" = "black", "2" = "green")))

#Plot heatmap
set.seed(1017)
Heatmap(cna_dri, name = "CNA scale", 
        show_row_names = TRUE, show_column_names = FALSE,
        cluster_columns = FALSE,show_column_dend = TRUE,
        show_row_dend = TRUE,top_annotation = ha,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        row_names_gp = gpar(fontsize = 8),
        column_order = colnames(cna_dri)
)
