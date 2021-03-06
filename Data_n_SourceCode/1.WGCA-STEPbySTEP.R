#library
install.packages(c("dynamicTreeCut","flashClust","Hmisc","WGCNA"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
library("dynamicTreeCut") # Module identification
library("flashClust")  #Fast implementation of hierarchical clustering
library("Hmisc")  #perform variables clustering
library("WGCNA")   #WGCNA tool
library("purrr")      #data processing
library("cluster")    # compute agglomerative coefficient

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads() ### Allowing parallel execution

#create clinical_exp1 for this step
clinical_exp1 = clinical_exp[,-c(3,5:6)] #remove OS_MONTHS, OS_STATUS and status columns

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(t(exp_dri_norm), powerVector = powers, verbose = 5,
                        networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.24,col="red")


softPower = 6;
adjacency = adjacency(t(exp_dri_norm), power = softPower,
                      type = "signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM
#hierichical clustering
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute agglomerative coefficient
set.seed(2583)
ac <- function(x) {
  agnes(exp_dri_norm, method = x)$ac
}

map_dbl(m, ac) # Agglomerative coefficient of each agglomeration method
# average    single  complete      ward 
# 0.4136737 0.3510251 0.4815638 0.5563910  

#assign gene names from adjacency to dissTOM
rownames(dissTOM) = rownames(adjacency) 
colnames(dissTOM) = colnames(adjacency) 
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "ward.D2");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,5)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     cex.lab = 1,cex.axis = 1, cex.main = 1.2);

#We set the minimum module size at 10:
minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# dynamicMods
# 1  2 
# 16 15 

# Convert numeric lables into colors
moduleColors = labels2colors(dynamicMods)
table(moduleColors)
# moduleColors
# blue turquoise 
# 15        16 
# Plot the dendrogram and colors underneath
sizeGrWindow(5,6)
plotDendroAndColors(geneTree, moduleColors, "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
