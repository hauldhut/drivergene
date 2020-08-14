#Define numbers of genes and samples
nGenes = ncol(t(exp_dri_norm));
nSamples = nrow(t(exp_dri_norm));
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(exp_dri_norm), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, clinical_exp1, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(20,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
colnames(clinical_exp1)[1] <-"lymph"
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical_exp1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-clinical feature relationships"))

####blue
# names (colors) of the modules
# Define variable stage containing the stage column of exp_dri_norm
stage = as.data.frame(clinical_exp1$stage);
names(stage) = "stage"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(t(exp_dri_norm), MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(t(exp_dri_norm), stage, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(stage), sep="");
names(GSPvalue) = paste("p.GS.", names(stage), sep="");

#Intramodular analysis: identifying genes with high GS and MM
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
x <- abs(geneModuleMembership[moduleGenes, column])
y <- abs(geneTraitSignificance[moduleGenes, 1])
limit <- range(c(x,y)) 
r <- round(cor(x,y),2)
plot(    x, 
         y, 
         ylim =range(0.00,0.25), xlim =range(0.1,0.7), 
         xlab =paste("Module Membership in", module, "module"), 
         ylab ="Gene significance for stage", col = module)
fit <-lm(y~x)
pintra= round(summary(fit)$coefficients[,4][[2]],2) #p-value
abline(fit, col='orange')
legend('topleft', col = "orange", lty = 1, box.lty = 1 ,legend = 'regression line', cex= 0.8) 
mtext(paste('correlation = ', r, ",", "P-value = ", pintra))
title(paste("Module membership vs. gene significance\n"))
text(y~x, labels=rownames(exp_dri_norm)[moduleColors=="blue"],data=exp_dri_norm, cex=0.8, font=4, pos = 2) #add label

#identifying genes with high GS and MM
intra_modular_analysis=data.frame(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]))
rownames(intra_modular_analysis) = colnames(t(exp_dri_norm))[moduleColors=="blue"] #only the blue module
View(intra_modular_analysis)

#high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the turquoise module to all other genes in the turquoise)
connectivity=intramodularConnectivity(adjacency, moduleColors)
connectivity = connectivity[colnames(t(exp_dri_norm))[moduleColors=="blue"],] #only the turquoise module
order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
connectivity = connectivity[order.kWithin,] #order rows following kWithin
connectivity = connectivity[1:5,] #top 5 genes that have a high connectivity to other genes in the turquoise module
View(connectivity)
