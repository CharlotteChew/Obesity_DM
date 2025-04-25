# Author:zwy
# Desc:WGCNA was performed on the lipidomic data to screen for key lipid

rm(list = ls())
library(WGCNA)
library(dplyr)

# lipidomic data
exp_tpm <- read.csv(file = 'data/liposome/liposome_ori_data.csv', row.names = 1)
exp_tpm <- as.data.frame(t(exp_tpm))
gsg = goodSamplesGenes(exp_tpm, verbose = 3)
datExpr <- exp_tpm

# phe
phenotype <- read.csv('final_phenotype_data.csv', row.names = 1) # fill data
phenotype <- phenotype %>% 
  select_if(~!any(is.na(.)))
colnames(phenotype)

sampleTree <- hclust(dist(exp_tpm), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
abline(h = 700000, col = "red")
clust <- cutreeStatic(sampleTree, cutHeight = 600000, minSize = 10)
table(clust)
keepSamples <- clust == 1
exp_tpm <- exp_tpm[keepSamples, ]
nGenes = ncol(exp_tpm)
nSamples = nrow(exp_tpm)
datExpr <- exp_tpm %>% filter(rownames(exp_tpm) %in% rownames(phenotype))
phenotype <- phenotype[rownames(datExpr),]

# choose soft threshold
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,RsquaredCut = 0.80)
sft$powerEstimate <- 12
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
power <- NA
if(is.na(power)){
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

# construct network
nGenes = ncol(datExpr)
cor <- WGCNA::cor
net <- blockwiseModules(datExpr, power = 10, 
                        maxBlockSize = nGenes, TOMType = "unsigned", 
                        minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = F, verbose = 3)
table(net$colors)
moduleColors <- labels2colors(net$colors)
table(moduleColors)

jpeg('output/WGCNA/Lipid_module_result.jpg', res = 600, height = 5600, width = 7000)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    main = 'lipidomic Cluster Dendrogram',
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLables <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
head(MEs)[1:5, 1:5]
geneTree <- net$dendrograms[[1]]
MEList <-  moduleEigengenes(datExpr, colors = moduleColors)
MEs0 <- MEList$eigengenes
head(MEs0)[1:5, 1:5]
MEs <- orderMEs(MEs0)
head(MEs)[1:5, 1:5]
moduleTraitCor <- cor(MEs, phenotype , use = "p");
head(moduleTraitCor)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)

jpeg('output/WGCNA/Lipid_module_relation_heatmap_phenotype.jpg', res = 600, height = 4800, width = 9600)
library(viridis)
colors = viridis(length(moduleTraitCor))
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(phenotype),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-0.4, 0.4), 
               main = paste("Relationships between lipids module and trait"))
dev.off()

## brown
module = "brown";
probes = colnames(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
head(modProbes)
write.csv(modProbes,"output/Lipid_brown_module_genes.csv")
modProbes <- read.csv(file = 'output/Lipid_brown_module_genes.csv',row.names = 1)
modProbes <- modProbes$x
## turquoise
module = "turquoise";
probes = colnames(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
head(modProbes)
write.csv(modProbes,"output/Lipid_turquoise_module_genes.csv")
modProbes <- read.csv(file = 'output/Lipid_turquoise_module_genes.csv',row.names = 1)
modProbes <- modProbes$x
