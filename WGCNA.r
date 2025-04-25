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























# 根据性别分组
sex_group_1 <- phenotype %>% filter(Sex == 1)  # 男性
sex_group_2 <- phenotype %>% filter(Sex == 2)  # 女性
# 获取对应性别的 protein 数据
datExpr_male <- datExpr[rownames(datExpr) %in% rownames(sex_group_1), ]
datExpr_female <- datExpr[rownames(datExpr) %in% rownames(sex_group_2), ]
# 分别选择对应性别的表型数据
phenotype_male <- phenotype %>% filter(rownames(phenotype) %in% rownames(datExpr_male))
phenotype_female <- phenotype %>% filter(rownames(phenotype) %in% rownames(datExpr_female))
# 使用 dplyr 的 select 函数去除 'Sex' 列
phenotype_male <- phenotype_male[, !colnames(phenotype_male) %in% "Sex"]
phenotype_female <- phenotype_female[, !colnames(phenotype_female) %in% "Sex"]
# 确保行名匹配  
datExpr_male <- datExpr_male[rownames(datExpr_male) %in% rownames(phenotype_male), ]
datExpr_female <- datExpr_female[rownames(datExpr_female) %in% rownames(phenotype_female), ]
phenotype_male <- phenotype_male[rownames(datExpr_male), ]
phenotype_female <- phenotype_female[rownames(datExpr_female), ]

# 保存结果
write.csv(datExpr_male, "E:/ye/课题组/项目/诊断模型/DM/output/Lipid_exp_male.csv")
write.csv(datExpr_female, "E:/ye/课题组/项目/诊断模型/DM/output/Lipid_exp_female.csv")
write.csv(phenotype_male, "E:/ye/课题组/项目/诊断模型/DM/output/Lipid_phenotype_male.csv")
write.csv(phenotype_female, "E:/ye/课题组/项目/诊断模型/DM/output/Lipid_phenotype_female.csv")

exp_tpm_male <- datExpr_male
exp_tpm_female <- datExpr_female

# 根据性别分别进行树形聚类和其他分析
sampleTree_male <- hclust(dist(exp_tpm_male), method = "average")
sampleTree_female <- hclust(dist(exp_tpm_female), method = "average")
# 为男性和女性样本分别绘制聚类图
plot(sampleTree_male, main = "Male Sample clustering to detect outliers", sub = "", xlab = "")
abline(h = 600000, col = "red")
plot(sampleTree_female, main = "Female Sample clustering to detect outliers", sub = "", xlab = "")
abline(h = 600000, col = "red")

# 你可以继续为性别分组做进一步分析
clust_male <- cutreeStatic(sampleTree_male, cutHeight = 600000, minSize = 10)
clust_female <- cutreeStatic(sampleTree_female, cutHeight = 600000, minSize = 10)
table(clust_male)
table(clust_female)
keepSamples_male <- clust_male == 1
keepSamples_female <- clust_female == 1
exp_tpm_male <- exp_tpm_male[keepSamples_male, ]
exp_tpm_female <- exp_tpm_female[keepSamples_female, ]
nGenes_male = ncol(exp_tpm_male)
nSamples_male = nrow(exp_tpm_male)
nGenes_female = ncol(exp_tpm_female)
nSamples_female = nrow(exp_tpm_female)
phenotype_male <- phenotype_male[rownames(exp_tpm_male), ]
phenotype_female <- phenotype_female[rownames(exp_tpm_female), ]





# 两个子集：exp_tpm_male, phenotype_male 和 exp_tpm_female, phenotype_female 来进行后续建模。

# 先对male子集分析
# 选择最佳阈值
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(exp_tpm_male, powerVector = powers, verbose = 5,RsquaredCut = 0.80)
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
  nSamples=nrow(exp_tpm_male)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

# 构建网络
nGenes = ncol(exp_tpm_male)
cor <- WGCNA::cor
net <- blockwiseModules(exp_tpm_male, power = 10, 
                        maxBlockSize = nGenes, TOMType = "unsigned", 
                        minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = F, verbose = 3)
table(net$colors)
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# png('output/Lipid_module_result.png',res = 600,,height = 5600,width = 7000)
jpeg('E:/ye/课题组/项目/诊断模型/DM/output/Lipid_module_result_male.jpg', res = 600, height = 5600, width = 7000)
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
MEList <-  moduleEigengenes(exp_tpm_male, colors = moduleColors)
MEs0 <- MEList$eigengenes
head(MEs0)[1:5, 1:5]
MEs <- orderMEs(MEs0)
head(MEs)[1:5, 1:5]
moduleTraitCor <- cor(MEs, phenotype_male , use = "p");
head(moduleTraitCor)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
# png('output/Lipid_module_relation_heatmap_phenotype.png',res = 600, height = 4800,width = 9600)
jpeg('E:/ye/课题组/项目/诊断模型/DM/output/Lipid_module_relation_heatmap_phenotype_male.jpg', res = 600, height = 4800, width = 9600)
library(viridis)
colors = viridis(length(moduleTraitCor))
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(phenotype_male),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               # 这里大家要注意，一般会设置为zlim = c(-1, 1)，怪我怪我选的这个例子不好，性状和模块之间相关性太低啦！
               # 不过如果你的相关性也很低，咱们也可以适当调整哟！
               zlim = c(-0.4, 0.4), 
               main = paste("Relationships between lipids module and trait"))
dev.off()

# 提取green和turquoise模块内的hub genes并且保存
## 提取green
module = "green";
probes = colnames(exp_tpm_male)
inModule = (moduleColors==module);
modProbes = probes[inModule];
head(modProbes)##之后就可以GO分析
write.csv(modProbes,"E:/ye/课题组/项目/诊断模型/DM/output/Lipid_green_module_genes_male.csv")
modProbes <- read.csv(file = 'E:/ye/课题组/项目/诊断模型/DM/output/Lipid_green_module_genes_male.csv',row.names = 1)
modProbes <- modProbes$x
## 提取turquoise
module = "turquoise";
probes = colnames(exp_tpm_male)
inModule = (moduleColors==module);
modProbes = probes[inModule];
head(modProbes)##之后就可以GO分析
write.csv(modProbes,"E:/ye/课题组/项目/诊断模型/DM/output/Lipid_turquoise_module_genes_male.csv")
modProbes <- read.csv(file = 'E:/ye/课题组/项目/诊断模型/DM/output/Lipid_turquoise_module_genes_male.csv',row.names = 1)
modProbes <- modProbes$x



# 再对female子集分析
# 选择最佳阈值
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(exp_tpm_female, powerVector = powers, verbose = 5,RsquaredCut = 0.80)
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
  nSamples=nrow(exp_tpm_female)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

# 构建网络
nGenes = ncol(exp_tpm_female)
cor <- WGCNA::cor
net <- blockwiseModules(exp_tpm_female, power = 10, 
                        maxBlockSize = nGenes, TOMType = "unsigned", 
                        minModuleSize = 20, reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = F, verbose = 3)
table(net$colors)
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# png('output/Lipid_module_result.png',res = 600,,height = 5600,width = 7000)
jpeg('E:/ye/课题组/项目/诊断模型/DM/output/Lipid_module_result_female.jpg', res = 600, height = 5600, width = 7000)
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
MEList <-  moduleEigengenes(exp_tpm_female, colors = moduleColors)
MEs0 <- MEList$eigengenes
head(MEs0)[1:5, 1:5]
MEs <- orderMEs(MEs0)
head(MEs)[1:5, 1:5]
moduleTraitCor <- cor(MEs, phenotype_female , use = "p");
head(moduleTraitCor)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
# png('output/Lipid_module_relation_heatmap_phenotype.png',res = 600, height = 4800,width = 9600)
jpeg('E:/ye/课题组/项目/诊断模型/DM/output/Lipid_module_relation_heatmap_phenotype_female.jpg', res = 600, height = 4800, width = 9600)
library(viridis)
colors = viridis(length(moduleTraitCor))
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(phenotype_female),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colors,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               # 这里大家要注意，一般会设置为zlim = c(-1, 1)，怪我怪我选的这个例子不好，性状和模块之间相关性太低啦！
               # 不过如果你的相关性也很低，咱们也可以适当调整哟！
               zlim = c(-0.4, 0.4), 
               main = paste("Relationships between lipids module and trait"))
dev.off()

# 提取blue和turquoise模块内的hub genes并且保存
## 提取blue
module = "blue";
probes = colnames(exp_tpm_male)
inModule = (moduleColors==module);
modProbes = probes[inModule];
head(modProbes)##之后就可以GO分析
write.csv(modProbes,"E:/ye/课题组/项目/诊断模型/DM/output/Lipid_blue_module_genes_female.csv")
modProbes <- read.csv(file = 'E:/ye/课题组/项目/诊断模型/DM/output/Lipid_blue_module_genes_female.csv',row.names = 1)
modProbes <- modProbes$x
## 提取turquoise
module = "turquoise";
probes = colnames(exp_tpm_male)
inModule = (moduleColors==module);
modProbes = probes[inModule];
head(modProbes)##之后就可以GO分析
write.csv(modProbes,"E:/ye/课题组/项目/诊断模型/DM/output/Lipid_turquoise_module_genes_female.csv")
modProbes <- read.csv(file = 'E:/ye/课题组/项目/诊断模型/DM/output/Lipid_turquoise_module_genes_female.csv',row.names = 1)
modProbes <- modProbes$x

