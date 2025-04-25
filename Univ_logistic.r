# 2025年3月24日
# Author:zwy
# Desc:WGCNA筛选出的关键脂质进行单因素逻辑回归筛选

library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(autoReg)
library(rrtable)
rm(list = ls())
getwd()

# 数据预处理
modProbes1 <- read.csv(file = 'output/WGCNA/Lipid_brown_module_genes.csv',row.names = 1)
modProbes2 <- read.csv(file = 'output/WGCNA/Lipid_turquoise_module_genes.csv',row.names = 1)
inter_lipid <- c(modProbes1$x,modProbes2$x)
exp_tpm <- read.csv(file = 'E:/ye/课题组/项目/诊断模型/DM/data/liposome/liposome_ori_data.csv', row.names = 1, check.names = F)
exp_tpm <- as.data.frame(t(exp_tpm))
gsg = goodSamplesGenes(exp_tpm, verbose = 3)
lipid_expr <- exp_tpm
lipid <- lipid_expr[,inter_lipid]
sample <- rownames(lipid)

### 表型数据处理
phenotype <- read.csv('E:/ye_r/gen_obe/final_phenotype_data.csv', row.names = 1) #读取填补并且选取之后的表型数据
phenotype <- phenotype %>% 
  select_if(~!any(is.na(.)))
colnames(phenotype)
# rownames(lipid_expr) == rownames(phenotype_male)
common_cols <- intersect(rownames(lipid_expr), rownames(phenotype))
lipid$sample <- rownames(lipid)
phenotype$sample <- rownames(phenotype)
data <- merge(phenotype, lipid, by = 'sample')
data <- data %>% 
  column_to_rownames('sample')

colnames(data)[36:504] <- gsub("[[:punct:]]", "_", colnames(data)[36:504])

# 单因素逻辑回归（脂质）
Uni_glm_f <- function(x){
  FML <- as.formula(paste0("gen_obe~", x)) #拟合结局和变量
  glm1 <- glm(FML, data=data, family = binomial) #glm()逻辑回归
  glm2 <- summary(glm1) #提取所有回归结果放入glm2中
  OR <- round(exp(coef(glm1)),4) #计算OR
  SE <- glm2$coefficients[,2] #提取SE
  CI5 <- round(exp(coef(glm1)-1.96*SE),4) # #计算95%CI保留两位小数并合并
  CI95 <- round(exp(coef(glm1)+1.96*SE),4) #计算95%CI保留两位小数并合并
  #CI <- paste0(CI5,'-',CI95)
  P <- round(glm2$coefficients[,4],4) #提取P值
  # 将变量名、OR、CI、P合并为一个表，删去第一行
  Uni_glm <- data.frame('Characteristics'=x,
                        'OR' = OR,
                        #'CI' = CI,
                        'CI5'=CI5,
                        'CI95'=CI95,
                        'P' = P)[-1,]
  return(Uni_glm) #返回循环函数继续上述操作 
} 

data <- data %>% 
  dplyr::select(-WHtR,)
Uni_variable <- colnames(data)[-11]
Uni_glm <- lapply(Uni_variable, Uni_glm_f)
library(plyr)
Uni_glm_result <- ldply(Uni_glm, data.frame)
Uni_glm_result
Uni_variable
Uni_glm_result <- Uni_glm_result[34:502,]
readr::write_excel_csv(Uni_glm_result,"output/Uni_glm_result(369lipid).csv")

## 绘制单因素回归表格
data <- data %>% 
  dplyr::select(-slim,-BMI_biv_WST,-WHtR,-overweight,-ov_ob,-waistline)
colnames(data)[19:487] <- lipid_name_46
colnames(data)[19:487] <- gsub("[[:punct:]]", "_", colnames(data)[19:487])
Muti_uni <- glm(gen_obe ~.,
                data=data,
                family = binomial)
summary(Muti_uni)
### 直接用Uni_glm_result绘制就可以了，用里面的数据
ft <- autoReg(Muti_uni,uni = TRUE, multi = FALSE)%>%myft()
table2docx(ft,title="Univariate logistic regression for 469 lipids")

# 看一下有多少脂质显著
table(Uni_glm_result[,'P']<0.01)
# 发现在p值小于0.01的条件下，有345个脂质是显著的，因为太多了无法画森林图展示

# 保存一下这些脂质的名称
Uni_glm_result[Uni_glm_result[,'P']<0.01,'Characteristics']
write.csv(Uni_glm_result[Uni_glm_result[,'P']<0.01,'Characteristics'],'./output/Lipid_lasso_input.csv')













################### 女
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(autoReg)
library(rrtable)
rm(list = ls())
getwd()

# 数据预处理
modProbes1 <- read.csv(file = 'output/WGCNA/lipid/Lipid_blue_module_genes_female.csv',row.names = 1)
modProbes2 <- read.csv(file = 'output/WGCNA/lipid/Lipid_turquoise_module_genes_female.csv',row.names = 1)
inter_lipid <- c(modProbes1$x,modProbes2$x)
lipid_expr <- read.csv(file = 'output/data/Lipid_exp_female.csv', row.names = 1, check.names = FALSE)  # check.names = FALSE 这个参数防止修改数据格式
lipid <- lipid_expr[,inter_lipid]
sample <- rownames(lipid)
### 表型数据处理
phenotype_female <- read.csv('output/data/Lipid_phenotype_female.csv',row.names = 1)
rownames(lipid_expr) == rownames(phenotype_female)
common_cols <- intersect(rownames(lipid_expr), rownames(phenotype_female))
lipid$sample <- rownames(lipid)
phenotype_female$sample <- rownames(phenotype_female)
data <- merge(phenotype_female, lipid, by = 'sample')
data <- data %>% 
  column_to_rownames('sample')

# 单因素逻辑回归（脂质）
Uni_glm_f <- function(x){
  FML <- as.formula(paste0("cen_obe~", x)) #拟合结局和变量
  glm1 <- glm(FML, data=data, family = binomial) #glm()逻辑回归
  glm2 <- summary(glm1) #提取所有回归结果放入glm2中
  OR <- round(exp(coef(glm1)),4) #计算OR
  SE <- glm2$coefficients[,2] #提取SE
  CI5 <- round(exp(coef(glm1)-1.96*SE),4) # #计算95%CI保留两位小数并合并
  CI95 <- round(exp(coef(glm1)+1.96*SE),4) #计算95%CI保留两位小数并合并
  #CI <- paste0(CI5,'-',CI95)
  P <- round(glm2$coefficients[,4],4) #提取P值
  # 将变量名、OR、CI、P合并为一个表，删去第一行
  Uni_glm <- data.frame('Characteristics'=x,
                        'OR' = OR,
                        #'CI' = CI,
                        'CI5'=CI5,
                        'CI95'=CI95,
                        'P' = P)[-1,]
  return(Uni_glm) #返回循环函数继续上述操作 
} 
colnames(data)[32:509] <- gsub("[[:punct:]]", "_", colnames(data)[32:509])
data <- data %>% 
  dplyr::select(-WHtR,)
Uni_variable <- colnames(data)[-10]
Uni_glm <- lapply(Uni_variable, Uni_glm_f)
library(plyr)
Uni_glm_result <- ldply(Uni_glm, data.frame)
Uni_glm_result
Uni_variable
Uni_glm_result <- Uni_glm_result[30:509,]
readr::write_excel_csv(Uni_glm_result,"output/female_Uni_glm_result(480lipid).csv")

## 绘制单因素回归表格
data <- data %>% 
  dplyr::select(-slim,-BMI_biv_WST,-WHtR,-overweight,-ov_ob,-waistline)
colnames(data)[19:487] <- lipid_name_46
colnames(data)[19:487] <- gsub("[[:punct:]]", "_", colnames(data)[19:487])
Muti_uni <- glm(cen_obe ~.,
                data=data,
                family = binomial)
summary(Muti_uni)
### 直接用Uni_glm_result绘制就可以了，用里面的数据
ft <- autoReg(Muti_uni,uni = TRUE, multi = FALSE)%>%myft()
table2docx(ft,title="Univariate logistic regression for 469 lipids")

# 看一下有多少脂质显著
table(Uni_glm_result[,'P']<0.01)
# 发现在p值小于0.01的条件下，有345个脂质是显著的，因为太多了无法画森林图展示

# 保存一下这些脂质的名称
Uni_glm_result[Uni_glm_result[,'P']<0.01,'Characteristics']
write.csv(Uni_glm_result[Uni_glm_result[,'P']<0.01,'Characteristics'],'./output/female_Lipid_lasso_input.csv')

