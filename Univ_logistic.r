# Author:zwy
# Desc:WGCNA

library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(autoReg)
library(rrtable)
rm(list = ls())
getwd()


modProbes1 <- read.csv(file = 'output/WGCNA/Lipid_brown_module_genes.csv',row.names = 1)
modProbes2 <- read.csv(file = 'output/WGCNA/Lipid_turquoise_module_genes.csv',row.names = 1)
inter_lipid <- c(modProbes1$x,modProbes2$x)
exp_tpm <- read.csv(file = 'data/liposome/liposome_ori_data.csv', row.names = 1, check.names = F)
exp_tpm <- as.data.frame(t(exp_tpm))
gsg = goodSamplesGenes(exp_tpm, verbose = 3)
lipid_expr <- exp_tpm
lipid <- lipid_expr[,inter_lipid]
sample <- rownames(lipid)


phenotype <- read.csv('final_phenotype_data.csv', row.names = 1)
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


Uni_glm_f <- function(x){
  FML <- as.formula(paste0("gen_obe~", x)) 
  glm1 <- glm(FML, data=data, family = binomial) 
  glm2 <- summary(glm1)
  OR <- round(exp(coef(glm1)),4) 
  SE <- glm2$coefficients[,2] 
  CI5 <- round(exp(coef(glm1)-1.96*SE),4) 
  CI95 <- round(exp(coef(glm1)+1.96*SE),4) 
  #CI <- paste0(CI5,'-',CI95)
  P <- round(glm2$coefficients[,4],4) 
  Uni_glm <- data.frame('Characteristics'=x,
                        'OR' = OR,
                        #'CI' = CI,
                        'CI5'=CI5,
                        'CI95'=CI95,
                        'P' = P)[-1,]
  return(Uni_glm) 
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


data <- data %>% 
  dplyr::select(-slim,-BMI_biv_WST,-WHtR,-overweight,-ov_ob,-waistline)
colnames(data)[19:487] <- lipid_name_46
colnames(data)[19:487] <- gsub("[[:punct:]]", "_", colnames(data)[19:487])
Muti_uni <- glm(gen_obe ~.,
                data=data,
                family = binomial)
summary(Muti_uni)

ft <- autoReg(Muti_uni,uni = TRUE, multi = FALSE)%>%myft()
table2docx(ft,title="Univariate logistic regression for 469 lipids")


table(Uni_glm_result[,'P']<0.01)

Uni_glm_result[Uni_glm_result[,'P']<0.01,'Characteristics']
write.csv(Uni_glm_result[Uni_glm_result[,'P']<0.01,'Characteristics'],'./output/Lipid_lasso_input.csv')
