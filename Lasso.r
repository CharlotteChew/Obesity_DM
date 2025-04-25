# Author:zwy
# Desc:Lasso reg

rm(list = ls())
getwd()
library(dplyr)
library(tibble)

## protein data
protein_expr <- read.csv(file = 'data/proteomics/NA_Processed1710167100.9991_MissForest.csv',row.names = 1)
protein_treated <- read.csv('data/proteomics/proteomics_treated_data.csv',row.names = 1)
protein_treated[,6:1971] <- NULL
protein_treated$sample <- rownames(protein_treated)
protein_matrix <- as.data.frame(t(protein_expr))
protein_matrix$sample <- rownames(protein_matrix)
protein_matrix <- protein_matrix %>% 
  dplyr::select(sample,everything())
colnames(protein_matrix)[1] <- 'code1'
protein <- merge(protein_treated,protein_matrix,by = 'code1')
protein <- protein %>% 
  column_to_rownames('sample')
protein <- protein %>% 
  dplyr::select(-c(1,2,3,4,5))
sample <- rownames(protein)
protein_lasso_input <- read.csv('output/uni/Protein_lasso_input.csv')
protein_lasso <- protein[,protein_lasso_input]

exp_tpm <- read.csv(file = 'data/liposome/liposome_ori_data.csv', row.names = 1, check.names = F)
exp_tpm <- as.data.frame(t(exp_tpm))
gsg = goodSamplesGenes(exp_tpm, verbose = 3)
lipid_expr <- exp_tpm
sample <- rownames(lipid)
colnames(lipid_expr) <- gsub("[[:punct:]]", "_", colnames(lipid_expr))
lipid_expr$sample <- rownames(lipid_expr)
Lipid_lasso_input <- read.csv('output/uni/Lipid_lasso_input.csv')
Lipid_lasso_input <- Lipid_lasso_input$x
lipid_lasso <- lipid_expr[,Lipid_lasso_input]

rownames(lipid_lasso) == rownames(protein_lasso)  
common_rows <- intersect(rownames(lipid_lasso), rownames(protein_lasso))
lipid_lasso_subset <- lipid_lasso[common_rows, ]
protein_lasso_subset <- protein_lasso[common_rows, ]
rownames(lipid_lasso_subset) == rownames(protein_lasso_subset)

lasso_input_expr <- cbind(protein_lasso_subset,lipid_lasso_subset)

phenotype <- read.csv('final_phenotype_data.csv',row.names = 1,check.names = F) 
phenotype <- phenotype[rownames(lasso_input_expr),]
table(phenotype$gen_obe)
phenotype$samples <- rownames(phenotype)
phenotype <- phenotype[,c('samples','gen_obe')]

lasso_input_expr$samples <- rownames(lasso_input_expr)
lasso_input_final <- merge(phenotype,lasso_input_expr,by = 'samples')
lasso_input_final <- lasso_input_final %>% 
  column_to_rownames('samples')
write.csv(lasso_input_final, file = 'output/lasso_input_final.csv')

library(glmnet)
set.seed(1234)
lasso_data <- lasso_input_final
x <- as.matrix(lasso_data[ , c(2:ncol(lasso_data))])
y <- as.matrix(lasso_data$gen_obe)
alpha1_fit <- glmnet(x, y, alpha = 1, family = "binomial", nlambda = 100)

png('./output/Lasso_regression_coef_path.png',res = 600,width = 6400,height = 4800)
plot(alpha1_fit, xvar = "lambda", label = TRUE)
dev.off()

set.seed(1234)
alpha1.fit.cv <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")

png('./output/Lasso_cross_valiedation.png',res = 600,width = 6400,height = 4800)
plot(alpha1.fit.cv)
dev.off()

print(alpha1.fit.cv)
coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)


feature_all <- as.data.frame(as.matrix(coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)))
colnames(feature_all) <- "coff"
feature_opt <-  feature_all %>% filter(abs(coff) > 0)
rownames(feature_opt)


write.csv(rownames(feature_opt)[-1],'output/Lasso_output.csv')