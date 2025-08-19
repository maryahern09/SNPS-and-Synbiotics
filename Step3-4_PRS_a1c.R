############################################################
################### SNPS Pipeline ##########################
#### Step 3 and 4: PRS Calc and Final Regression A1c########
############# created by Mary 7-28-2025 ####################
############################################################

#STEPS
#-1 Get data ready for PLINK
#0 Run SNP QC in PLINK 
#1.Data curation in R
#2 Regularization for coefficient reduction and feature selection using glmnet in R
#3 PRS calculation using PLINK through bigsnpr
#4 linear models using nlme or lme4 in R

#this was run on the HPC

#install packages
install.packages("bigsnpr")
install.packages("bigreadr")  # useful for reading GWAS files

library(bigsnpr)
library(bigreadr)

#load data
data_a1c = fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/Data_final_qc_a1c_ENfiltered.txt")
OUT_EN_a1c = fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/EN_a1c_OUT_filtered.txt")

#create several PRS, start at 20-70 by every 10
top_snps_a1c20 <- OUT_EN_a1c[order(abs(OUT_EN_a1c$Coef), decreasing = TRUE), ][1:20, ]
top_snps_a1c30 <- OUT_EN_a1c[order(abs(OUT_EN_a1c$Coef), decreasing = TRUE), ][1:30, ]
top_snps_a1c40 <- OUT_EN_a1c[order(abs(OUT_EN_a1c$Coef), decreasing = TRUE), ][1:40, ]
top_snps_a1c50 <- OUT_EN_a1c[order(abs(OUT_EN_a1c$Coef), decreasing = TRUE), ][1:50, ]
top_snps_a1c60 <- OUT_EN_a1c[order(abs(OUT_EN_a1c$Coef), decreasing = TRUE), ][1:60, ]
top_snps_a1c70 <- OUT_EN_a1c[order(abs(OUT_EN_a1c$Coef), decreasing = TRUE), ][1:70, ]

top_snps_a1c20 <- OUT_EN_a1c[order(OUT_EN_a1c$Coef), ][1:20, ]
top_snps_a1c30 <- OUT_EN_a1c[order(OUT_EN_a1c$Coef), ][1:30, ]
top_snps_a1c40 <- OUT_EN_a1c[order(OUT_EN_a1c$Coef), ][1:40, ]
top_snps_a1c50 <- OUT_EN_a1c[order(OUT_EN_a1c$Coef), ][1:50, ]
top_snps_a1c60 <- OUT_EN_a1c[order(OUT_EN_a1c$Coef), ][1:60, ]
top_snps_a1c70 <- OUT_EN_a1c[order(OUT_EN_a1c$Coef), ][1:70, ]

#Top 20
top_20 <- top_snps_a1c20$SNP
matrix_top20 <- as.matrix(data_a1c[,..top_20])
beta_coeff_top20 <-as.numeric(top_snps_a1c20$Coef)
PRS_top20 <- as.numeric(as.matrix(matrix_top20) %*% beta_coeff_top20)
data_a1c$PRS_top20 <- scale(PRS_top20) #scaling to make more interpretable

#Top 30
top_30 <- top_snps_a1c30$SNP
matrix_top30 <- as.matrix(data_a1c[,..top_30])
beta_coeff_top30 <-as.numeric(top_snps_a1c30$Coef)
PRS_top30 <- as.numeric(as.matrix(matrix_top30) %*% beta_coeff_top30)
data_a1c$PRS_top30 <- scale(PRS_top30)

#Top 40
top_40 <- top_snps_a1c40$SNP
matrix_top40 <- as.matrix(data_a1c[,..top_40])
beta_coeff_top40 <-as.numeric(top_snps_a1c40$Coef)
PRS_top40 <- as.numeric(as.matrix(matrix_top40) %*% beta_coeff_top40)
data_a1c$PRS_top40 <- scale(PRS_top40)

#Top 50
top_50 <- top_snps_a1c50$SNP
matrix_top50 <- as.matrix(data_a1c[,..top_50])
beta_coeff_top50 <-as.numeric(top_snps_a1c50$Coef)
PRS_top50 <- as.numeric(as.matrix(matrix_top50) %*% beta_coeff_top50)
data_a1c$PRS_top50 <- scale(PRS_top50)

#Top 60
top_60 <- top_snps_a1c60$SNP
matrix_top60 <- as.matrix(data_a1c[,..top_60])
beta_coeff_top60 <-as.numeric(top_snps_a1c60$Coef)
PRS_top60 <- as.numeric(as.matrix(matrix_top60) %*% beta_coeff_top60)
data_a1c$PRS_top60 <- scale(PRS_top60)

#Top 70
top_70 <- top_snps_a1c70$SNP
matrix_top70 <- as.matrix(data_a1c[,..top_70])
beta_coeff_top70 <-as.numeric(top_snps_a1c70$Coef)
PRS_top70 <- as.numeric(as.matrix(matrix_top70) %*% beta_coeff_top70)
data_a1c$PRS_top70 <- scale(PRS_top70)

#test association- 3 PC's (see which one performs the best)
model_a1c20 <- lm(hgb_a1c_glyco ~ PRS_top20 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 , data = data_a1c)
summary(model_glucose20)

model_a1c30 <- lm(hgb_a1c_glyco ~ PRS_top30 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 , data = data_a1c)
summary(model_a1c30)

model_a1c40 <- lm(hgb_a1c_glyco ~ PRS_top40 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 , data = data_a1c)
summary(model_a1c40)

model_a1c50 <- lm(hgb_a1c_glyco ~ PRS_top50 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 , data = data_a1c)
summary(model_a1c50)

model_a1c60 <- lm(hgb_a1c_glyco ~ PRS_top60 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 , data = data_a1c)
summary(model_a1c60)

model_a1c70 <- lm(hgb_a1c_glyco ~ PRS_top70 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 , data = data_a1c)
summary(model_a1c70)

#test association- 10 PC's
model_a1c20 <- lm(hgb_a1c_glyco ~ PRS_top20 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_a1c)
summary(model_glucose20)

model_a1c30 <- lm(hgb_a1c_glyco ~ PRS_top30 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_a1c)
summary(model_a1c30)

model_a1c40 <- lm(hgb_a1c_glyco ~ PRS_top40 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_a1c)
summary(model_a1c40)

model_a1c50 <- lm(hgb_a1c_glyco ~ PRS_top50 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_a1c)
summary(model_a1c50)

model_a1c60 <- lm(hgb_a1c_glyco ~ PRS_top60 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_a1c)
summary(model_a1c60)

model_a1c70 <- lm(hgb_a1c_glyco ~ PRS_top70 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_a1c)
summary(model_a1c70)


