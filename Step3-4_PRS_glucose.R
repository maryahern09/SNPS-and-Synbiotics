####################################################################
######################## SNPS Pipeline #############################
#### Step 3 and 4: PRS Calc and Final Regression Glucose ###########
################# created by Mary 7-28-2025 ########################
####################################################################

#STEPS
#-1 Get data ready for PLINK
#0 Run SNP QC in PLINK 
#1.Data curation in R
#2 Regularization for coefficient reduction and feature selection using glmnet in R
#3 PRS calculation using PLINK through bigsnpr
#4 linear models using nlme or lme4 in R

#this was run on the HPC

#install packages ####

#install packages
library(data.table)

#load data
data_glucose = fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/Data_final_qc_glucose_ENfiltered.txt")
OUT_EN_glucose = fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/EN_Glucose_OUT_filtered.txt")

#create several PRS, start at 20-70 by every 10
top_snps_glucose20 <- OUT_EN_glucose[order(abs(OUT_EN_glucose$Coef), decreasing = TRUE), ][1:20, ]
top_snps_glucose30 <- OUT_EN_glucose[order(abs(OUT_EN_glucose$Coef), decreasing = TRUE), ][1:30, ]
top_snps_glucose40 <- OUT_EN_glucose[order(abs(OUT_EN_glucose$Coef), decreasing = TRUE), ][1:40, ]
top_snps_glucose50 <- OUT_EN_glucose[order(abs(OUT_EN_glucose$Coef), decreasing = TRUE), ][1:50, ]
top_snps_glucose60 <- OUT_EN_glucose[order(abs(OUT_EN_glucose$Coef), decreasing = TRUE), ][1:60, ]
top_snps_glucose70 <- OUT_EN_glucose[order(abs(OUT_EN_glucose$Coef), decreasing = TRUE), ][1:70, ]

#Top 20
top_20 <- top_snps_glucose20$SNP
matrix_top20 <- as.matrix(data_glucose[,..top_20])
beta_coeff_top20 <-as.numeric(top_snps_glucose20$Coef)
PRS_top20 <- matrix_top20 %*% beta_coeff_top20
data_glucose$PRS_top20 <- scale(PRS_top20) #scaling to make more interpretable
hist(data_glucose$PRS_top70)


#Top 30
top_30 <- top_snps_glucose30$SNP
matrix_top30 <- as.matrix(data_glucose[,..top_30])
beta_coeff_top30 <-as.numeric(top_snps_glucose30$Coef)
PRS_top30 <- matrix_top30 %*% beta_coeff_top30
data_glucose$PRS_top30 <- scale(PRS_top30)

#Top 40
top_40 <- top_snps_glucose40$SNP
matrix_top40 <- as.matrix(data_glucose[,..top_40])
beta_coeff_top40 <-as.numeric(top_snps_glucose40$Coef)
PRS_top40 <- matrix_top40 %*% beta_coeff_top40
data_glucose$PRS_top40 <- scale(PRS_top40)

#Top 50
top_50 <- top_snps_glucose50$SNP
matrix_top50 <- as.matrix(data_glucose[,..top_50])
beta_coeff_top50 <-as.numeric(top_snps_glucose50$Coef)
PRS_top50 <- matrix_top50 %*% beta_coeff_top50
data_glucose$PRS_top50 <- scale(PRS_top50)

#Top 60
top_60 <- top_snps_glucose60$SNP
matrix_top60 <- as.matrix(data_glucose[,..top_60])
beta_coeff_top60 <-as.numeric(top_snps_glucose60$Coef)
PRS_top60 <- matrix_top60 %*% beta_coeff_top60
data_glucose$PRS_top60 <- scale(PRS_top60)

#Top 70
top_70 <- top_snps_glucose70$SNP
matrix_top70 <- as.matrix(data_glucose[,..top_70])
beta_coeff_top70 <-as.numeric(top_snps_glucose70$Coef)
PRS_top70 <- matrix_top70 %*% beta_coeff_top70
data_glucose$PRS_top70 <- scale(PRS_top70)

#test association
#null model for comparison
null_model <- lm(glucose ~ gender_2 + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
anova(null_model, model_glucose20)
summary(null_model)

model_glucose20 <- lm(glucose ~ PRS_top20 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
summary(model_glucose20)

model_glucose30 <- lm(glucose ~ PRS_top30 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
summary(model_glucose30)

model_glucose40 <- lm(glucose ~ PRS_top40 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
summary(model_glucose40)

model_glucose50 <- lm(glucose ~ PRS_top50 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
summary(model_glucose50)

model_glucose60 <- lm(glucose ~ PRS_top60 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
summary(model_glucose60)

model_glucose70 <- lm(glucose ~ PRS_top70 + as.factor(gender_2) + age_2 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10, data = data_glucose)
summary(model_glucose70)

