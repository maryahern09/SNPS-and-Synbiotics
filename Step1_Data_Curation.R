############################################
############ SNPS Pipeline #################
###### Step 1: Data Curation in R ##########
###### created by Mary 8-28-2025 ###########
############################################

#STEPS
#-1 Get data ready for PLINK
#0 Run SNP QC in PLINK 
#1.Data curation in R and Creation of final cleaned dataset for glycogenes analysis
#2 Regularization for coefficient reduction and feature selection using glmnet in R
#3 PRS calculation using PLINK through bigsnpr
#4 linear models using nlme or lme4 in R


#install packages
library(data.table)
library(readxl)
library(tidyverse)

setwd("C:/Users/marya/OneDrive/Documents/Research/Projects/SNPs and synbiotics/Analysis/R")

#load data
data_SNPs = fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/GlycogenesSNPsOUT.raw")
data_phenotype = read_excel("~/Research/Projects/SNPs and synbiotics/Analysis/R code/Data Phenotype for Transfer 1.15.25-Sri.xlsx")
data_pca <- fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/PCA_wBasicQC (2).eigenvec", header = FALSE)
  setnames(data_pca, c("FID", "IID", paste0("PCA", 1:(ncol(data_pca) - 2))))
pruned_SNPs <- read.table("~/Research/Projects/SNPs and synbiotics/Analysis/R code/plink.prune.in", header = FALSE)
  
#prep to merge datasets by IID
data_SNPs2 <- data_SNPs %>%
  select(-c("FID", "MAT", "PAT", "PHENOTYPE", "SEX"))%>%
  rename_with(~ str_remove(.x, "_[^_]+$"))

# Subset list of SNPs names that made it through PLINK QC
pruned_SNPs <- pruned_SNPs$V1

data_SNPs3 <- data_SNPs2 %>%
  select(1, all_of(pruned_SNPs))

data_phenotype2 <- data_phenotype %>%
  rename(IID = biobank_id)%>%
  rename(gender_2 = "gender_2 (Female is 0)")%>%
  select(-c("record_id", "AT_number", "Height (inches)"))

data_pca2 <- data_pca %>%
  select(-"FID")

#merge datasets
data_merged <- merge(data_phenotype2, data_pca2, by = "IID")
final_dataset_merged <- merge(data_merged, data_SNPs3, by = "IID") #dataset now has 1026 individuals and 7577 variables

### Data/phenotype Quality control 
#make sure important numeric columns are read as numeric
final_dataset_merged$glucose <- as.numeric(final_dataset_merged$glucose)
final_dataset_merged$hgb_a1c_glyco <- as.numeric(final_dataset_merged$hgb_a1c_glyco)

#make sure there are no missing values in variables that will be in the models
#glucose and hba1c will be done later for their respective datasets
#age and gender will be done now
#confirmed there is no missing for age using JMP, but 6 missing for gender

summary(final_dataset_merged$age_2) #confirmed no NA
table(final_dataset_merged$gender_2) #confirmed 6

final_dataset_merged$gender_2 <- as.numeric(final_dataset_merged$gender_2)

data_final <- final_dataset_merged %>%
  drop_na(gender_2)

data_final$gender_2 <- as.factor(data_final$gender_2)

#sanity check
summary(data_final$age_2) #confirmed no NA
table(data_final$gender_2) #confirmed no NA

####data_final is 1020 observations and 7,577 variables

#remove outliers above 95% and below 5%
lower <- quantile(data_final$glucose, 0.05, na.rm = TRUE) #lower is 82
upper <- quantile(data_final$glucose, 0.95, na.rm = TRUE) #upper is 305

data_filtered <- data_final %>%
  filter(glucose >82 & glucose <305)

###final cleaned dataset (data_filtered)
#877 observations
#7577 variables

##### Now we will do complete cases on A1c and glucose individually and create binary variables for logistic regression ####
# filter out rows where hgb_a1c_glyco is NA for complete cases
summary(data_filtered$hgb_a1c_glyco) #confirmed 0 NA

#filter out rows where glucose is NA
summary(data_filtered$glucose) #confirmed 0 NA

# There are no missing in both A1c and glucose so we will add the new variables and only export 1 dataset
# A1c variable creation

#create DM categories
data_hbA1c <- data_filtered %>%
  mutate(DM = case_when(
    hgb_a1c_glyco <5.7 ~ "Normal",
    hgb_a1c_glyco >=5.7 & hgb_a1c_glyco <=6.4 ~ "Pre-DM",
    hgb_a1c_glyco >=6.5 ~ "Diabetes"
    ))

#check numbers/ sanity check
table(data_hbA1c$DM) #Normal- 179 , PreDM- 192 , DM- 506, no NA
class(data_hbA1c$DM) #character

#create binary variable
data_hbA1c <- data_hbA1c %>%
  mutate(DM_binary = case_when(
    DM %in% c("Normal","Pre-DM") ~ 0,
    DM == "Diabetes" ~ 1
  ))

#check numbers/ sanity check 
table(data_hbA1c$DM_binary) #0- 371, 1- 506, no NA
class(data_hbA1c$DM_binary) #numeric
data_hbA1c$DM_binary <- as.factor(data_hbA1c$DM_binary)
class(data_hbA1c$DM_binary) #factor

# Glucose variable creation
#create DM categories
data_glucose <- data_hbA1c %>%
  mutate(glucose_cat = case_when(
    glucose < 100 ~ "Normal",
    glucose >= 100 & glucose <= 125 ~ "Pre-DM",
    glucose >= 126 ~ "Diabetes"
    ))

#check numbers
table(data_glucose$glucose_cat)#Normal- 290 , PreDM- 208 , DM- 379
class(data_glucose$glucose_cat)#character

#create binary variable
data_glucose <- data_glucose %>%
  mutate(glucose_binary = case_when(
    glucose_cat %in% c("Normal","Pre-DM") ~ 0,
    glucose_cat == "Diabetes" ~ 1
  ))

#check numbers
table(data_glucose$glucose_binary)#0- 498 , 1- 379 
class(data_glucose$glucose_binary)#numeric
data_glucose$glucose_binary <- as.factor(data_glucose$glucose_binary)
class(data_glucose$glucose_binary)#factor

#final variable count of 281,503 with 877 observations

#rearrange dataset so that all phenotype and PCA info is first and all SNPs are next
data_glucose <- data_glucose %>%
  relocate(DM:glucose_binary, .after = insulin_concentration)

#export final dataset 
write.csv(data_glucose, file = "GlycogenesSNPs_qc_PLINKandR.csv")

