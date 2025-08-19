#SNPs logistic regression prep and mode for HbA1c
# finalized 6-26-25

#install packages
library(data.table)
library(readxl)
library(tidyverse)

#load data
data_SNPs = fread("GlycogenesSNPsOUT.raw")
data_phenotype = read_excel("Data Phenotype for Transfer 1.15.25-Sri (1).xlsx")
data_pca <- fread("PCA_wBasicQC (1).eigenvec", header = FALSE)
setnames(data_pca, c("FID", "IID", paste0("PCA", 1:(ncol(data_pca) - 2))))

#prep to merge datasets by IID
data_SNPs2 <- data_SNPs %>%
  select(-c("FID", "MAT", "PAT", "PHENOTYPE", "SEX"))

data_phenotype2 <- data_phenotype %>%
  rename(IID = biobank_id)%>%
  rename(gender_2 = "gender_2 (Female is 0)")%>%
  select(-c("record_id", "AT_number", "Height (inches)"))

data_pca2 <- data_pca %>%
  select(-"FID")

#merge datasets
data_merged <- merge(data_phenotype2, data_pca2, by = "IID")
data_merged <- merge(data_merged, data_SNPs2, by = "IID")

####merged dataset is 1026 observations and 294,472 variables

### SNPs Quality control 
#remove individuals with more than 2% missing SNPs
for(i in 1:1026){
  data_merged$percentmissing[i] <-sum(is.na(data_merged[i,])/294472)*100
  print(i)
}

summary(data_merged$percentmissing)#none were above 2% so the dataset is good

#remove percent missing to make it more stream line
data_merged <- data_merged %>%
  select(-"percentmissing")

#make sure important numeric columns are read as numeric
data_merged$glucose <- as.numeric(data_merged$glucose)
data_merged$hgb_a1c_glyco <- as.numeric(data_merged$hgb_a1c_glyco)

#make sure there are no missing values in variables that will be in the models
#glucose and hba1c will be done later for their respective datasets
#age and gender will be done now
#confirmed there is no missing for age using JMP, but 6 missing for gender

summary(data_merged$age_2) #confirmed no NA
table(data_merged$gender_2) #confirmed 6

data_merged$gender_2 <- as.numeric(data_merged$gender_2)

data_final <- data_merged %>%
  drop_na(gender_2)

data_final$gender_2 <- as.factor(data_final$gender_2)

#sanity check
summary(data_final$age_2) #confirmed no NA
table(data_final$gender_2) #confirmed no NA

####final dataset is 1020 observations and 294,472 variables

#clean up SNPs
#remove all SNP columns that have all of the same values
# Identify columns (SNPs) with only one unique value
columns_to_remove <- sapply(data_final[, 1:294472], function(x) length(unique(na.omit(x))) < 2)

# Remove those SNPs from the dataset
data_removed <- data_final[, !columns_to_remove]

# Check how many SNPs were removed
sum(columns_to_remove)
#n= 12,973 SNPS removed due to this

#remove SNPs where results are all NA
data_SNPs_allNA <- data_removed %>%
  select(where(~ !all(is.na(.)))) #no columns with all NAs

###final cleaned dataset (data removed)
#1020 observations
#281499 variables


### create DM variable ###
##categorical variable
#filter out rows where hgb_a1c_glyco is NA for complete cases
summary(data_removed$hgb_a1c_glyco) #confirmed X NA

data_hbA1c <- data_removed %>%
  filter(!is.na(hgb_a1c_glyco))

summary(data_hbA1c$hgb_a1c_glyco) #confirmed no NA

#now there are 982 complete cases for hbA1c and 281,499 variables

#create DM categories
data_hbA1c <- data_hbA1c %>%
  mutate(DM = case_when(
    hgb_a1c_glyco <5.7 ~ "Normal",
    hgb_a1c_glyco >=5.7 & hgb_a1c_glyco <=6.4 ~ "Pre-DM",
    hgb_a1c_glyco >=6.5 ~ "Diabetes"
  ))

#check numbers/ sanity check
table(data_hbA1c$DM) #Normal- 216 , PreDM- 195 , DM- 571, no NA
class(data_hbA1c$DM) #character

#create binary variable
data_hbA1c <- data_hbA1c %>%
  mutate(DM_binary = case_when(
    DM %in% c("Normal","Pre-DM") ~ 0,
    DM == "Diabetes" ~ 1
  ))

#check numbers/ sanity check 
table(data_hbA1c$DM_binary) #0- 411, 1- 571, no NA
class(data_hbA1c$DM_binary) #numeric
data_hbA1c$DM_binary <- as.factor(data_hbA1c$DM_binary)
class(data_hbA1c$DM_binary) #factor

#should have final variable count of 281,501 with 982 observations
#rearrange dataset so that all phenotype and PCA info is first and all SNPs are next
data_hbA1c <- data_hbA1c %>%
  relocate(DM:DM_binary, .after = insulin_concentration)

#export final dataset for complete cases using HbA1c
write.csv(data_hbA1c, file = "GlycogenesSNPs_qc_HbA1c.csv")


#### Logistic Regression Model ####
# Create blank OUT file
OUT_log_HbA1c <- data.frame(SNP = character(),
                            Coef = numeric(),
                            StdErr = numeric(),
                            OR = numeric(),
                            Pval = numeric(),
                            stringsAsFactors = FALSE)

# Loop through all SNPs
for (i in 39:281501) {  # Assuming first column is ID
  SNP <- data_SNPs_HbA1c[[i]]
  model <- glm(DM_binary ~ SNP + gender_2 + age_2 + PCA1 + PCA2, data = data_SNPs_HbA1c, family = binomial(link = "logit"))
  coef_vals <- summary(model)$coef
  OUT_log_HbA1c[i-38, ] <- c(colnames(data_SNPs_HbA1c)[i],  # SNP name
                             coef_vals[2, 1],                      # Coefficient
                             coef_vals[2, 2],                      # Standard Error
                             exp(coef_vals[2, 1]),                 # Odds Ratio
                             coef_vals[2, 4])                      # P-value
}

# Save OUT file
write.table(OUT_log_HbA1c, "SNP_Logistic Regression_HbA1c_OUT.txt", quote = FALSE, row.names = FALSE, sep = "\t")
