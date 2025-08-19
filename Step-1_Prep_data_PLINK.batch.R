############################################################
################### SNPS Pipeline ##########################
#### Step -1: Prep SNPs datasets for PLINK #################
############# created by Mary 8-8-2025 #####################
############################################################

#this was run on the HPC
#to run this on your local computer I suggest subsetting out a section
#Goal is to take the raw dataset and create two PLINK ready dataset, the PED and the MAP files

#install packages
library(data.table)
library(tidyverse)
library(readxl)

#load data
SNPs_location = fread("~/PLINK/SNP_Information_Glycogene_SNPS_matched_rows_with_header.pvar")
data_SNPs = fread("~/PLINK/GlycogenesSNPsOUT.raw")

getwd()
SNPs_location = fread("~/Research/Projects/SNPs and synbiotics/Analysis/R code/SNP_Information_Glycogene_SNPS_matched_rows_with_header.pvar")

#get data into 2 columns per allele- then remove header
#subset data_SNPs into just genotypes
data_geno <- data_SNPs[,7:ncol(data_SNPs)]

# Create an empty list to store allele-separated data
allele_list <- list()

for (i in seq_along(colnames(data_geno))) {
  snp <- colnames(data_geno)[i]
  genotypes <- data_geno[[snp]]
  
  ref <- SNPs_location$REF[i]
  alt <- SNPs_location$ALT[i]
  
  # Create vectors for each allele
  allele1 <- rep(NA_character_, length(genotypes))
  allele2 <- rep(NA_character_, length(genotypes))
  
  allele1[genotypes == 0] <- alt
  allele2[genotypes == 0] <- alt
  
  allele1[genotypes == 1] <- ref
  allele2[genotypes == 1] <- alt
  
  allele1[genotypes == 2] <- ref
  allele2[genotypes == 2] <- ref
  
  # Add to list with appropriate names
  allele_list[[paste0(snp, "_1")]] <- allele1
  allele_list[[paste0(snp, "_2")]] <- allele2
}

# Combine into a new data frame
allele_data <- as.data.frame(allele_list)

# add family info back to finalize PED file
data_family <-data_SNPs[,1:6]
SNPs_PED <- as.data.frame(cbind(data_family, allele_data))

#sanity check 
head(SNPs_PED[,1:100])
colnames(SNPs_PED) <- NULL

#save as a tab deliminated file
write.table(SNPs_PED,
            file = "~/PLINK/SNPs_PED.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
)

#next step is just to check and make sure the order is still the same

#SNPs_location data prep
SNPs_location2 <- SNPs_location %>%
  rename("Chromosome" = "#CHROM", "Position"= "POS", "rs" = "ID")%>%
  select("Chromosome":"rs") %>%
  relocate("rs", .after = "Chromosome")

SNPs_location2 <- as.data.frame(SNPs_location2)
colnames(SNPs_location2) <- NULL


#save as a tab deliminated file
write.table(
  SNPs_location2,
  file = "~/Research/Projects/SNPs and synbiotics/Analysis/R code/SNPs.map",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
