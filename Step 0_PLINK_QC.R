############################################################
################### SNPS Pipeline ##########################
############ Step -0.5: SNPs QC in PLINK ###################
############# created by Mary 8-8-2025 #####################
############################################################

#This code was run in PLINK on the command line and stored in R for reproducibility
#to run code, make sure plink and the PED and MAP files you are working with are in the same folder

#open command line and set working directory to that folder
cd "file path"

#check to make sure you are in the right folder by looking at it
dir

#need to make sure files are listed correctly. Both files should have the same name but one should be .map and the other .ped
#rename files and remove .txt that will follow the dataset if you followed the code in Step -1
ren SNPs.map.txt SNPs.map
ren SNPs.ped.txt SNPs.ped

#now check missingness of SNPs
plink --file SNPs --mind 0.1 --make-bed 

#now you can run one code to include all of the QC at once
plink --file SNPs --maf 0.05 --hwe 0.001 --indep-pairwise 50 5 0.5 --make-bed #all of these numbers are the default conditions

#it will automatically create output files for you. The one you are looking for is pruned in

#results
#156 variants removed with HWE
#243340 variants removed for MAF
#none removed for missingness of SNP
#43399 variants removed with pruning

#Total left after QC was 1025 individuals and 7,541 SNPs