#Bi-directional Two-sample MR analyses of psychiatric phenotypes on all UKB immunological protein GWAS 
#This script:
#Imports all psychiatric data GWAS
#Extracts genome-wide significant and independent SNPs
#Depending on number of available independent SNPs applies two thresholds for extraction: if SNPs>5 p<=5e-08 or if SNPs <5 p<=5e-07
#Extracts from outcome, harmonises & performs two-sample MR
#Saves IVW output for supplementary material

#Load the necessary packages

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(tidyr)

gwasvcf::set_bcftools("/hpcpath/bcftools")

#List all exposure (psychiatric) and outcome (protein) data

exposure_files <- list.files(path="/path/psychiatric_data", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE) 

outcome_files<- list.files(path="/path/UKBGWAS_local_directory/", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

#Create functions for changing the threshold depending on availabilty of independent SNPs ie less or more than 5

typical_threshold <- function(k){
df<- k[which(k$P<=5e-08),]
exposure_i<- TwoSampleMR::format_data(df, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", 
				id_col="pheno", samplesize_col="NTOT", , ncase_col="NCAS", ncontrol_col="NCON", 
				phenotype_col="pheno", chr_col="CHR", pos_col="BP")
exposure_i$rsid<- exposure_i$SNP
exposure_i$pval<- exposure_i$pval.exposure

exposure_clumped<- ieugwasr::ld_clump(
  dat = exposure_i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/path/REF_PANEL/EUR",
  plink_bin ="/hpcpath/plink-1.90/plink")
  return(list(exposure_clumped))}
				
relaxed_threshold <- function(l){
df<- l[which(l$P<=5e-07),]
exposure_i<- TwoSampleMR::format_data(df, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", 
				id_col="pheno", samplesize_col="NTOT", , ncase_col="NCAS", ncontrol_col="NCON", 
				phenotype_col="pheno", chr_col="CHR", pos_col="BP")
exposure_i$rsid<- exposure_i$SNP
exposure_i$pval<- exposure_i$pval.exposure

exposure_clumped<- ieugwasr::ld_clump(
  dat = exposure_i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/path/REF_PANEL/EUR",
  plink_bin ="/hpcpath/plink-1.90/plink")
  return(list(exposure_clumped))}

#Load two-sample mr package

library(TwoSampleMR)

#Loop to analyse all psychiatric phenotypes on all proteins

lapply(list(exposure_files), function(x){
  for(x in exposure_files) {
  
tryCatch({
file<- vroom::vroom(x)
file_df<- as.data.frame(file)
file_i<- file_df[which(file_df$P<=5e-08),]
i<- TwoSampleMR::format_data(file_i, type = "exposure", snps = NULL, header = TRUE, 
                snp_col = "SNP", beta_col = "logOR", se_col = "SE", effect_allele_col = "A1", 
                other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", 
				id_col="pheno", samplesize_col="NTOT", , ncase_col="NCAS", ncontrol_col="NCON", 
				phenotype_col="pheno", chr_col="CHR", pos_col="BP")
i$rsid<- i$SNP
i$pval<- i$pval.exposure
clumped_test<- ieugwasr::ld_clump(
  dat = i,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/path/REF_PANEL/EUR",
  plink_bin ="/hpcpath/plink-1.90/plink")
print(nrow(clumped_test))

#Applying the functions to relax the threshold depending on availability of SNPs for analyses

final_exposure<- ifelse((nrow(clumped_test))>5, typical_threshold(file_df), relaxed_threshold(file_df))
head(final_exposure)
final_exposure<- as.data.frame(final_exposure)

setwd("/path/UKB/REVERSE/instruments")
if(file.exists("psychiatric_instruments.txt") == TRUE) {
  write.table(final_exposure,"psychiatric_instruments.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(final_exposure, "psychiatric_instruments.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#MR analyses across all available proteins

iterations= 1:741

for (var in iterations){
outcome <- TwoSampleMR::read_outcome_data(snps = final_exposure$SNP, filename = outcome_files[[var]], 
                        sep = "\t", snp_col = "rsid", beta_col = "BETA", 
                        se_col = "SE", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", eaf_col="A1FREQ", chr_col= "CHROM", pos_col="POS38", id_col="olink_assay_ukb", 
						samplesize_col="N", phenotype_col="HGNC.symbol", gene_col = "ensembl_id")
dat <- TwoSampleMR::harmonise_data(final_exposure, outcome)
dat<- dat[which(dat$mr_keep=="TRUE"),]

results<- mrpipeline::do_mr(dat, f_cutoff = 0)
setwd("/path/UKB/REVERSE/results")
if(file.exists("psychiatric_results.txt") == TRUE) {
  write.table(results,"psychiatric_results.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(results, "psychiatric_results.txt", sep="\t", quote=F, row.names=F, append=FALSE)}
  
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})


