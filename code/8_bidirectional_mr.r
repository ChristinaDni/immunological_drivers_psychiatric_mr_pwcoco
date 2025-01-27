#Reverse direction analyses
#Updated- across all exposure phenotypes

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(tidyr)

gwasvcf::set_bcftools("/dir/bcftools-1.20/bcftools")

exposure_files <- list.files(path="/dir/proteomewide/psychiatric_data/clean", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE) 

outcome_files<- list.files(path="/dir/data/ukb/clean_cd/gwas/inflammation", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

#Create the functions to change the thresholds if necessary

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
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink")
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
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink")
  return(list(exposure_clumped))}


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
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink")
print(nrow(clumped_test))

final_exposure<- ifelse((nrow(clumped_test))>5, typical_threshold(file_df), relaxed_threshold(file_df))
print(head(final_exposure))
final_exposure<- as.data.frame(final_exposure)
print(head(final_exposure))

setwd("/dir/proteomewide/reverse/instruments")
if(file.exists("psychiatric_instruments.txt") == TRUE) {
  write.table(final_exposure,"psychiatric_instruments.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(final_exposure, "psychiatric_instruments.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

iterations= 1:743

for (var in iterations){
outcome <- TwoSampleMR::read_outcome_data(snps = final_exposure$SNP, filename = outcome_files[[var]], 
                        sep = "\t", phenotype_col = "file", snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "HGNC_symbol",
  chr_col = "CHR",
  pos_col = "POS37",
  gene_col = "Ensembl_ID",
  eaf_col= "EAF"
)

dat <- TwoSampleMR::harmonise_data(final_exposure, outcome)
dat<- dat[which(dat$mr_keep=="TRUE"),]

results<- mrpipeline::do_mr(dat)
results<- results[which(results$method=="Inverse variance weighted"),]
setwd("/dir/proteomewide/reverse/results")
if(file.exists("psychiatric_results.txt") == TRUE) {
  write.table(results,"psychiatric_results.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(results, "psychiatric_results.txt", sep="\t", quote=F, row.names=F, append=FALSE)} 
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})


