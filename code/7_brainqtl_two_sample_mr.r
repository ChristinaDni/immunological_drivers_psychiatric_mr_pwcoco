##Immunomewide MR project
#Exposure data: Brain pQTLs
#This is an updated script in line with suggested revisions
#Date: 27/12/2024
#This script covers two sample mr


library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(tidyr)
gwasvcf::set_bcftools("/dir/bcftools-1.20/bcftools")

#First identify which genes are in the ukb inflammation panel

#Import the file that has the inflammation data

base<- read.table("/dir/data/ukb/clean_cd/metadata/inflammation_panel_synids.txt")

#Extract the HGNC gene from the file

base$hgnc<- sapply(strsplit(base$V2, "_"), "[", 1)

#Import the brainqtl data

brainqtl<- vroom::vroom("/dir/data/brain_pqtl_wingo/clean_cd/ROSMAP_Wingo_2021_pqtl_clean.txt")
brainqtl<- as.data.frame(brainqtl)

#Do the extraction

brainqtl_inflammation<- brainqtl[which(brainqtl$hgnc%in%base$hgnc),] #314 markers available
brainqtl_inflammation$pheno<- "BrainQTL_ROSMAP"
brainqtl_inflammation$N<- 400

#Split by gene

genes<- brainqtl_inflammation %>% split(x=brainqtl_inflammation, f=brainqtl_inflammation$hgnc)

#Create the loop for the analyses

lapply(list(genes), function(x){
  
  for(x in genes) {
  tryCatch({
 
#Taking a not of the gene, useful for the log file
log_hgnc_symbol<- x$hgnc[1]

  #Extract genome wide significant
  
  variants<- x[which(x$P<=5e-08),]
  
#Here keeping info for the log file as we move along the script
log_nvariants_significant<- nrow(variants)

if(nrow(variants)<1){
print("no variants")
log_ncis<- NA
log_ntrans<- NA
log_nvariants_after_clumping<- NA
log_nvariantsf<- NA
log_ninstruments<- NA

#Save the log file for information
setwd("/dir/proteomewide/brainqtl/inflammation/data")

log<- cbind(log_hgnc_symbol, log_nvariants_significant, log_ncis, log_ntrans, log_nvariants_after_clumping, log_nvariantsf, log_ninstruments)
if(file.exists("log_brainqtl.txt") == TRUE) {
write.table(log,"log_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}else{  
variants_formated<- format_data(
  dat=variants,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "pheno",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "hgnc",
  chr_col = "Chr",
  pos_col = "BP",
  gene_col = "ensemblid"
)

#Create two new cols so that the clumping format can work
variants_formated$rsid<- variants_formated$SNP
variants_formated$pval<- variants_formated$pval.exposure

try({variants_clumped<- ieugwasr::ld_clump(
  dat = variants_formated,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/dir/REF_PANEL/EUR",
  plink_bin ="/dir/plink_linux_x86_64_20240818/plink")})


if(!exists("variants_clumped")){
print("no variants retained during clumping")

log_nvariants_after_clumping<- 0
log_ncis<- NA
log_ntrans<- NA
log_nvariantsf<- NA
log_ninstruments<- NA

setwd("/dir/proteomewide/brainqtl/inflammation/data")

log<- cbind(log_hgnc_symbol, log_nvariants_significant, log_ncis, log_ntrans, log_nvariants_after_clumping, log_nvariantsf, log_ninstruments)
if(file.exists("log_brainqtl.txt") == TRUE) {
write.table(log,"log_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}else{  
log_nvariants_after_clumping<- nrow(variants_clumped)

#Anotate and estimate F statistic

exposure_instruments<-  mrpipeline::calc_f_stat(variants_clumped) %>%
  mrpipeline::cis_trans(values_col="gene.exposure", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "ensembl_gene_id", build = "grch37") #the dataset is in 37

#Here remove CHR denoted observations generated during the cistrans annotation
exposure_instruments<-  exposure_instruments %>% 
  filter(!grepl('CHR', chromosome_name))
  
#Info for log file

log_nvariants_cis<- exposure_instruments[which(exposure_instruments$cistrans=="C"),]
log_ncis<- nrow(log_nvariants_cis)
if(isTRUE(log_ncis)&&nrow(log_ncis<1)){
log_ncis<- 0}

log_nvariants_trans<- exposure_instruments[which(exposure_instruments$cistrans=="T"),]
log_ntrans<- nrow(log_nvariants_trans)
if(isTRUE(log_ntrans)&&nrow(log_ntrans<1)){
log_ntrans<- 0}

log_nvariants_f<- exposure_instruments[which(exposure_instruments$f.stat.exposure<10),]
log_nvariantsf<- nrow(log_nvariants_f)
if(isTRUE(log_nvariantsf)&&nrow(log_nvariantsf<1)){
log_nvariantsf<-0}

log_nvariants_i<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]
log_ninstruments<- nrow(log_nvariants_i)
if(isTRUE(log_ninstruments)&&nrow(log_ninstruments<1)){
log_ninstruments<-0}

setwd("/dir/proteomewide/brainqtl/inflammation/data")

log<- cbind(log_hgnc_symbol, log_nvariants_significant, log_ncis, log_ntrans, log_nvariants_after_clumping, log_nvariantsf, log_ninstruments)
if(file.exists("log_brainqtl.txt") == TRUE) {
write.table(log,"log_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}


#Entering in the analyses only instruments that pass the f statistic 10 threshold
exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

#Final file formating to save
exposure_instruments_save<- exposure_instruments%>%dplyr::select(-exposure, -id.exposure, -maf.exposure, -pve.exposure,
																	-chromosome_name, -start_position, -end_position,
																	-rsid, -pval, -id)
#Try to keep same naming convention across files
colnames(exposure_instruments_save)[which(names(exposure_instruments_save) == "gene.exposure")] <- "ensemblid"

#Save the file for the supplement
setwd("/dir/proteomewide/brainqtl/inflammation/data/instruments")
if(file.exists("instruments_rosmap.txt") == TRUE) {
  write.table(exposure_instruments_save,"instruments_rosmap.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(exposure_instruments_save, "instruments_rosmap.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }})


#Here MR analyses after instrument extraction
#Indicate the location of the outcome files
psychiatric_files<- list.files(path="/dir/proteomewide/psychiatric_data/clean", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

#Import the file that has the instruments
setwd("/dir/proteomewide/brainqtl/inflammation/data/instruments")
instruments<- read.table("instruments_rosmap.txt", header=T)

#First part of script: only cis variants

cis<- instruments[which(instruments$cistrans=="C"),]

#Split by gene

primary<- cis %>% split(x=cis, f=cis$hgnc_symbol)
setwd("/dir/proteomewide/brainqtl/inflammation/results/mr")

#Create the loop for the analyses

lapply(list(primary), function(x){
  
  for(x in primary) {
  tryCatch({
  
  #Create indicators with information on the SNPs
  x$id.exposure<- paste0(x$cistrans, ":", x$hgnc_symbol, ":", x$ensemblid, ":", "analysis", ":", "cis_only")
  x$exposure<- paste0(x$cistrans, ":", x$hgnc_symbol, ":", x$ensemblid, ":", "analysis", ":", "cis_only")
  iterations= 1:12

for (var in iterations){
try ({
outcome_primary <- TwoSampleMR::read_outcome_data(snps = x$SNP, filename = psychiatric_files[[var]], 
                        sep = "\t", snp_col = "SNP", beta_col = "logOR", 
                        se_col = "SE", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", chr_col= "CHR", pos_col="BP", id_col="pheno", 
						ncase_col="NCAS", ncontrol_col="NCON", samplesize_col="NTOT")

dat_primary <- mrpipeline::harmonise(x, outcome_primary, action = 2)
dat_primary<- dat_primary[which(dat_primary$mr_keep=="TRUE"),]

res_primary<-mrpipeline::do_mr(dat_primary, f_cutoff = 10)

#Here printing only IVW if multiple cis instruments exist
if(any(res_primary$method=="Inverse variance weighted")){

ivw_res_primary<- res_primary[which(res_primary$method=="Inverse variance weighted"),]

#Formating the snp info because it will be useful to merge the pwcoco results
ivw_res_primary<- ivw_res_primary %>% 
    mutate(snp = strsplit(as.character(snp), ",")) %>% 
    unnest(snp)
ivw_res_primary$snp<- gsub("\\s+", "", ivw_res_primary$snp)
	
ivw_res_primary<- as.data.frame(ivw_res_primary)
#Remove the pleiotropy columns from the output because they will be estimated separately
ivw_res_primary<- ivw_res_primary%>%dplyr::select(-egger_intercept,-se.egger,-pval.egger)

#Adding the heterogeneity

heterogeneity<- mr_heterogeneity(dat_primary, method_list = c("mr_ivw"))

heterogeneity<- heterogeneity%>%dplyr::select(id.exposure, Q, Q_df, Q_pval)

#Adding the pleiotropy

pleiotropy<- mr_pleiotropy_test(dat_primary)
pleiotropy<- pleiotropy%>%dplyr::select(id.exposure, egger_intercept, se, pval)
names(pleiotropy)[names(pleiotropy) == 'se'] <- 'se_egger'
names(pleiotropy)[names(pleiotropy) == 'pval'] <- 'pval_egger'

#Merge the tests

tests_primary<- merge(heterogeneity,pleiotropy, by="id.exposure")

#Merge the tests to the final dataframe
final_ivw_res_primary<- merge(ivw_res_primary,tests_primary, by="id.exposure", all.x=T)
print(final_ivw_res_primary)

setwd("/dir/proteomewide/brainqtl/inflammation/results/mr")
if(file.exists("ivw_primary_rosmap.txt") == TRUE) {
  write.table(final_ivw_res_primary,"ivw_primary_rosmap.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(final_ivw_res_primary, "ivw_primary_rosmap.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}


}

else{
wald_res_primary<- res_primary
print(wald_res_primary)
if(file.exists("wald_primary_rosmap.txt") == TRUE) {
  write.table(wald_res_primary,"wald_primary_rosmap.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(wald_res_primary, "wald_primary_rosmap.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
}

})
}
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }})
  
#Second part of script cis & trans analyses
#For formating reasons genes that have cis variants only will be excluded~ since they are covered in analyses above

secondary<- instruments %>% split(x=instruments, f=instruments$hgnc_symbol)
lapply(list(secondary), function(i){
  
  for(i in secondary) {
  tryCatch({

if(all(i$cistrans=="C")){
print("all variants cis")

}else{

#Create an indicator of cis & trans instruments included

i$cistrans_indicator<- paste0(i$cistrans, collapse = ":")

print(i)

  #Create indicators with information on the SNPs
  i$id.exposure<- paste0(i$cistrans_indicator,":", i$hgnc_symbol, ":", i$ensemblid, ":", "analysis", ":", "cis_trans")
  i$exposure<- paste0(i$cistrans_indicator, ":", i$hgnc_symbol, ":", i$ensemblid, ":", "analysis", ":", "cis_trans")
  iterations= 1:12

for (var in iterations){
try ({
outcome_secondary <- TwoSampleMR::read_outcome_data(snps = i$SNP, filename = psychiatric_files[[var]], 
                        sep = "\t", snp_col = "SNP", beta_col = "logOR", 
                        se_col = "SE", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", chr_col= "CHR", pos_col="BP", id_col="pheno", 
						ncase_col="NCAS", ncontrol_col="NCON", samplesize_col="NTOT")

dat_secondary <- mrpipeline::harmonise(i, outcome_secondary, action = 2)
dat_secondary<- dat_secondary[which(dat_secondary$mr_keep=="TRUE"),]

res_secondary<-mrpipeline::do_mr(dat_secondary, f_cutoff = 10)



#Here printing only IVW if multiple cis instruments exist
if(any(res_secondary$method=="Inverse variance weighted")){

ivw_res_secondary<- res_secondary[which(res_secondary$method=="Inverse variance weighted"),]

#Formating the snp info because it will be useful to merge the pwcoco results
ivw_res_secondary<- ivw_res_secondary %>% 
    mutate(snp = strsplit(as.character(snp), ",")) %>% 
    unnest(snp)
ivw_res_secondary$snp<- gsub("\\s+", "", ivw_res_secondary$snp)
	
ivw_res_secondary<- as.data.frame(ivw_res_secondary)
#Remove the pleiotropy columns from the output because they will be estimated separately
ivw_res_secondary<- ivw_res_secondary%>%dplyr::select(-egger_intercept,-se.egger,-pval.egger)

#Adding the heterogeneity

heterogeneity<- mr_heterogeneity(dat_secondary, method_list = c("mr_ivw"))

heterogeneity<- heterogeneity%>%dplyr::select(id.exposure, Q, Q_df, Q_pval)

#Adding the pleiotropy

pleiotropy<- mr_pleiotropy_test(dat_secondary)
pleiotropy<- pleiotropy%>%dplyr::select(id.exposure, egger_intercept, se, pval)
names(pleiotropy)[names(pleiotropy) == 'se'] <- 'se_egger'
names(pleiotropy)[names(pleiotropy) == 'pval'] <- 'pval_egger'

#Merge the tests

tests_secondary<- merge(heterogeneity,pleiotropy, by="id.exposure")

#Merge the tests to the final dataframe
final_ivw_res_secondary<- merge(ivw_res_secondary,tests_secondary, by="id.exposure", all.x=T)
print(final_ivw_res_secondary)

if(file.exists("ivw_secondary_rosmap.txt") == TRUE) {
  write.table(final_ivw_res_secondary,"ivw_secondary_rosmap.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(final_ivw_res_secondary, "ivw_secondary_rosmap.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}
else{
wald_res_secondary<- res_secondary
print(wald_res_secondary)
if(file.exists("wald_secondary_rosmap.txt") == TRUE) {
  write.table(wald_res_secondary,"wald_secondary_rosmap.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(wald_res_secondary, "wald_secondary_rosmap.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
}

})
}
}
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }})

