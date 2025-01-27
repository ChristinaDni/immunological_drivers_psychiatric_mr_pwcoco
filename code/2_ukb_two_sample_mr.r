##Immunomewide MR project
#Exposure data: UKB
#This is an updated script in line with suggested revisions
#Date: 27/12/2024
#This script covers two sample mr


library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(tidyr)
gwasvcf::set_bcftools("/dir/bcftools-1.20/bcftools")

#Indicate the location of the outcome files
psychiatric_files<- list.files(path="/dir/proteomewide/psychiatric_data/clean", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

#Import the file that has the instruments
setwd("/dir/proteomewide/ukb/inflammation/data/instruments")
instruments<- read.table("instruments_ukb.txt", header=T)

#First part of script: only cis variants

cis<- instruments[which(instruments$cistrans=="C"),]

#Split by gene

primary<- cis %>% split(x=cis, f=cis$exposure)

#Create the loop for the analyses

lapply(list(primary), function(x){
  
  for(x in primary) {
  tryCatch({
  
  #Create indicators with information on the SNPs
  x$id.exposure<- paste0(x$exposure, ":", x$hgnc_symbol, ":", x$ensemblid, ":", "analysis", ":", "cis_only")
  x$exposure<- paste0(x$exposure, ":", x$hgnc_symbol, ":", x$ensemblid, ":", "analysis", ":", "cis_only")
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

#Here create a dataframe with the SNP and cistrans information

cistrans_info_primary<- dat_primary%>%dplyr::select(SNP,cistrans)

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

#Here merge the cistrans dataframe with the results so that this info is retained

ivw_res_primary<- merge(ivw_res_primary,cistrans_info_primary, by.x="snp",by.y="SNP")

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

setwd("/dir/proteomewide/ukb/inflammation/results/mr")
if(file.exists("ivw_primary_ukb.txt") == TRUE) {
  write.table(final_ivw_res_primary,"ivw_primary_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(final_ivw_res_primary, "ivw_primary_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}


}

else{
setwd("/dir/proteomewide/ukb/inflammation/results/mr")
wald_res_primary<- res_primary
wald_res_primary<- merge(wald_res_primary, cistrans_info_primary, by.x="snp", by.y="SNP")
print(wald_res_primary)
if(file.exists("wald_primary_ukb.txt") == TRUE) {
  write.table(wald_res_primary,"wald_primary_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(wald_res_primary, "wald_primary_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
}

})
}
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }})
  
#Second part of script cis & trans analyses
#For formating reasons genes that have cis variants only will be excluded~ since they are covered in analyses above

secondary<- instruments %>% split(x=instruments, f=instruments$exposure)
lapply(list(secondary), function(i){
  
  for(i in secondary) {
  tryCatch({

if(all(i$cistrans=="C")){
print("all variants cis")

}else{


  #Create indicators with information on the SNPs
  i$id.exposure<- paste0(i$exposure, ":", i$hgnc_symbol, ":", i$ensemblid, ":", "analysis", ":", "cis_trans")
  i$exposure<- paste0(i$exposure, ":", i$hgnc_symbol, ":", i$ensemblid, ":", "analysis", ":", "cis_trans")
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

#Here create a dataframe with the SNP and cistrans information

cistrans_info_secondary<- dat_secondary%>%dplyr::select(SNP,cistrans)

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

ivw_res_secondary<- merge(ivw_res_secondary,cistrans_info_secondary, by.x="snp",by.y="SNP")


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

if(file.exists("ivw_secondary_ukb.txt") == TRUE) {
  write.table(final_ivw_res_secondary,"ivw_secondary_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(final_ivw_res_secondary, "ivw_secondary_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}
else{
wald_res_secondary<- res_secondary
wald_res_secondary<- merge(wald_res_secondary, cistrans_info_secondary, by.x="snp", by.y="SNP")
print(wald_res_secondary)
if(file.exists("wald_secondary_ukb.txt") == TRUE) {
  write.table(wald_res_secondary,"wald_secondary_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(wald_res_secondary, "wald_secondary_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
}

})
}
}
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }})

