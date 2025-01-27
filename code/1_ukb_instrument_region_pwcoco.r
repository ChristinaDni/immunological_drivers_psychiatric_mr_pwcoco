#Immunomewide MR project
#Exposure data: UKB
#This is an updated script in line with suggested revisions
#Date: 27/12/2024
#In this script, instrument, exposure& outcome region extraction, colocalisation happen all at once

###

#Script starts here

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(data.table)
gwasvcf::set_bcftools("/dir/bcftools-1.20/bcftools")

setwd("/dir/data/ukb/clean_cd/gwas/inflammation")
ukb<- list.files(path = "/dir/data/ukb/clean_cd/gwas/inflammation", recursive = TRUE)

#Initiate loop to run over gwas data

lapply(list(ukb), function(k){
  
  for(k in ukb) {
  
  tryCatch({
  
setwd("/dir/data/ukb/clean_cd/gwas/inflammation")

data<- vroom::vroom(k)
data<- as.data.frame(data)

#Keeping gene info for log file

log_hgnc_symbol<- data$HGNC_symbol[1]

variants<- data[which(data$P<=5e-08),]

#Keeping number of genome wide significant variants for log file

log_nvariants_significant<- nrow(variants)

if(nrow(variants)<1){
log_ncis<- NA
log_ntrans<- NA
log_nvariants_after_clumping<- NA
log_nvariantsf<- NA
log_ninstruments<- NA

setwd("/dir/proteomewide/ukb/inflammation/data")

log<- cbind(log_hgnc_symbol, log_nvariants_significant, log_ncis, log_ntrans, log_nvariants_after_clumping, log_nvariantsf, log_ninstruments)
if(file.exists("log_ukb.txt") == TRUE) {
write.table(log,"log_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}else{  
variants_formated<- format_data(
  dat=variants,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "file",
  snp_col = "SNP",
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

#Keeping number of clumped variants for log file

if(!exists("variants_clumped")){
print("no variants retained during clumping")
log_nvariants_after_clumping<- 0
log_ncis<- NA
log_ntrans<- NA
log_nvariantsf<- NA
log_ninstruments<- NA

setwd("/dir/proteomewide/ukb/inflammation/data")

log<- cbind(log_hgnc_symbol, log_nvariants_significant, log_ncis, log_ntrans, log_nvariants_after_clumping, log_nvariantsf, log_ninstruments)
if(file.exists("log_ukb.txt") == TRUE) {
write.table(log,"log_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}else{  

log_nvariants_after_clumping<- nrow(variants_clumped)

#Estimate F and annotate cis/trans

variants_final<-  mrpipeline::calc_f_stat(variants_clumped) %>%
  mrpipeline::cis_trans(values_col="gene.exposure", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "ensembl_gene_id", build = "grch37")
 
#Here remove CHR denoted observations generated during the cistrans annotation
variants_final<-  variants_final %>% 
  filter(!grepl('CHR', chromosome_name))

#Keeping number of cis/trans and f passing/not passing variants for log file
  
log_nvariants_cis<- variants_final[which(variants_final$cistrans=="C"),]
log_ncis<- nrow(log_nvariants_cis)
if(isTRUE(log_ncis)&&nrow(log_ncis<1)){
log_ncis<- 0}

log_nvariants_trans<- variants_final[which(variants_final$cistrans=="T"),]
log_ntrans<- nrow(log_nvariants_trans)
if(isTRUE(log_ntrans)&&nrow(log_ntrans<1)){
log_ntrans<- 0}

log_nvariants_f<- variants_final[which(variants_final$f.stat.exposure<10),]
log_nvariantsf<- nrow(log_nvariants_f)
if(isTRUE(log_nvariantsf)&&nrow(log_nvariantsf<1)){
log_nvariantsf<-0}

log_nvariants_i<- variants_final[which(variants_final$f.stat.exposure>=10),]
log_ninstruments<- nrow(log_nvariants_i)
if(isTRUE(log_ninstruments)&&nrow(log_ninstruments<1)){
log_ninstruments<-0}

setwd("/dir/proteomewide/ukb/inflammation/data")

log<- cbind(log_hgnc_symbol, log_nvariants_significant, log_ncis, log_ntrans, log_nvariants_after_clumping, log_nvariantsf, log_ninstruments)
if(file.exists("log_ukb.txt") == TRUE) {
write.table(log,"log_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}


#Entering in the analyses only instruments that pass the f statistic 10 threshold
variants_final<- variants_final[which(variants_final$f.stat.exposure>=10),]

#Annotate those in the MHC region- they will not be removed but the information will be retained
variants_final$MHC<- ifelse(variants_final$chr.exposure=="6"&(variants_final$pos.exposure<34000000 & variants_final$pos.exposure>25000000), "caution_MHC", "no_MHC")


#Final file formating to save

variants_final_save<- variants_final%>%dplyr::select(-id,-id.exposure, -maf.exposure, -pve.exposure,
																	-chromosome_name, -start_position, -end_position, -rsid, -pval)

colnames(variants_final_save)[which(names(variants_final_save) == "gene.exposure")] <- "ensemblid"


#Save the file for the supplement
setwd("/dir/proteomewide/ukb/inflammation/data/instruments")
if(file.exists("instruments_ukb.txt") == TRUE) {
write.table(variants_final_save,"instruments_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(variants_final_save, "instruments_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#Region extraction
#The approach for extraction and subsequent colocalisation is prioritise cis variants
#If no cis available then the trans with the smallest p will be extracted

if(any(variants_final_save$cistrans=="C")){

instruments<- variants_final_save[which(variants_final_save$cistrans=="C"),]
instruments$exposure<- gsub(":", "_", instruments$exposure) #Here the replacement is happening because the extraction of the chromosome with PWCOCO will be difficult~ it is lenghtier than the usual.

map<- instruments %>% split(x=instruments, f=instruments$SNP)

  lapply(list(map), function(l){
  
  for(l in map) {

region<- l%>%dplyr::select("exposure","SNP","chr.exposure","pos.exposure", "ensemblid", "hgnc_symbol", "cistrans", "MHC")

#Here we create a list with eligible SNPs so that it can be comparable with the regions
region$TopSNP<- paste0((l$exposure),
                     ":",
                     (l$SNP), ":", (l$chr.exposure),
                     ":", (l$pos.exposure), ":", (l$ensemblid),
					 ":",
                     (l$hgnc_symbol), ":", (l$cistrans), ":", (l$MHC))


setwd("/dir/proteomewide/ukb/inflammation/data/regions")
if(file.exists("topregions_ukb_inflammation.txt") == TRUE) {
  write.table(region,"topregions_ukb_inflammation.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_ukb_inflammation.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#The following line creates an indicator which will be useful for saving the results
snp_list<- paste0(region$exposure, ":", region$SNP,":",region$chr.exposure, ":", region$pos.exposure, ":", region$ensemblid, ":", region$hgnc_symbol, ":", region$cistrans, ":", region$MHC)
snp_list<- as.data.frame(snp_list)

#Here using the full dataset that was used initially so that regions can be extracted

#Define the extraction window
window<- 500000

all_chr<- data[which(data$CHR==region$chr.exposure),]
all_pos<- all_chr[which(all_chr$POS37>=max(region$pos.exposure-window, 0) & all_chr$POS37<=region$pos.exposure+window),]

variants_region <- all_pos %>% 
    dplyr::select("SNP", "A1", "A2", "EAF", "BETA", "SE", "P", "N") 
	
variants_region<- as.data.frame(variants_region)
	
print(head(variants_region))
print(nrow(variants_region))

#Saving the output
setwd("/dir/proteomewide/ukb/inflammation/data/regions/exposure")

snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}
}})

}else{
instruments<- variants_final_save[which.min(variants_final_save$pval.exposure),]
instruments$exposure<- gsub(":", "_", instruments$exposure) #Here the replacement is happening because the extraction of the chromosome with PWCOCO will be difficult~ it is lenghtier than the usual.

region<- instruments%>%dplyr::select("exposure","SNP","chr.exposure","pos.exposure", "ensemblid", "hgnc_symbol", "cistrans", "MHC")

#Here we create a list with eligible SNPs so that it can be comparable with the regions
region$TopSNP<- paste0((instruments$exposure),
                     ":",
                     (instruments$SNP), ":", (instruments$chr.exposure),
                     ":", (instruments$pos.exposure), ":", (instruments$ensemblid),
					 ":",
                     (instruments$hgnc_symbol), ":", (instruments$cistrans), ":", (instruments$MHC))


setwd("/dir/proteomewide/ukb/inflammation/data/regions")
if(file.exists("topregions_ukb_inflammation.txt") == TRUE) {
  write.table(region,"topregions_ukb_inflammation.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_ukb_inflammation.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#The following line creates an indicator which will be useful for saving the results
snp_list<- paste0(region$exposure, ":", region$SNP,":",region$chr.exposure, ":", region$pos.exposure, ":", region$ensemblid, ":", region$hgnc_symbol, ":", region$cistrans, ":", region$MHC)
snp_list<- as.data.frame(snp_list)

#Here using the full dataset that was created initially so that regions can be extracted

#Define the extraction window
window<- 500000

all_chr<- data[which(data$CHR==region$chr.exposure),]
all_pos<- all_chr[which(all_chr$POS37>=max(region$pos.exposure-window, 0) & all_chr$POS37<=region$pos.exposure+window),]

variants_region <- all_pos %>% 
    dplyr::select("SNP", "A1", "A2", "EAF", "BETA", "SE", "P", "N") 
	
variants_region<- as.data.frame(variants_region)
	
print(head(variants_region))
print(nrow(variants_region))

#Saving the output
setwd("/dir/proteomewide/ukb/inflammation/data/regions/exposure")

snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")
write.table(variants_region, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)}
}
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})


#Next part of the script: extraction from psychiatric data
 
region_files<- list.files(path="/dir/proteomewide/ukb/inflammation/data/regions/exposure", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

psychiatric_files<- list.files(path="/dir/proteomewide/psychiatric_data/clean", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

lapply(list(region_files), function(m){
  for(m in region_files) {

tryCatch({

region<- read.table(m, header=T, sep="\t")
print(nrow(region))

#Create a variable that indicates the name of the file~ this is useful for naming the psychiatric files later

name<- basename(m)
print(name)

iterations= 1:12

for (var in iterations){

psychiatric_data<- vroom::vroom(psychiatric_files[[var]])
psychiatric_data<- as.data.frame(psychiatric_data)

#Select the variants from the molecular region in the psychiatric data
psychiatric_region<- psychiatric_data[which(psychiatric_data$SNP%in%region$SNP),]

#Create name for the files that will be exported
name_pheno<- psychiatric_region$pheno[1]
psychiatric_region_file_name<- paste0(name, ":", name_pheno, ".txt")
print(psychiatric_region_file_name)

#Keep the information for the pwcoco
psychiatric_region<- psychiatric_region%>%dplyr::select(SNP,A1,A2,FRQ,logOR,SE,P,NTOT,NCAS)
print(nrow(psychiatric_region))

#Save each file with the filename that was created
setwd("/dir/proteomewide/ukb/inflammation/data/regions/outcome")
write.table(psychiatric_region, psychiatric_region_file_name, sep="\t", quote=F, row.names=F, col.names=T)

}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}})

#Next part of script: pwcoco analyses

mol_region_files<- list.files(path="/dir/proteomewide/ukb/inflammation/data/regions/exposure", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

psych_region_files<- list.files(path="/dir/proteomewide/ukb/inflammation/data/regions/outcome", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

lapply(list(mol_region_files), function(p){
  for(p in mol_region_files) {

tryCatch({

#Creating the commands for running pwcoco

name<- basename(p)
print(name)

chr=sapply(strsplit(name, ":"), "[", 3) #Check here whether replacement of : with _ extracted correctly the chromosome
print(chr)

maf=0.01

#Here creating a pattern variable so that files for a given region

pattern=paste0("*",name,"*")
print(pattern)

#List all files for a give Region

mol_trait <- list.files(path = "/dir/proteomewide/ukb/inflammation/data/regions/exposure", pattern=pattern, recursive = TRUE)
print(mol_trait)
pheno<- list.files(path = "/dir/proteomewide/ukb/inflammation/data/regions/outcome", pattern=pattern, recursive = TRUE)
print(pheno)

#Assign to the identified files the full directory so that pwcoco can work 
mol_trait_dir<- paste0("/dir/proteomewide/ukb/inflammation/data/regions/exposure/", mol_trait)
print(mol_trait_dir)
pheno_dir<- paste0("/dir/proteomewide/ukb/inflammation/data/regions/outcome/", pheno)
print(pheno_dir)

setwd("/dir/pwcoco_update/pwcoco/build")

#Build the commands
bim_cmd<- paste("--bfile", paste0("/dir/REF_PANEL/refpanel_alspac/CHR", chr))
print(bim_cmd)
chr_cmd<- paste("--chr", chr)
print(chr_cmd)
maf_cmd<- paste("--maf", maf)
print(maf_cmd)

lapply(list(mol_trait_dir), function(r){
  
for(r in mol_trait_dir) {

sumstats1<- paste("--sum_stats1", r)
print(sumstats1)

lapply(list(pheno_dir), function(y){
  
  for(y in pheno_dir) {

sumstats2<- paste("--sum_stats2", y)
print(sumstats2)

out_cmd<- paste("--out", paste0("/dir/proteomewide/ukb/inflammation/results/colocalisation/ukb"))
print(out_cmd)

cmd <- paste("./pwcoco", bim_cmd,
				sumstats1,  
				sumstats2,
				out_cmd,
				chr_cmd, maf_cmd)
print(cmd)
system(cmd)


}})

}})

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}})

