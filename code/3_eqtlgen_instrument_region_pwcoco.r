#Immunomewide MR project
#Exposure data: EQTLGEN
#This is an updated script in line with suggested revisions
#Date: 27/12/2024
#In this script, instrument, exposure& outcome region extraction, colocalisation happen all at once

###

#Script starts here

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
gwasvcf::set_bcftools("/dir/bcftools-1.20/bcftools")

#Import the file that has the inflammation data

base<- read.table("/dir/data/ukb/clean_cd/metadata/inflammation_panel_synids.txt")

#Extract the HGNC gene from the file

base$hgnc<- sapply(strsplit(base$V2, "_"), "[", 1)

#Convert the hgnc to ensembl to match the eqtlgen data

mart.gene <- biomaRt::useMart(biomart = "ensembl",
                                  host="https://www.ensembl.org",
                                  dataset = "hsapiens_gene_ensembl")
								  
ensembl <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filters = "hgnc_symbol",
                            values = unique(base$hgnc),
                            mart = mart.gene)
							
base_final <- merge(base, ensembl, by.x="hgnc", by.y="hgnc_symbol")	

#Create the path to each file so that data can be extracted

base_final$path_name <- sub("^", "/dir/data/public/eqtl-a-", base_final$ensembl_gene_id)

files<- list.files(path=base_final$path_name, 
pattern="*.vcf.gz$", full.names=TRUE, recursive=FALSE) #559 files available~ from total 736 available in inflammation panel

###

#The instrument extraction, annotation and region extraction starts here:

lapply(list(files), function(x){
  for(x in files) {

tryCatch({

exposure_instruments<- read_exposure(x, pval = 5e-8, plink="/dir/plink_linux_x86_64_20240818/plink", 
bfile= "/dir/REF_PANEL/EUR", clump_r2=0.001) #Update after revision~ clumping threshold

#Keeping info for log file
if(nrow(exposure_instruments)<1){
log_id<- x
log_nvariants_significant<- 0
log_ncis<- NA
log_ntrans<- NA
log_nvariantsf<- NA
log_ninstruments<- NA

setwd("/dir/proteomewide/eqtlgen/inflammation/data")

log<- cbind(log_id, log_nvariants_significant, log_ncis, log_ntrans, log_nvariantsf, log_ninstruments)
if(file.exists("log_eqtlgen.txt") == TRUE) {
write.table(log,"log_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

}else{ 

#Format the df to keep necessary columns: 
#Uiversal format (for all outputs across analyses)
exposure_instruments<- exposure_instruments%>%dplyr::select(-ncase.exposure, -ncontrol.exposure, 
														-rsid, -pval, -id, -file.exposure)

exposure_instruments$ensemblid<- gsub('eqtl-a-','', exposure_instruments$exposure)

#Estimate F and annotate cis/trans
exposure_instruments<-  mrpipeline::calc_f_stat(exposure_instruments) %>%
  mrpipeline::cis_trans(values_col="ensemblid", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "ensembl_gene_id", build = "grch37") #eqtlgen is in hg37

#Here remove CHR denoted observations generated during the cistrans annotation
exposure_instruments<-  exposure_instruments %>% 
  filter(!grepl('CHR', chromosome_name))

#Keeping info for log file  
log_id<- exposure_instruments$hgnc_symbol[1]
log_nvariants_significant<- nrow(exposure_instruments)
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

setwd("/dir/proteomewide/eqtlgen/inflammation/data")

log<- cbind(log_id, log_nvariants_significant, log_ncis, log_ntrans, log_nvariantsf, log_ninstruments)
if(file.exists("log_eqtlgen.txt") == TRUE) {
write.table(log,"log_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
write.table(log, "log_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#Entering in the analyses only instruments that pass the f statistic 10 threshold
exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

#Annotate those in the MHC region- they will not be removed but the information will be retained
exposure_instruments$MHC<- ifelse(exposure_instruments$chr.exposure=="6"&(exposure_instruments$pos.exposure<34000000 & exposure_instruments$pos.exposure>25000000), "caution_MHC", "no_MHC")


#Final file formating to save
exposure_instruments_save<- exposure_instruments%>%dplyr::select(-exposure, -id.exposure, -maf.exposure, -pve.exposure,
																	-chromosome_name, -start_position, -end_position)

#Save the file for the supplement
setwd("/dir/proteomewide/eqtlgen/inflammation/data/instruments")
if(file.exists("instruments_eqtlgen.txt") == TRUE) {
  write.table(exposure_instruments_save,"instruments_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(exposure_instruments_save, "instruments_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#Region extraction
#The approach for extraction and subsequent colocalisation is prioritise cis variants
#If no cis available then the trans with the smallest p will be extracted

#Here addition for correct region extraction
#Ie in previous iteration positions <500000 gave negative lower limit and therefore were not extracted
#Now I have adapted the code from Jamies mr pipeline github to do correct extraction

#Denoting the window
window<- 500000

#Starting the extraction

if(any(exposure_instruments_save$cistrans=="C")){

instruments<- exposure_instruments_save[which(exposure_instruments_save$cistrans=="C"),]

map<- instruments %>% split(x=instruments, f=instruments$SNP)

lapply(list(map), function(i){
  
  for(i in map) {

region<- i[c("SNP","chr.exposure","pos.exposure", "ensemblid", "hgnc_symbol", "cistrans", "MHC")]

#Here we create a list with eligible SNPs so that it can be comparable with the regions
region$TopSNP<- paste0(
                     (i$SNP), ":", (i$chr.exposure),
                     ":", (i$pos.exposure), ":", (i$ensemblid),
					 ":",
                     (i$hgnc_symbol), ":", (i$cistrans), ":", (i$MHC))

setwd("/dir/proteomewide/eqtlgen/inflammation/data/regions")
if(file.exists("topregions_eqtlgen_inflammation.txt") == TRUE) {
  write.table(region,"topregions_eqtlgen_inflammation.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_eqtlgen_inflammation.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#The following line creates an indicator which will be useful for saving the results
snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$pos.exposure, ":", region$ensemblid, ":", region$hgnc_symbol, ":", region$cistrans, ":", region$MHC)
snp_list<- as.data.frame(snp_list)

#Creating chrpos indicator, useful for the extraction of regions through gwasvcf
chrpos <- paste0(as.character(region$chr.exposure),
                     ":",
                     max(as.numeric(region$pos.exposure) - window, 0),
                     "-",
                     as.numeric(region$pos.exposure) + window)

variants_region <- gwasvcf::query_gwas(x, chrompos = chrpos)

variants_region_tibble <- variants_region %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
    dplyr::select(dplyr::any_of(c("ID", "ALT", "REF", "AF", "ES", "SE", "LP", "SS"))) 
    
variants_region_final<- variants_region_tibble %>% dplyr::rename(
      SNP = ID,
      A1 = ALT,
      A2 = REF,
      freq = AF,
      b = ES,
      se = SE,
      p = LP,
      N = SS,
    ) #Note here that all files will follow the same naming convention

#Transforming the pvalue as is stored in vcf format	
variants_region_final$p <- 10^(-variants_region_final$p)

variants_region_final<- as.data.frame(variants_region_final)

print(nrow(variants_region_final))

#Saving the output
setwd("/dir/proteomewide/eqtlgen/inflammation/data/regions/exposure")

snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region_final, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}

}
})
}else{
instruments<- exposure_instruments_save[which.min(exposure_instruments_save$pval.exposure),]
region<- instruments[c("SNP","chr.exposure","pos.exposure", "ensemblid", "hgnc_symbol", "cistrans", "MHC")]

#Here we create a list with eligible SNPs so that it can be comparable with the regions
region$TopSNP<- paste0(
                     (instruments$SNP), ":", (instruments$chr.exposure),
                     ":", (instruments$pos.exposure), ":", (instruments$ensemblid),
					 ":",
                     (instruments$hgnc_symbol), ":", (instruments$cistrans), ":", (instruments$MHC))

setwd("/dir/proteomewide/eqtlgen/inflammation/data/regions")
if(file.exists("topregions_eqtlgen_inflammation.txt") == TRUE) {
  write.table(region,"topregions_eqtlgen_inflammation.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_eqtlgen_inflammation.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}


#The following line creates an indicator which will be useful for saving the results
snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$pos.exposure, ":", region$ensemblid, ":", region$hgnc_symbol, ":", region$cistrans, ":", region$MHC)
snp_list<- as.data.frame(snp_list)

#Creating chrpos indicator, useful for the extraction of regions through gwasvcf
chrpos <- paste0(as.character(region$chr.exposure),
                     ":",
                     max(as.numeric(region$pos.exposure) - window, 0),
                     "-",
                     as.numeric(region$pos.exposure) + window)

variants_region <- gwasvcf::query_gwas(x, chrompos = chrpos)

variants_region_tibble <- variants_region %>% gwasvcf::vcf_to_granges() %>% dplyr::as_tibble() %>%
    dplyr::select(dplyr::any_of(c("ID", "ALT", "REF", "AF", "ES", "SE", "LP", "SS"))) 
    
variants_region_final<- variants_region_tibble %>% dplyr::rename(
      SNP = ID,
      A1 = ALT,
      A2 = REF,
      freq = AF,
      b = ES,
      se = SE,
      p = LP,
      N = SS,
    ) #Note here that all files will follow the same naming convention

#Transforming the pvalue as is stored in vcf format	
variants_region_final$p <- 10^(-variants_region_final$p)

variants_region_final<- as.data.frame(variants_region_final)

print(nrow(variants_region_final))

#Saving the output
setwd("/dir/proteomewide/eqtlgen/inflammation/data/regions/exposure")

snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")

write.table(variants_region_final, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)}
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}})

#Next part of the script: extraction from psychiatric data
 
region_files<- list.files(path="/dir/proteomewide/eqtlgen/inflammation/data/regions/exposure", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

psychiatric_files<- list.files(path="/dir/proteomewide/psychiatric_data/clean", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

lapply(list(region_files), function(x){
  for(x in region_files) {

tryCatch({

region<- read.table(x, header=T, sep="\t")
print(nrow(region))

#Create a variable that indicates the name of the file~ this is useful for naming the psychiatric files later

name<- basename(x)
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
setwd("/dir/proteomewide/eqtlgen/inflammation/data/regions/outcome")
write.table(psychiatric_region, psychiatric_region_file_name, sep="\t", quote=F, row.names=F, col.names=T)

}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}})

#Next part of script: pwcoco analyses

mol_region_files<- list.files(path="/dir/proteomewide/eqtlgen/inflammation/data/regions/exposure", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

psych_region_files<- list.files(path="/dir/proteomewide/eqtlgen/inflammation/data/regions/outcome", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

lapply(list(mol_region_files), function(x){
  for(x in mol_region_files) {

tryCatch({

#Creating the commands for running pwcoco

name<- basename(x)
print(name)

chr=sapply(strsplit(name, ":"), "[", 2)
print(chr)

maf=0.01

#Here creating a pattern variable so that files for a given region

pattern=paste0("*",name,"*")
print(pattern)

#List all files for a give Region

mol_trait <- list.files(path = "/dir/proteomewide/eqtlgen/inflammation/data/regions/exposure", pattern=pattern, recursive = TRUE)
print(mol_trait)
pheno<- list.files(path = "/dir/proteomewide/eqtlgen/inflammation/data/regions/outcome", pattern=pattern, recursive = TRUE)
print(pheno)

#Assign to the identified files the full directory so that pwcoco can work 
mol_trait_dir<- paste0("/dir/proteomewide/eqtlgen/inflammation/data/regions/exposure/", mol_trait)
print(mol_trait_dir)
pheno_dir<- paste0("/dir/proteomewide/eqtlgen/inflammation/data/regions/outcome/", pheno)
print(pheno_dir)

setwd("/dir/pwcoco_update/pwcoco/build")

#Build the commands
bim_cmd<- paste("--bfile", paste0("/dir/REF_PANEL/refpanel_alspac/CHR", chr))
print(bim_cmd)
chr_cmd<- paste("--chr", chr)
print(chr_cmd)
maf_cmd<- paste("--maf", maf)
print(maf_cmd)

lapply(list(mol_trait_dir), function(x){
  
for(x in mol_trait_dir) {

sumstats1<- paste("--sum_stats1", x)
print(sumstats1)

lapply(list(pheno_dir), function(i){
  
  for(i in pheno_dir) {
print(i)

sumstats2<- paste("--sum_stats2", i)
print(sumstats2)

out_cmd<- paste("--out", paste0("/dir/proteomewide/eqtlgen/inflammation/results/colocalisation/eqtlgen"))
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


