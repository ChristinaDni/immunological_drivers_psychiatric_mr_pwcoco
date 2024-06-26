#UKB pQTL data region extraction for subsequent colocalisation analyses
#This script:
#Imports all GWAS inflammation panel datasets from UKB PPP
#Extracts genome-wide significant and independent pQTLs
#Saves them for supplementary material
#Prioritises for reqion extraction cis pQTLs 
#If no cis pQTLs exist then prioritises trans pQTLs with the smallest p-value

setwd("/path/UKBGWAS_local_directory")

#Load all immunological protein GWAS from UKB

files<- list.files(path="/path/UKBGWAS_local_directory/", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE)

#Load necessary packages

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
gwasvcf::set_bcftools("/hpcpath/bcftools")

mart.gene <- biomaRt::useMart(biomart = "ensembl",
                                  host="https://www.ensembl.org",
                                  dataset = "hsapiens_gene_ensembl")

#Loop to extract regions for all available pQTLs

lapply(list(files), function(x){
  for(x in files) {
  
  tryCatch({

file<- vroom::vroom(x)
file_df<- as.data.frame(file)

#Extraction of genome-wide significant pQTLs

file_i<- file_df[which(file_df$P<=5e-08),]

#Formating necessary to use clumping function

exposure<- format_data(
  file_i,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "olink_target_fullname",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  id_col = "HGNC.symbol",
  chr_col = "CHROM",
  pos_col = "POS38",
  gene_col = "ensembl_id",
  eaf_col= "A1FREQ"
)

exposure$rsid<- exposure$SNP
exposure$pval<- exposure$pval.exposure

#Clumping

exposure_clumped<- ieugwasr::ld_clump(
  dat = exposure,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 0.99,
  pop = "EUR",
  bfile = "/path/REF_PANEL/EUR",
  plink_bin ="/hpcpath/plink-1.90/plink")
  
#Anotation of cis and trans variants and estimation of f statistic
  
exposure_instruments<-  mrpipeline::calc_f_stat(exposure_clumped) %>%
  mrpipeline::cis_trans(values_col="id.exposure", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "hgnc_symbol", build = "grch38") #UKB selected column is in hg38

exposure_instruments<-  exposure_instruments %>% 
  filter(!grepl('CHR', chromosome_name))

#Keeping variants that pass F statistic

exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

exposure_instruments_final<- merge(exposure_instruments, file_df, by.x="SNP", by.y="rsid")

#Addition here- removing duplicates generated by mr pipeline

exposure_instruments_final$toremove<- ifelse(exposure_instruments_final$gene.exposure!=exposure_instruments_final$ensembl_gene_id, "Remove", "Keep")
exposure_instruments_final<- exposure_instruments_final[which(exposure_instruments_final$toremove!="Remove"),]

#Keeping information that the supplementary file needs to have

exposure_instruments_final<- exposure_instruments_final%>%dplyr::select("id.exposure", "gene.exposure", "SNP", "chr.exposure", "pos.exposure",
								"other_allele.exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", 
								"pval.exposure", "eaf.exposure", "samplesize.exposure", "f.stat.exposure",
								"cistrans", "exposure", "olink_assay_ukb", "UniProt")

setwd("/path/UKB/MR/instruments")
if(file.exists("instruments_ukb.txt") == TRUE) {
  write.table(exposure_instruments_final,"instruments_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(exposure_instruments_final, "instruments_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#Until here we extract instruments

#From here onwards we extract regions

#Prioritising extraction and therefore subsequent colocalisation analyses for biomarkers with at least one cis pQTL

if(any(exposure_instruments_final$cistrans=="C")){
instruments<- exposure_instruments_final[which(exposure_instruments_final$cistrans=="C"),]

map<- instruments %>% split(x=instruments, f=instruments$SNP)
  
  lapply(list(map), function(i){
  
  for(i in map) {

region<- i[c("SNP","chr.exposure","pos.exposure", "id.exposure", "cistrans", "olink_assay_ukb", "gene.exposure")]

#Here we make an indicator with information that is useful for pwcoco loop to work

region$TopSNP<- paste0((i$SNP),
                     ":",
                     (i$chr.exposure), ":", (i$id.exposure),
                     ":", (i$gene.exposure), ":", (i$cistrans),
					 ":",
                     (i$olink_assay_ukb))


setwd("/path/UKB/PWCOCO/regions")
if(file.exists("topregions_ukb.txt") == TRUE) {
  write.table(region,"topregions_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#This indicator here is useful for naming the files- again for pwcoco loop to work and for retaining important information

snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$id.exposure, ":", region$gene.exposure, ":", region$cistrans, ":", region$olink_assay_ukb)
snp_list<- as.data.frame(snp_list)

#Region extraction

all_chr<- file_df[which(file_df$CHROM==region$chr.exposure),]
all_pos<- all_chr[which(all_chr$POS38>=region$pos.exposure-500000 & all_chr$POS38<=region$pos.exposure+500000),]

variants_region <- all_pos %>% 
    dplyr::select("rsid", "A1", "A2", "A1FREQ", "BETA", "SE", "P", "N") 
    
variants_region_final<- variants_region %>% dplyr::rename(
      SNP = rsid,
      A1 = A1,
      A2 = A2,
      freq = A1FREQ,
      b = BETA,
      se = SE,
      p = P,
      N = N,
    ) 
	
variants_region_final<- as.data.frame(variants_region_final)

print(nrow(variants_region_final))

setwd("/path/UKB/PWCOCO/regions/exposure")

snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region_final, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}
  }
})}

#Here creating second condition in case the biomarker has trans pQTLs only
#The process that follows is identical with above with the exception that here the trans variant with the smallest p-value is prioritised

else{
instruments<- exposure_instruments_final[which.min(exposure_instruments_final$pval.exposure),]
region<- instruments[c("SNP","chr.exposure","pos.exposure", "id.exposure", "cistrans", "olink_assay_ukb", "gene.exposure")]

region$TopSNP<- paste0((instruments$SNP),
                     ":",
                     (instruments$chr.exposure), ":", (instruments$id.exposure),
                     ":", (instruments$gene.exposure), ":", (instruments$cistrans),
					 ":",
                     (instruments$olink_assay_ukb))


setwd("/path/UKB/PWCOCO/regions")
if(file.exists("topregions_ukb.txt") == TRUE) {
  write.table(region,"topregions_ukb.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_ukb.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$id.exposure, ":", region$gene.exposure, ":", region$cistrans, ":", region$olink_assay_ukb)
snp_list<- as.data.frame(snp_list)

all_chr<- file_df[which(file_df$CHROM==region$chr.exposure),]
all_pos<- all_chr[which(all_chr$POS38>=region$pos.exposure-500000 & all_chr$POS38<=region$pos.exposure+500000),]

variants_region <- all_pos %>% 
    dplyr::select("rsid", "A1", "A2", "A1FREQ", "BETA", "SE", "P", "N") 
    
variants_region_final<- variants_region %>% dplyr::rename(
      SNP = rsid,
      A1 = A1,
      A2 = A2,
      freq = A1FREQ,
      b = BETA,
      se = SE,
      p = P,
      N = N,
    ) 
	
variants_region_final<- as.data.frame(variants_region_final)

print(nrow(variants_region_final))

setwd("/path/UKB/PWCOCO/regions/exposure")
snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")

write.table(variants_region_final, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)}


  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }})




