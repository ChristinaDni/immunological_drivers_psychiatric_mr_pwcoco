#EQTLGEN eQTL data region extraction for subsequent colocalisation analyses
#This script:
#Extracts all EQTLGEN data based on the immunological marker panel
#Extracts genome-wide significant and independent eQTLs
#Saves them for supplementary material
#Prioritises for reqion extraction cis eQTLs 
#If no cis eQTLs exist then prioritises trans eQTLs with the smallest p-value

#Load all packages

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
gwasvcf::set_bcftools("/hpcpath/bcftools")

mart.gene <- biomaRt::useMart(biomart = "ensembl",
                                  host="https://www.ensembl.org",
                                  dataset = "hsapiens_gene_ensembl")
								  
#Load the panel to aid extraction

base<-read.table("/path/panel/UKB_panel_base.txt", header=F)
colnames(base)<- c("olinkid", "UniProtID", "HGNC", "Ensemblid")

#Note that the eqtlgen data were available in the mrc ieu repository
#All eqtlgen files have been named using the prefix eqtl-a- followed by the Ensemblid
#For this reason creating a column in the panel file using the path that the data are stored and the prefix so that the data can be extracted

base$path_name <- sub("^", "/hpcpath/public/eqtl-a-", base$Ensemblid)

#Extract the files of interest based on this column

files<- list.files(path=base$path_name, 
pattern="*.vcf.gz$", full.names=TRUE, recursive=FALSE)

#For the following loop functions available in the gwasvcf package are being used
								 
lapply(list(files), function(x){
  for(x in files) {

tryCatch({

#Here genome-wide and independent eQTLs are being extracted 

exposure_instruments<- read_exposure(x, pval = 5e-8, plink="/hpcpath/plink-1.90/plink", 
bfile= "/path/REF_PANEL/EUR", clump_r2=0.01)

exposure_instruments$exposure<- gsub('eqtl-a-','', exposure_instruments$exposure)

#Estimation of F statistic and cis/trans annotation

exposure_instruments<-  mrpipeline::calc_f_stat(exposure_instruments) %>%
  mrpipeline::cis_trans(values_col="exposure", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "ensembl_gene_id", build = "grch37") 

exposure_instruments<-  exposure_instruments %>% 
  filter(!grepl('CHR', chromosome_name))

exposure_instruments$basename<- basename(exposure_instruments$file.exposure)

#Here merging with the base file in order to retain information to be included in the supplement

exposure_instruments<- merge(exposure_instruments, base, by.x="exposure", by.y="Ensemblid")

#Keeping instruments that pass the f statistic

exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

#Keeping information that the supplementary file needs to have

exposure_instruments_final<- exposure_instruments%>%dplyr::select("exposure", "HGNC", "SNP", "chr.exposure", "pos.exposure",
								"other_allele.exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", 
								"pval.exposure", "eaf.exposure", "samplesize.exposure", "f.stat.exposure",
								"cistrans")

setwd("/path/EQTLGEN/MR/instruments")
if(file.exists("instruments_eqtlgen.txt") == TRUE) {
  write.table(exposure_instruments_final,"instruments_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(exposure_instruments_final, "instruments_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#Until here we extract instruments

#From here onwards we extract regions
#Prioritising extraction and therefore subsequent colocalisation analyses for biomarkers with at least one cis eQTL

if(any(exposure_instruments_final$cistrans=="C")){
instruments<- exposure_instruments_final[which(exposure_instruments_final$cistrans=="C"),]

map<- instruments %>% split(x=instruments, f=instruments$SNP)

lapply(list(map), function(i){
  
  for(i in map) {

region<- i[c("SNP","chr.exposure","pos.exposure", "exposure", "HGNC", "cistrans")]

#Here we make an indicator with information that is useful for pwcoco loop to work

region$TopSNP<- paste0((i$SNP),
                     ":",
                     (i$chr.exposure), ":", (i$HGNC),
                     ":", (i$exposure), ":", (i$cistrans))

setwd("/path/EQTLGEN/PWCOCO/regions/")
if(file.exists("topregions_eqtlgen.txt") == TRUE) {
  write.table(region,"topregions_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#This indicator here is useful for naming the files- again for pwcoco loop to work and for retaining important information

snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$HGNC, ":", region$exposure, ":", region$cistrans)
snp_list<- as.data.frame(snp_list)

#Using functions of the gwasvcf package to extract regions

chrpos <- paste0(as.character(region$chr.exposure),
                     ":",
                     (as.numeric(region$pos.exposure) - 500000),
                     "-",
                     as.numeric(region$pos.exposure) + 500000)

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
    ) 

#Transforming the p-value as p-values have been originally stored as -log10 
	
variants_region_final$p <- 10^(-variants_region_final$p)

variants_region_final<- as.data.frame(variants_region_final)

print(nrow(variants_region_final))

setwd("/path/EQTLGEN/PWCOCO/regions/exposure")

snp<- snp_list$snp_list

file_out_snp<- paste0(snp, ".txt")
for(o in 1:length(snp)) {

write.table(variants_region_final, file_out_snp[o], sep="\t", quote=F, row.names=F, col.names=T)}

}
})
}

#Here creating second condition in case the biomarker has trans eQTLs only
#The process that follows is identical with above with the exception that here the trans variant with the smallest p-value is prioritised

else{
instruments<- exposure_instruments_final[which.min(exposure_instruments_final$pval.exposure),]
region<- instruments[c("SNP","chr.exposure","pos.exposure", "exposure", "HGNC", "cistrans")]

region$TopSNP<- paste0((instruments$SNP),
                     ":",
                     (instruments$chr.exposure), ":", (instruments$HGNC),
                     ":", (instruments$exposure), ":", (instruments$cistrans))

setwd("/path/EQTLGEN/PWCOCO/regions/")
if(file.exists("topregions_eqtlgen.txt") == TRUE) {
  write.table(region,"topregions_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(region, "topregions_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

snp_list<- paste0(region$SNP,":",region$chr.exposure, ":", region$HGNC, ":", region$exposure, ":", region$cistrans)
snp_list<- as.data.frame(snp_list)

chrpos <- paste0(as.character(region$chr.exposure),
                     ":",
                     (as.numeric(region$pos.exposure) - 500000),
                     "-",
                     as.numeric(region$pos.exposure) + 500000)

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
    ) 
	
variants_region_final$p <- 10^(-variants_region_final$p)

variants_region_final<- as.data.frame(variants_region_final)

print(nrow(variants_region_final))

setwd("/path/EQTLGEN/PWCOCO/regions/exposure")

snp<- snp_list$snp_list
file_out_snp<- paste0(snp, ".txt")

write.table(variants_region_final, file_out_snp, sep="\t", quote=F, row.names=F, col.names=T)}


}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}})

