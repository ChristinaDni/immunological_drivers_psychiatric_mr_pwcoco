#Two-sample MR analyses of brain pQTLs on psychiatric outcomes
#This script:
#Imports all brain pQTL data- please note that regions are not provided so only two sample MR analyses for these data
#Extracts available genome-wide significant pQTLs
#Extracts from outcome, harmonises & performs two-sample MR
#Saves harmonised datasets as well as IVW & Wald ratio output for supplementary material

#load necessary packages
library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
gwasvcf::set_bcftools("/hpcpath/bcftools")

#Load the markers included in inflammation panels I & II so that they can be extracted from BrainQTL

base<- read.table("/path/panel/UKB_panel_base.txt")
colnames(base)<- c("olinkid", "UniProtID", "HGNC", "Ensemblid")

#Load all brainQTL data

brainqtl<- read.table("/path/BRAINQTL/data.txt", header=T, sep="\t")

#Merge to extract the immunological proteins
files<- merge(brainqtl, base, by="UniProtID")

files<- files%>%dplyr::select(UniProtID, SNP.x, CHR, POS, Gene_name, Protein_name, other_allele.exposure, effect_allele.exposure,
						eaf.exposure, beta.exposure, se.exposure, samplesize.exposure, pval.exposure)
						
files <- split(files, files$UniProtID)

outcome_files <- list.files(path="/path/psychiatric_data", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE) 

#Loop to analyse all brain pQTLs on all psychiatric outcomes

lapply(list(files), function(x){
  for(x in files) {
  
  tryCatch({

#Formating for mr analyses to run

exposure<- format_data(
  x,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  phenotype_col = "UniProtID",
  snp_col = "SNP.x",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  gene_col = "Gene_name",
  chr_col = "CHR",
  pos_col = "POS",
  eaf_col="eaf.exposure",
  id_col="Protein_name"
)

#Estimating f statistic 

exposure_instruments<-  mrpipeline::calc_f_stat(exposure) %>%
  mrpipeline::cis_trans(values_col="gene.exposure", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "hgnc_symbol", build = "grch37") 
  
  exposure_instruments<-  exposure_instruments %>% 
  filter(!grepl('CHR', chromosome_name))

#Keeping variants with an F>10
 
  exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

setwd("/path/BRAINQTL/instruments")
if(file.exists("instruments_brainqtl.txt") == TRUE) {
  write.table(exposure_instruments,"instruments_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(exposure_instruments, "instruments_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

iterations= 1:8

#Looping over all psychiatric outcomes

for (var in iterations){
try ({
outcome<- read_outcome_data(snps = exposure_instruments$SNP, filename = outcome_files[[var]], 
                        sep = "\t", snp_col = "SNP", beta_col = "logOR", chr_col= "CHR", pos_col="BP",
                        se_col = "SE", eaf_col = "FRQ", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", id_col="pheno", samplesize_col="NTOT", phenotype_col="pheno")

#Harmonisation & keeping only those that pass it

dat <- mrpipeline::harmonise(exposure_instruments, outcome, action = 2)
dat<- dat[which(dat$mr_keep=="TRUE"),]

setwd("/path/BRAINQTL/supplement")
if(file.exists("harmonised_brainqtl.txt") == TRUE) {
  write.table(dat,"harmonised_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(dat, "harmonised_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#MR analyses

mr_res <- mrpipeline::do_mr(dat, f_cutoff = 10)

#Saving separately ivw and wald ratios- although in this case we had single instruments for each available biomarker

try ({mr_res_ivw<- mr_res[which(mr_res$method=="Inverse variance weighted"),]
mr_res_ivw$marker<- paste0((mr_res_ivw$snp), ":", (mr_res_ivw$exposure), ":", (mr_res_ivw$id.exposure))})

mr_res_wald<- mr_res[which(mr_res$method=="Wald ratio"),]
mr_res_wald$marker<- paste0((mr_res_wald$snp), ":", (mr_res_wald$exposure), ":", (mr_res_wald$id.exposure))

#Keeping all necessary information to be included in the supplement

map<- dat
map$marker<- paste0((map$SNP), ":", (map$exposure), ":", (map$id.exposure))

mr_wald<- merge(mr_res_wald, map, by="marker")
mr_wald<- mr_wald[,c("exposure.x", "outcome.x", "id.exposure.x", "id.outcome.x", "method", 
			"nsnp", "snp", "b", "se", "pval", "lo_ci", "up_ci", "or", "or_lci95", "or_uci95", 
			"egger_intercept", "se.egger", "pval.egger", "snp_r2.exposure", "snp_r2.outcome", 
			"correct_causal_direction", "steiger_pval", "steigerflag", "marker", "f.stat.exposure", 
			"chr.exposure", "pos.exposure", "gene.exposure", "ensembl_gene_id")]

colnames(mr_wald)[1]  <- "exposure"
colnames(mr_wald)[2]  <- "outcome"
colnames(mr_wald)[3]  <- "id.exposure"
colnames(mr_wald)[4]  <- "id.outcome"
colnames(mr_wald)[10]  <- "pval"


setwd("/path/BRAINQTL/results")
if(file.exists("ivw_results_brainqtl.txt") == TRUE) {
  write.table(mr_res_ivw,"ivw_results_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(mr_res_ivw, "ivw_results_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

if(file.exists("wald_results_brainqtl.txt") == TRUE) {
  write.table(mr_wald,"wald_results_brainqtl.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(mr_wald, "wald_results_brainqtl.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

})

}


}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})



}
})
