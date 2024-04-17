#Two-sample MR analyses of METABRAIN eQTLs on psychiatric outcomes
#This script:
#Imports all Metabrain data
#Extracts genome-wide significant and independent immunological eQTLs
#Extracts from outcome, harmonises & performs two-sample MR
#Saves harmonised datasets as well as IVW & Wald ratio output for supplementary material 

#Load metabrain cortex data 

library(readr)
all_metabrain<- read_table("/path/METABRAIN/2021-07-23-cortex-EUR-80PCs-chr_all.txt")
all_metabrain_df<- as.data.frame(all_metabrain)

#Import the markers included in inflammation panels I & II so that markers can be extracted

base<- read.table("/path/panel/UKB_panel_base.txt")
colnames(base)<- c("olinkid", "UniProtID", "HGNC", "Ensemblid")

#Merge the two files so that immunological markers from metabrain can be extracted

metabrain_im<- merge(all_metabrain_df, base, by.x="GeneSymbol", by.y="HGNC")

library(dplyr)
metabrain_im<- metabrain_im%>%select(Gene, GeneChr, GenePos, GeneSymbol, SNPChr, SNPPos, SNPEffectAlleleFreq, MetaP, MetaBeta, MetaSE, RSID, A1, A2)
metabrain_im$N<- 2970
metabrain_im$ensemblid<- gsub("\\..*", "", metabrain_im$Gene)

files<- metabrain_im[which(metabrain_im$MetaP<=5e-08),]

files <- split(files, files$Gene)

#Load all packages

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(tidyr)
gwasvcf::set_bcftools("/hpcpath/bcftools")


mart.gene <- biomaRt::useMart(biomart = "ensembl",
                                  host="https://www.ensembl.org",
                                  dataset = "hsapiens_gene_ensembl")

#List all psychiatric GWAS data
								  
outcome_files <- list.files(path="/path/psychiatric_data", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE) 

#Loop to conduct mr analyses across all biomarkers and psychiatric outcomes

lapply(list(files), function(x){
  for(x in files) {
  
  tryCatch({

#Formating necessary for clumping and subsequent mr

exposure<- format_data(
  x,
  type = "exposure", snps = NULL, header = TRUE, phenotype_col = "GeneSymbol",
  snp_col = "RSID", beta_col = "MetaBeta", se_col = "MetaSE", effect_allele_col = "A1",
  other_allele_col = "A2", pval_col = "MetaP", samplesize_col = "N", id_col = "GeneSymbol",
  chr_col = "SNPChr", pos_col = "SNPPos", gene_col = "ensemblid", eaf_col= "SNPEffectAlleleFreq")

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

#Estimation of F statistic and cis/trans annotation (note METABRAIN offers cis only)

exposure_instruments<-  mrpipeline::calc_f_stat(exposure_clumped) %>%
  mrpipeline::cis_trans(values_col="exposure", chr_col="chr.exposure", 
  pos_col="pos.exposure", filter = "hgnc_symbol", build = "grch38") #metabrain is in hg38

exposure_instruments<-  exposure_instruments %>% 
  filter(!grepl('CHR', chromosome_name))

#Keeping variants passing the Fstatistic

exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

iterations= 1:8

#Looping over psychiatric outcomes

for (var in iterations){

try ({
outcome <- TwoSampleMR::read_outcome_data(snps = exposure_instruments$SNP, filename = outcome_files[[var]], 
                        sep = "\t", snp_col = "SNP", beta_col = "logOR", 
                        se_col = "SE", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", chr_col= "CHR", pos_col="BP", id_col="pheno", 
						ncase_col="NCAS", ncontrol_col="NCON", samplesize_col="NTOT")

#Harmonisation & keeping those that pass it

dat <- mrpipeline::harmonise(exposure_instruments, outcome, action = 2)
dat<- dat[which(dat$mr_keep=="TRUE"),]

#The following steps until saving the harmonised file are necessary to save all infomation considered necessary for the supplement

dat$marker<- paste0((dat$SNP), ":", (dat$chr.exposure), ":", (dat$exposure), ":", (dat$gene.exposure), ":", (dat$cistrans), ":", (dat$id.outcome))

dat_final<- dat%>%dplyr::select("marker","exposure", "gene.exposure", "SNP", "chr.exposure", "pos.exposure",
								"other_allele.exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", 
								"pval.exposure", "eaf.exposure", "samplesize.exposure", "id.exposure", "f.stat.exposure",
								"cistrans", "chr.outcome", "pos.outcome",
								"other_allele.outcome", "effect_allele.outcome", "beta.outcome", "se.outcome", 
								"pval.outcome", "eaf.outcome", "samplesize.outcome", "ncase.outcome", "ncontrol.outcome",
								"id.outcome", "remove", "palindromic", "ambiguous", "mr_keep")

setwd("/path/METABRAIN/MR/supplement")
if(file.exists("harmonised_metabrain.txt") == TRUE) {
  write.table(dat_final,"harmonised_metabrain.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(dat_final, "harmonised_metabrain.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#MR analyses

mr_res <- mrpipeline::do_mr(dat, f_cutoff = 10)

#Saving separately available IVW and individual Wald ratios

try ({mr_res_ivw<- mr_res[which(mr_res$method=="Inverse variance weighted"),]

mr_res_ivw_f<- mr_res_ivw %>% 
    mutate(snp = strsplit(as.character(snp), ",")) %>% 
    unnest(snp)
mr_res_ivw_f$snp<- gsub("\\s+", "", mr_res_ivw_f$snp)
	
mr_res_ivw_f<- as.data.frame(mr_res_ivw_f)

map<- dat%>%dplyr::select("marker", "SNP", "cistrans")

mr_res_ivw_ff<- merge(mr_res_ivw_f, map, by.x="snp", by.y="SNP")

mr_res_ivw_ff$marker<- paste0((mr_res_ivw_ff$marker), ":", "IVW")

mr_res_ivw_ff$notes<- ifelse(all(mr_res_ivw_ff$cistrans=="C"), "ALLCIS", 
						ifelse(all(mr_res_ivw_ff$cistrans=="T"), "ALLTRANS", "ATL1CIS"))
})


mr_res_wald<- mr_res[which(mr_res$method=="Wald ratio"),]

map<- dat%>%dplyr::select("marker", "SNP", "cistrans")
mr_res_wald_f<- merge(mr_res_wald, map, by.x="snp", by.y="SNP")
mr_res_wald_f$marker<- paste0((mr_res_wald_f$marker), ":", "WALD")

mr_res_wald_f$notes<- ifelse(nrow(mr_res_ivw)!=0, "PrioritiseIVW", "Thisisit")


setwd("/path/METABRAIN/MR/results")
if(file.exists("ivw_results_metabrain.txt") == TRUE) {
  write.table(mr_res_ivw_ff,"ivw_results_metabrain.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(mr_res_ivw_ff, "ivw_results_metabrain.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

if(file.exists("wald_results_metabrain.txt") == TRUE) {
  write.table(mr_res_wald_f,"wald_results_metabrain.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(mr_res_wald_f, "wald_results_metabrain.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
})
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


}})
