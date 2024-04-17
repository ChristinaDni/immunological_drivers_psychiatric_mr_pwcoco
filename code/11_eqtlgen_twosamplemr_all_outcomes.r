#Two-sample MR analyses of EQTLGEN eQTLs on psychiatric outcomes
#This script:
#Extracts all EQTLGEN data based on the immunological marker panel
#Extracts genome-wide significant and independent eQTLs
#Extracts from outcome, harmonises & performs two-sample MR
#Saves harmonised datasets as well as IVW & Wald ratio output for supplementary material

#Load all packages

library(mrpipeline)
library(dplyr)
library(biomaRt)
library(TwoSampleMR)
library(tidyr)
gwasvcf::set_bcftools("/hpcpath/bin/bcftools")

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

outcome_files <- list.files(path="/path/psychiatric_data", 
pattern="*.txt$", full.names=TRUE, recursive=FALSE) 

								 
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

#Keeping variants that pass the f statistic

exposure_instruments<- exposure_instruments[which(exposure_instruments$f.stat.exposure>=10),]

iterations= 1:8

for (var in iterations){

try ({
outcome <- TwoSampleMR::read_outcome_data(snps = exposure_instruments$SNP, filename = outcome_files[[var]], 
                        sep = "\t", snp_col = "SNP", beta_col = "logOR", 
                        se_col = "SE", effect_allele_col = "A1", 
                        other_allele_col = "A2", pval_col = "P", eaf_col="FRQ", chr_col= "CHR", pos_col="BP", id_col="pheno", 
						ncase_col="NCAS", ncontrol_col="NCON", samplesize_col="NTOT")

#Harmonisation and keeping those that pass it

dat <- mrpipeline::harmonise(exposure_instruments, outcome, action = 2)
dat<- dat[which(dat$mr_keep=="TRUE"),]

#Merging with panel file here to retain infomation for the supplement

dat<- merge(dat, base, by.x="exposure", by.y="Ensemblid")

dat$marker<- paste0((dat$SNP), ":", (dat$chr.exposure), ":", (dat$HGNC), ":", (dat$exposure), ":", (dat$cistrans), ":", (dat$id.outcome))

dat_final<- dat%>%dplyr::select("marker","exposure", "HGNC", "SNP", "chr.exposure", "pos.exposure",
								"other_allele.exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", 
								"pval.exposure", "eaf.exposure", "samplesize.exposure", "id.exposure", "f.stat.exposure",
								"cistrans", "chr.outcome", "pos.outcome",
								"other_allele.outcome", "effect_allele.outcome", "beta.outcome", "se.outcome", 
								"pval.outcome", "eaf.outcome", "samplesize.outcome", "ncase.outcome", "ncontrol.outcome",
								"id.outcome", "remove", "palindromic", "ambiguous", "mr_keep")

setwd("/path/EQTLGEN/MR/supplement")
if(file.exists("harmonised_eqtlgen.txt") == TRUE) {
  write.table(dat_final,"harmonised_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(dat_final, "harmonised_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

#MR analysis

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


setwd("/path/EQTLGEN/MR/results")
if(file.exists("ivw_results_eqtlgen.txt") == TRUE) {
  write.table(mr_res_ivw_ff,"ivw_results_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(mr_res_ivw_ff, "ivw_results_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}

if(file.exists("wald_results_eqtlgen.txt") == TRUE) {
  write.table(mr_res_wald_f,"wald_results_eqtlgen.txt", sep="\t", quote=F, row.names=F, col.names=F, append=TRUE)
} else { 
  write.table(mr_res_wald_f, "wald_results_eqtlgen.txt", sep="\t", quote=F, row.names=F, append=FALSE)
}
})
}
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


}})
