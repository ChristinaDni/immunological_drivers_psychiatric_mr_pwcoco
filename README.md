# immunological_drivers_psychiatric_mr_pwcoco
Code accompanying manuscript: Insights into causal immunological biomarkers for major neuropsychiatric conditions through multi-tissue Mendelian randomization and genetic colocalisation

#Contents:

#1. Extraction of regions for UKB immunological blood pQTLs for subsequent colocalisation analyses

#2. Extraction of regions from the psychiatric outcomes and PWCOCO analyses with UKB pQTL regions

#3. Two-sample MR analyses of UKB immunological blood pQTLs on neuropsychiatric outcomes

#4. Bi-directional MR analyses of psychiatric phenotypes on UKB immunological protein levels

#5. Two-sample MR analyses of BrainQTL immunological brain pQTLs on neuropsychiatric outcomes (note that full data were not available and therefore colocalisation could not be performed)

#6. Extraction of regions for metabrain immunological brain eQTLs for subsequent colocalisation analyses

#7. Extraction of regions from the psychiatric outcomes and PWCOCO analyses with metabrain eQTL regions

#8. Two-sample MR analyses of metabrain immunological brain eQTLs on neuropsychiatric outcomes

#9. Extraction of regions for eqtlgen immunological blood eQTLs for subsequent colocalisation analyses

#10. Extraction of regions from the psychiatric outcomes and PWCOCO analyses with eqtlgen eQTL regions

#11. Two-sample MR analyses of eqtlgen immunological blood eQTLs on neuropsychiatric outcomes

#Files required:

#1. In order to extract the biomarkers of interest we need the protein annotation data file from the metadata section of the UKB PPP: https://www.synapse.org/#!Synapse:syn52364558 
#This file is necessary to extract the relevant data (immunological biomarkers) from UKB as well as rest of QTL datasets. 
#Two files with all the synids for UKB immunological biomarkers, Olink ID, UniProt ID, HGNC symbol & ensembl ids are provided in the data section of the repository.

#2. UKB GWASs for proteins included in Olink explore 3072 Inflammation Panels I & II available at http://ukb-ppp.gwas.eu

#3. Brain pQTL data for proteins included in Olink explore 3072 Inflammation Panels I & II available at https://www.synapse.org/#!Synapse:syn24172458

#4. EQTLGEN data for genes of proteins included in Olink explore 3072 Inflammation Panels I & II available at https://www.eqtlgen.org/phase1.html

#5. MetaBrain brain cortex data for genes of proteins included in Olink explore 3072 Inflammation Panels I & II available at https://www.metabrain.nl/

#6. Outcome data: GWASs of neuropsychiatric conditions
#ADHD, autism, anxiety, bipolar disorder, and schizophrenia can be accessed at: https://pgc.unc.edu/for-researchers/download-results/ 
#depression can be accessed at: https://ipsych.dk/en/research/downloads/
#Alzheimer’s disease can be accessed at: https://ctg.cncr.nl/software/summary_statistics

#Environment & software required to run the scripts:

#R version: R 4.1.2
#Plink: 1.90
#Bcftools: 1.9

#Analyses were carried out using the computational facilities of the Advanced Computing Research Centre of the University of Bristol (http://www.bris.ac.uk/acrc/). 
#Blood blood cell derived eQTL data were extracted and processed using the gwasvcf package version 1.0 in R (https://github.com/MRCIEU/gwasvcf). 
#Two-sample MR, Steiger filtering, and bi-directional MR analyses were conducted using functions from the TwoSampleMR R package version 0.5.6 (https://github.com/MRCIEU/TwoSampleMR) and the mrpipeline R package (https://github.com/jwr-git/mrpipeline).
#The PWCoCo algorithm was implemented using the Pair-Wise Conditional analysis and Colocalisation analysis package v1.0 (https://github.com/jwr-git/pwcoco)
