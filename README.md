# immunological_drivers_psychiatric_mr_pwcoco
Code accompanying manuscript: Immunological Drivers and Potential Novel Drug Targets for Major Psychiatric, Neurodevelopmental, and Neurodegenerative Conditions 
Contents:
1_ genetic instrument and region extraction from ukb data, pwcoco analyses
2_ two-sample mr analyses using ukb data
3_ genetic instrument and region extraction from eqtlgen data, pwcoco analyses
4_ two-sample mr analyses using eqtlgen data
5_ genetic instrument and region extraction from metabrain data, pwcoco analyses
6_ two-sample mr analyses using metabrain data
7_ two-sample mr analyses using brain pQTL data
8_ bi-directional two-sample mr analyses

~

Data availability:
Across all analyses published summary-level data were used and no patient identifiable information was included. UKB blood pQTL data can be accessed through the portal: http://ukb-ppp.gwas.eu. deCODE blood pQTL data can be accessed through the platform: https://www.decode.com/summarydata/. Brain pQTL data can be accessed through the https://adknowledgeportal.synapse.org/. Blood eQTL data can be accessed through https://www.eqtlgen.org/phase1.html. Brain cortex eQTL data can be accessed through the MetaBrain platform: https://www.metabrain.nl/. GWAS data on ADHD, autism, anxiety, bipolar disorder, and schizophrenia can be accessed at: https://pgc.unc.edu/for-researchers/download-results/. GWAS data on depression can be accessed at: https://ipsych.dk/en/research/downloads/. GWAS data on Alzheimer’s disease can be accessed at: https://ctg.cncr.nl/software/summary_statistics. 

~

Software:
Analyses were carried out using the computational facilities of the Advanced Computing Research Centre of the University of Bristol (http://www.bris.ac.uk/acrc/). Blood plasma pQTL data and blood cell derived eQTL data were extracted and processed using the gwasvcf package version 1.0 in R (https://github.com/MRCIEU/gwasvcf). The summary data from MetaBrain were lifted over  from GRCh38 to GRCh37 using the UCSC liftover tool to match the build of the rest of the data. Two-sample MR, Steiger filtering, and bi-directional MR analyses were conducted using functions from the TwoSampleMR R package version 0.5.6 (https://github.com/MRCIEU/TwoSampleMR) and the mrpipeline R package (https://github.com/jwr-git/mrpipeline). The PWCoCo algorithm was implemented using the Pair-Wise Conditional analysis and Colocalisation analysis package (https://github.com/jwr-git/pwcoco). 

~
