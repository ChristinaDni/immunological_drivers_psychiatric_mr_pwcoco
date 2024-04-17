#!/bin/bash

#Region extraction from psychiatric phenotypes and pwcoco analyses
#This script uses as an example ADHD
#The same code is applied to all psychiatric phenotypes used in the study considering that they have the same format- i.e. headers:
#SNP     CHR     BP      A1      A2      FRQ     logOR   SE      P       NTOT    NCAS    NCON    pheno   (INFO)

cd /path/UKB/PWCOCO/regions/outcome
mkdir adhd
cd /path/UKB/PWCOCO/results
mkdir adhd


BASEDIR="/path/UKB/PWCOCO/regions/exposure"
OUTDIR="/path/UKB/PWCOCO/regions/outcome/adhd"
OUTDAT="/path/psychiatric_data/C_adhd_2022.txt"
TOPHITSDIR="/path/UKB/PWCOCO/regions/topregions_ukb.txt" #This file was created during region extraction for the pQTLs- it helps loop over pwcoco

cd $BASEDIR

#Region extraction from the psychiatric GWAS- ADHD in this case

for f in *.txt; do
awk '{ print $1 }'  $BASEDIR/$f > $OUTDIR/tmp_"$f"
grep -wFf $OUTDIR/tmp_"$f" $OUTDAT > $OUTDIR/tmp2_"$f"
awk '{ print $1,$4,$5,$6,$7,$8,$9,$10,$11 }' $OUTDIR/tmp2_"$f" > $OUTDIR/"$f"
rm $OUTDIR/tmp_"$f" $OUTDIR/tmp2_"$f"
done


module load /hpcpath/gcc/9.3.0
cd /path/pwcoco-master/build

set -f
IFS='
'
set -- $(cat $TOPHITSDIR | tr -s ' ' | cut -f8) #extracting file name and chromosome information so that pwcoco can work in loop
for j in $(cat $TOPHITSDIR | tr -s ' ' | cut -f2) 
do

./pwcoco --bfile "/path/REF_PANEL/refpanel_alspac/CHR$j" \
--sum_stats1 "$BASEDIR/$1.txt" --sum_stats2 "$OUTDIR/$1.txt" \
--chr $j --maf 0.01 --out "/path/UKB/PWCOCO/results/adhd/adhd_ukb.out"
shift
done
