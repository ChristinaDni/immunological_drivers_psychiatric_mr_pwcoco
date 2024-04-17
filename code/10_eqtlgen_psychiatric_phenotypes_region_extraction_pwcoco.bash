#!/bin/bash

#Region extraction from psychiatric phenotypes and pwcoco analyses
#This script uses as an example ADHD
#The same code is applied to all psychiatric phenotypes used in the study considering that they have the same format- i.e. headers:
#SNP     CHR     BP      A1      A2      FRQ     logOR   SE      P       NTOT    NCAS    NCON    pheno   (INFO)


cd /path/EQTLGEN/PWCOCO/regions/outcome
mkdir adhd
cd /path/EQTLGEN/PWCOCO/results
mkdir adhd


BASEDIR="/path/EQTLGEN/PWCOCO/regions/exposure"
OUTDIR="/path/EQTLGEN/PWCOCO/regions/outcome/adhd"
OUTDAT="/path/psychiatric_data/C_adhd_2022.txt"
TOPHITSDIR="/path/EQTLGEN/PWCOCO/regions/topregions_eqtlgen.txt" #This file was created during region extraction for the eQTLs- it helps loop over pwcoco

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
set -- $(cat $TOPHITSDIR | tr -s ' ' | cut -f7) #extracting file name and chromosome information so that pwcoco can work in loop
for j in $(cat $TOPHITSDIR | tr -s ' ' | cut -f2) 
do

./pwcoco --bfile "/path/REF_PANEL/refpanel_alspac/CHR$j" \
--sum_stats1 "$BASEDIR/$1.txt" --sum_stats2 "$OUTDIR/$1.txt" \
--chr $j --maf 0.01 --out "/path/EQTLGEN/PWCOCO/results/adhd/adhd_EQTLGEN.out"
shift
done
