#!/bin/bash

cd /dir/proteomewide/metabrain

mkdir inflammation
cd inflammation
mkdir data results
cd data
mkdir instruments regions
cd regions
mkdir exposure outcome
cd ..
cd ..
cd results
mkdir mr colocalisation

module load gcc/12.3.0
module load languages/R/4.3.3

cd /dir/proteomewide/metabrain/scripts
Rscript 5_metabrain_instrument_region_pwcoco.r