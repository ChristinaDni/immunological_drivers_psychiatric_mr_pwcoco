#!/bin/bash

cd /dir/proteomewide/ukb

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

cd /dir/proteomewide/ukb/scripts
Rscript 1_ukb_instrument_region_pwcoco.r