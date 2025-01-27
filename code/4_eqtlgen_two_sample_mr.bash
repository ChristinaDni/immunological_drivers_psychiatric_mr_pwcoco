#!/bin/bash

module load gcc/12.3.0
module load languages/R/4.3.3
cd /dir/proteomewide/eqtlgen/scripts
Rscript 4_eqtlgen_two_sample_mr.r