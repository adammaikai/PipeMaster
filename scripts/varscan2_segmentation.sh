#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

## INPUTS: $1=sample  $2=R_version  $3=CNV_DIR

module load R/$2

# Segmentation
Rscript ~/.virtualenvs/pm/omics_pipe/PipeMaster/scripts/varscan2_segmentation.R $3/$1 $1

# Merge segments
perl ~/.virtualenvs/pm/omics_pipe/PipeMaster/scripts/mergeSegments.pl $3/$1/$1.varscan_cnv.copynumber.seg --ref-arm-sizes 
/data/database/chrom_lengths_hg19_varscan2.txt --output-basename $3/$1/$1.varscan_cnv

exit 0
