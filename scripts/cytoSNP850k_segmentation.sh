#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

## INPUTS: $1=sample  $2=R_version  $3=CNV_DIR

module load R/$2

# Segmentation
Rscript ~/.virtualenvs/pm/omics_pipe/PipeMaster/scripts/cytoSNP850k_segmentation.R $3/$1 $1

# Edit output file from segmentation step to remove `chr0` and `chrXY` values
for sample in $(ls /data/storage/patients/cnv/*/*seg | awk -F '/' '{print $6}'); do awk '{print $1, "chr"$2, $3, $4, $5, $6, $7, $8, $9, $10}' 
/data/storage/patients/cnv/$sample/$sample.cytoSNP850k.seg | grep -v XY | grep -v "chr0" > 
/data/storage/patients/cnv/$sample/$sample.cytoSNP850k.seg.updated; done

# Merge segments
perl ~/.virtualenvs/pm/omics_pipe/PipeMaster/scripts/mergeSegments.pl $3/$1/$1.cytoSNP850k.seg.updated --ref-arm-sizes 
/data/database/chrom_lengths_hg19_varscan2.txt --output-basename $3/$1/$1.cytoSNP850k

exit 0
