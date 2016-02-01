#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load R/$3

####INPUTS: $1=sample $2=VARIANT_DIR $3=R_VERSION 

Rscript ~/.virtualenvs/op2/omics_pipe/omics_pipe/scripts/somaticAnnotation.R $2/$1/$1 $2/$1/$1_varscan_somatic.vcf.gz $2/$1/$1_mutect.filt.vcf.gz

exit 0