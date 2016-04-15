#!/bin/bash

## Script to run illumina's copy number caller, Canvas on Pipemaster.
## Note: The .bed file for Canvas must contain the header [Regions] in the line preceeding the column names-Chromosome, Start, End, Name, Length and Orientation.


#### NOTE: update the truseq .bed file to the medexome .bed file ####

set -x

# source modules for current shell
source $MODULESHOME/init/bash

## INPUTS : $1=sample  $2=CNV_DIR    $3=ALIGNMENT_DIR    $4=VARIANT_DIR
mkdir -p $2/canvas/$1

# Run the module for tumor-normal samples
mono /opt/software/canvas/1.3.9/canvas-1.3.9_x64/Canvas.exe Tumor-normal-enrichment -b $3/$1-T1-DNA/$1-T1-DNA_gatk_recal.bam --normal-bam=$3/$1-B1-DNA/$1-B1-DNA_gatk_recal.bam --b-allele-vcf=$4/$1-B1-DNA/$1-B1-DNA_merged.vcf.gz --reference=/data/database/canvas/kmer.fa --manifest=/data/database/canvas/truseq.exome.targeted.regions.hg19.canvas.bed -g /data/database/canvas/ -n $1 -f /data/database/canva/filter13.bed -o $2/canvas/$1

# Rename final VCF output file with patient sample ID in the filename
mv $2/canvas/$1/CNV.vcf.gz $2/canvas/$1/$1.CNV.vcf.gz 

