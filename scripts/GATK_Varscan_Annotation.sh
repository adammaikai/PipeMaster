#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load R/$4
module load vcflib/$5
module load vcftools/$6

####INPUTS: $1=sample $2=VARIANT_DIR $3=GENOME $4=R_VERSION $5=VCFLIB_VERSION $6=VCFTOOLS_VERSION

if [ -f $2/$1/$1\_vqsr.vcf.gz ]; then
	vcfintersect -r $3 -i $2/$1/$1\_vqsr.vcf.gz $2/$1/$1\_varscan.4.1.vcf.gz > $2/$1/$1\_merged.vcf
else
	vcfintersect -r $3 -i $2/$1/$1\_gatk.vcf.gz $2/$1/$1\_varscan.4.1.vcf.gz > $2/$1/$1\_merged.vcf
fi

bgzip -f $2/$1/$1_merged.vcf
tabix -f -p vcf $2/$1/$1_merged.vcf.gz

Rscript ~/.virtualenvs/op2/omics_pipe/omics_pipe/scripts/annotateVariantsFromVcf.R $2/$1/$1\_merged.vcf.gz FALSE FALSE

exit 0
