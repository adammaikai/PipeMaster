#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load varscan/$2
module load R/$3
module load vcftools/$4
module load vcflib/$9

####INPUTS $1=sample $2=VARSCAN_VERSION $3=R_VERSION $4=VCFTOOLS_VERSION $5=RESULTS $6=ALIGNMENT_DIR $7=DNA_NORMAL_EXT $8=DNA_TUMOR_EXT
mkdir -p "$5"

normal=$(echo $7 | sed -e 's/-DNA//' | cut -c 2-)
tumor=$(echo $8 | sed -e 's/-DNA//' | cut -c 2-)

## Varscan2 call variants
VarScan.v$2.jar somatic $6/$1$7/$1$7.flt.pileup $6/$1$8/$1$8.flt.pileup $5/$1/$1\_$tumor\_$normal\_varscan_somatic --min-var-freq 0.01 --output-vcf 1

## Process output
VarScan.v$2.jar processSomatic \
$5/$1/$1\_$tumor\_$normal\_varscan_somatic.snp.vcf

VarScan.v$2.jar processSomatic \
$5/$1/$1\_$tumor\_$normal\_varscan_somatic.indel.vcf

vcfcombine $5/$1/$1\_$tumor\_$normal\_varscan_somatic.snp.Somatic.hc.vcf $5/$1/$1\_$tumor\_$normal\_varscan_somatic.indel.Somatic.hc.vcf > $5/$1/$1\_$tumor\_$normal\_varscan_somatic.vcf
bgzip -f $5/$1/$1\_$tumor\_$normal\_varscan_somatic.vcf
tabix -f -p vcf $5/$1/$1\_$tumor\_$normal\_varscan_somatic.vcf.gz

## Copynumber
VarScan.v2.3.9.jar \
copynumber \
$6/$1$7/$1$7.flt.pileup \
$6/$1$8/$1$8.flt.pileup \
$5/$1/$1\_$tumor\_$normal\.varscan_cnv

VarScan.v2.3.9.jar \
copyCaller \
$5/$1/$1\_$tumor\_$normal\.varscan_cnv.copynumber \
--output-file $5/$1/$1\_$tumor\_$normal\.varscan_cnv.copynumber.called \
--output-homdel-file $5/$1/$1\_$tumor\_$normal\.varscan_cnv.copynumber.homdel

#rm $5/$1/$1\_$tumor\_$normal\_varscan_somatic.snp*
#rm $5/$1/$1\_$tumor\_$normal\_varscan_somatic.indel*

exit 0
