#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load mojo/$4
cpu=$(grep -c "processor" /proc/cpuinfo)

mkdir -p "$3/$1"

####INPUTS $1=sample $2=RAW_DATA $3=RESULTS $4=MOJO_VERSION $5=CONFIG $6=CPU $7=MEMORY

MOJO \
--config $5 \
--sample_name $1 \
--output_dir $3/$1 \
--fq1 $2/$1/$1_1.fastq.gz \
--fq2 $2/$1/$1_2.fastq.gz \
--cores $6 \
--mem 50

cat $3/$1/$1/$1\.fusions | awk -F '\t' '{print $1}' | sed -e 's/_/--/g' > $3/$1/$1/$1\.fusions.genes



exit 0