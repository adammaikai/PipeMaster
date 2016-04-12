#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

## INPUTS : $1=sample  $2=CNV_DIR
mkdir -p $2/$1/cnvkit

docker pull anu9109/cnvkit

docker run -v /data:/home/data -i anu9109/cnvkit bash -c "cnvkit.py 
batch -p 0 /home/data/storage/patients/alignments/$1-B1-DNA/$1-B1-DNA_gatk_recal.bam -r 
/home/$2/cnvkit-medexome/pooledref_medexome.cnn --output-dir 
/home/$2/$1/cnvkit --scatter"

docker run -v /data:/home/data -i anu9109/cnvkit bash -c "cnvkit.py call 
/home/$2/$1/cnvkit/$1-T1-DNA_gatk_recal.cns -o 
/home/$2/$1/cnvkit/$1-T1-DNA_gatk_recal.call.cns"

docker run -v /data:/home/data -i anu9109/cnvkit bash -c "cnvkit.py 
export bed /home/$2/$1/cnvkit/$1-T1-DNA_gatk_recal.call.cns 
--show-neutral -o /home/$2/$1/cnvkit/$1.out.bed"

#module load R/$2
#Rscript nanostringCNV.R $3 $1

exit 0
