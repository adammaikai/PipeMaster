#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

## INPUTS : $1=sample  $2=R_version  $3=CNV_DIR

module load R/$2

Rscript nanostringCNV.R $3 $1

exit 0
