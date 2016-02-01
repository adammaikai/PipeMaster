#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Loads specified python addons module for cutadapt
module load bbmap/$5
cpu=$(grep -c "processor" /proc/cpuinfo)

#Make assemblies output director if it doesn't exist
mkdir -p $6/$1$2

#Runs bbduk with $1=PATIENT, $2=EXT, $3=RAW_DATA_DIR, $4=DNA/RNA, $5=BBDUK_VERSION, $6=RESULTS, $7=REF_ADAPTOR, $8=BBDUK_PARAMETERS, $9=MEMORY, $10=CPU

mem=$(expr $(echo $9 | sed 's/gb//') - 4)

bbduk.sh \
	in=$3/$1/$4/$1$2_1.fastq.gz \
	in2=$3/$1/$4/$1$2_2.fastq.gz \
	out=$6/$1$2/$1$2_1.fastq.gz \
	out2=$6/$1$2/$1$2_2.fastq.gz \
	ref=$7 \
	stats=$6/$1$2/stats/$1$2.stats \
	bhist=$6/$1$2/stats/$1$2.bhist \
	qhist=$6/$1$2/stats/$1$2.qhist \
	aqhist=$6/$1$2/stats/$1$2.aqhist \
	lhist=$6/$1$2/stats/$1$2.lhist \
	gchist=$6/$1$2/stats/$1$2.gchist \
	gcbins=auto \
	threads=${10} \
	$8 \
	-Xmx${mem}g

exit 0
