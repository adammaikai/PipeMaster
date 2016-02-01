#!/bin/bash
set -x

# source modules for current shell
source $MODULESHOME/init/bash

module load bwa/$5
module load samtools/$6
module load sambamba/$7
module load samblaster/$8
module load gatk/$9
cpu=$(grep -c "processor" /proc/cpuinfo)

####INPUTS: $1=sample $2=TEMP_DIR $3=GENOME $4=BWA_INDEX $5=BWA_VERSION $6=SAMTOOLS_VERSION $7=SAMBAMBA_VERSION $8=SAMBLASTER_VERSION $9=GATK_VERSION $10=CAPTURE_KIT_BED $11=ALIGNMENT_DIR $12=FASTQ_PATH $13=LOG_PATH $14=DBSNP $15=MILLS $16=G1000

mkdir -p "${11}/$1"

RGR="@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:$1"

## align to hg19
bwa mem -M -t $cpu -R $RGR $4 ${12}/$1/$1\_1.fastq.gz ${12}/$1/$1\_2.fastq.gz \
| samblaster --splitterFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 16G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin) \
--discordantFile >(samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 16G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin) \
| samtools view -S -u /dev/stdin \
| sambamba sort -t $cpu -m 16G --tmpdir $2 -o $2/$1\_sorted.bam /dev/stdin

## index
samtools index $2/$1\_sorted.bam

## mark duplicates
if [ -f $2/$1\_sorted.bam ]; then
	sambamba markdup -t $cpu --overflow-list-size 1000000 --hash-table-size 1000000 $2/$1\_sorted.bam ${11}/$1/$1\_dedup.bam
fi

## index the bam file
if [ -f ${11}/$1/$1\_dedup.bam ]; then
	samtools index  ${11}/$1/$1\_dedup.bam
    rm $2/$1\_sorted*
fi

## GATK Realignment
java -Xmx28g -jar `which GenomeAnalysisTK.jar` -T RealignerTargetCreator -nt $cpu -R $3 -I ${11}/$1/$1\_dedup.bam -known ${15} -known ${16} -o $2/$1\.intervals

java -Xmx28g -jar `which GenomeAnalysisTK.jar` -T IndelRealigner -R $3 -I ${11}/$1/$1\_dedup.bam -targetIntervals $2/$1\.intervals -known ${15} -known ${16} -o ${11}/$1/$1\_realigned.bam

if [ -f ${11}/$1/$1\_realigned.bam ]; then
    rm ${11}/$1/$1\_dedup.ba*
fi

## BQSR
java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
-T BaseRecalibrator -nct $cpu \
-R $3 \
-I ${11}/$1/$1\_realigned.bam  \
-knownSites ${14} \
-knownSites ${15} \
-knownSites ${16} \
-o $2/$1\_recal_data.table

java -Xmx28g -jar `which GenomeAnalysisTK.jar` \
-T PrintReads -nct $cpu \
-R $3 \
-I ${11}/$1/$1\_realigned.bam \
-BQSR $2/$1\_recal_data.table \
-o ${11}/$1/$1\_gatk_recal.bam

## sort and index again
if [ -f ${11}/$1/$1\_gatk_recal.bam ]; then
    rm ${11}/$1/$1\_realigned.ba*
    rm $2/$1\_recal_data.table
    rm $2/$1.intervals
fi


exit 0




