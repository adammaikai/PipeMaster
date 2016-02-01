#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def WES_preprocess(sample, WES_preprocess_flag):
    '''Preprocessing steps for whole exome sequencing.
        
        input:
        .fastq.gz
        output:
        _realigned_sorted.bam
        citation:
        McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20:1297-303.
        link:
        http://www.broadinstitute.org/gatk/
        parameters from parameters file:
        
        TEMP_DIR:
        
        GENOME:
        
        BWA_VERSION:
        
        ABRA_VERSION:
        
        SAMTOOLS_VERSION:
        
        SAMBAMBA_VERSION:
        
        SAMBLASTER_VERSION:
    
        SAMTOOLS_VERSION:
        '''
    if p.DNA["TUMOR_EXT"]:
        dna_tumor_sample = sample + p.DNA["TUMOR_EXT"]
        spawn_job(jobname = 'WES_preprocess', SAMPLE = dna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.PREPROCESS["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.PREPROCESS["NODES"], ppn = p.PREPROCESS["CPU"], memory = p.PREPROCESS["MEMORY"], script = "/WES_preprocess.sh", args_list = [dna_tumor_sample, p.OMICSPIPE["TEMP_DIR"], p.PREPROCESS["GENOME"], p.PREPROCESS["BWA_INDEX"], p.PREPROCESS["BWA_VERSION"], p.PREPROCESS["SAMTOOLS_VERSION"], p.PREPROCESS["SAMBAMBA_VERSION"], p.PREPROCESS["SAMBLASTER_VERSION"], p.PREPROCESS["ABRA_VERSION"], p.CAPTURE_KIT_BED, p.PREPROCESS["ALIGNMENT_DIR"], p.PREPROCESS["FASTQ_PATH"], p.OMICSPIPE["LOG_PATH"], sample + p.DNA["TUMOR_EXT"]])
        job_status(jobname = 'WES_preprocess', resultspath = p.PREPROCESS["ALIGNMENT_DIR"] + "/" + dna_tumor_sample, SAMPLE = dna_tumor_sample,  outputfilename = dna_tumor_sample + "_realigned_sorted.bam", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.DNA["NORMAL_EXT"]:
        dna_normal_sample = sample + p.DNA["NORMAL_EXT"]
        spawn_job(jobname = 'WES_preprocess', SAMPLE = dna_normal_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.PREPROCESS["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.PREPROCESS["NODES"], ppn = p.PREPROCESS["CPU"], memory = p.PREPROCESS["MEMORY"], script = "/WES_preprocess.sh", args_list = [dna_normal_sample, p.OMICSPIPE["TEMP_DIR"], p.PREPROCESS["GENOME"], p.PREPROCESS["BWA_INDEX"], p.PREPROCESS["BWA_VERSION"], p.PREPROCESS["SAMTOOLS_VERSION"], p.PREPROCESS["SAMBAMBA_VERSION"], p.PREPROCESS["SAMBLASTER_VERSION"], p.PREPROCESS["ABRA_VERSION"], p.CAPTURE_KIT_BED, p.PREPROCESS["ALIGNMENT_DIR"], p.PREPROCESS["FASTQ_PATH"], p.OMICSPIPE["LOG_PATH"], sample + p.DNA["NORMAL_EXT"]])
        job_status(jobname = 'WES_preprocess', resultspath = p.PREPROCESS["ALIGNMENT_DIR"] + "/" + dna_normal_sample, SAMPLE = dna_normal_sample,  outputfilename = dna_normal_sample + "_realigned_sorted.bam", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    WES_preprocess(sample, WES_preprocess_flag)
    sys.exit(0)
