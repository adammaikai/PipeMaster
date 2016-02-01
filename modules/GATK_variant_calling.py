#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_variant_calling(sample, GATK_variant_calling_flag):
    '''Calling variants with GATK for whole exome sequencing.
        
        input:
        _WES_realigned_sorted.bam
        output:
        _vqsr.vcf
        citation:
        
        link:
        
        parameters from parameters file:
        
        TEMP_DIR:
        
        GENOME:
        
        GATK_VERSION:
        
        R_VERSION:
        
        CAPTURE_KIT_BED:
        
        ALIGNMENT_DIR:
        
        RESULTS:
        
        OMNI:
        
        HAPMAP:
        
        DBSNP:
        
        MILLS:
        
        G1000:
        
        KG:
        '''
    if p.DNA["TUMOR_EXT"]:
        dna_tumor_sample = sample + p.DNA["TUMOR_EXT"]
        spawn_job(jobname = 'GATK_variant_calling', SAMPLE = dna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.GATK["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.GATK["NODES"], ppn = p.GATK["CPU"], memory = p.GATK["MEMORY"], script = "/GATK_variant_calling.sh", args_list = [dna_tumor_sample, p.OMICSPIPE["TEMP_DIR"], p.GATK["GENOME"], p.GATK["VERSION"], p.GATK["R_VERSION"], p.CAPTURE_KIT_BED, p.GATK["ALIGNMENT_DIR"], p.GATK["RESULTS"], p.GATK["OMNI"], p.GATK["HAPMAP"], p.GATK["DBSNP"], p.GATK["MILLS"], p.GATK["G1000"], p.GATK["KG"]])
        job_status(jobname = 'GATK_variant_calling', resultspath = p.GATK["RESULTS"] + "/" + dna_tumor_sample, SAMPLE = dna_tumor_sample,  outputfilename = dna_tumor_sample + "_gatk.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.DNA["NORMAL_EXT"]:
        dna_normal_sample = sample + p.DNA["NORMAL_EXT"]
        spawn_job(jobname = 'GATK_variant_calling', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.GATK["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.GATK["NODES"], ppn = p.GATK["CPU"], memory = p.GATK["MEMORY"], script = "/GATK_variant_calling.sh", args_list = [dna_normal_sample, p.OMICSPIPE["TEMP_DIR"], p.GATK["GENOME"], p.GATK["VERSION"], p.GATK["R_VERSION"], p.CAPTURE_KIT_BED, p.GATK["ALIGNMENT_DIR"], p.GATK["RESULTS"], p.GATK["OMNI"], p.GATK["HAPMAP"], p.GATK["DBSNP"], p.GATK["MILLS"], p.GATK["G1000"], p.GATK["KG"]])
        job_status(jobname = 'GATK_variant_calling', resultspath = p.GATK["RESULTS"] + "/" + dna_normal_sample, SAMPLE = dna_normal_sample,  outputfilename = dna_normal_sample + "_gatk.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    GATK_variant_calling(sample, GATK_variant_calling_flag)
    sys.exit(0)
