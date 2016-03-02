#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def GATK_Varscan_Annotation(sample, GATK_Varscan_Annotation_flag):
    '''Intersect GATK and Varscan VCF files and annotating with MyVariant.info.
        
        input:
        _vqsr.vcf.gz, _varscan.4.1.vcf.gz
        output:
        _merged.annnotated.vcf.gz
        citation:
        
        link:
        
        parameters from parameters file:
        
        sample:
        
        VARIANT_DIR:
        '''
    if p.DNA["NORMAL_EXT"]:
        dna_normal_sample = sample + p.DNA["NORMAL_EXT"]
        spawn_job(jobname = 'GATK_Varscan_Annotation', SAMPLE = dna_normal_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.ANNOTATION["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.ANNOTATION["NODES"], ppn = p.ANNOTATION["CPU"], memory = p.ANNOTATION["MEMORY"], script = "/GATK_Varscan_Annotation.sh", args_list = [dna_normal_sample, p.ANNOTATION["VARIANT_DIR"], p.ANNOTATION["GENOME"], p.VARSCAN["R_VERSION"], p.VARSCAN["VCFLIB_VERSION"], p.VARSCAN["VCFTOOLS_VERSION"]])
        job_status(jobname = 'GATK_Varscan_Annotation', resultspath = p.ANNOTATION["VARIANT_DIR"] + "/" + dna_normal_sample, SAMPLE = dna_normal_sample,  outputfilename = dna_normal_sample + "_merged.annotated.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.DNA["TUMOR_EXT"]:
        dna_tumor_sample = sample + p.DNA["TUMOR_EXT"]
        spawn_job(jobname = 'GATK_Varscan_Annotation', SAMPLE = dna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.ANNOTATION["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.ANNOTATION["NODES"], ppn = p.ANNOTATION["CPU"], memory = p.ANNOTATION["MEMORY"], script = "/GATK_Varscan_Annotation.sh", args_list = [dna_tumor_sample, p.ANNOTATION["VARIANT_DIR"], p.ANNOTATION["GENOME"], p.VARSCAN["R_VERSION"], p.VARSCAN["VCFLIB_VERSION"], p.VARSCAN["VCFTOOLS_VERSION"]])
        job_status(jobname = 'GATK_Varscan_Annotation', resultspath = p.ANNOTATION["VARIANT_DIR"] + "/" + dna_tumor_sample, SAMPLE = dna_tumor_sample,  outputfilename = dna_tumor_sample + "_merged.annnotated.vcf.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
   GATK_Varscan_Annotation(sample, GATK_Varscan_Annotation_flag)
   sys.exit(0)

