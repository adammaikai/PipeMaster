#!/usr/bin/env python
# -*- coding: utf-8 -*-
from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def bbduk(sample, bbduk_flag):
    ''' 
    
    input: 
        fastq
    output: 
        
    citation: 
        
    link: 
        
    parameters from parameters file: 

        '''
    if p.DNA["TUMOR_EXT"]:
        dna_tumor_sample = sample + p.DNA["TUMOR_EXT"]
        spawn_job(jobname = 'bbduk', SAMPLE = dna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.BBDUK["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.BBDUK["NODES"], ppn = p.BBDUK["CPU"], memory = p.BBDUK["MEMORY"], script = "/bbduk.sh", args_list = [sample, p.DNA["TUMOR_EXT"], p.BBDUK["RAW_DATA"], "dna", p.BBDUK["VERSION"], p.BBDUK["RESULTS"], p.BBDUK["REF_ADAPTER"], p.BBDUK["OPTIONS"], p.BBDUK["MEMORY"], str(p.BBDUK["CPU"])])
        job_status(jobname = 'bbduk', resultspath = p.BBDUK["RESULTS"] + "/" + dna_tumor_sample + "/stats/", SAMPLE = dna_tumor_sample, outputfilename = dna_tumor_sample + ".stats", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.DNA["NORMAL_EXT"]:
        dna_normal_sample = sample + p.DNA["NORMAL_EXT"]
        spawn_job(jobname = 'bbduk', SAMPLE = dna_normal_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.BBDUK["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.BBDUK["NODES"], ppn = p.BBDUK["CPU"], memory = p.BBDUK["MEMORY"], script = "/bbduk.sh", args_list = [sample, p.DNA["NORMAL_EXT"], p.BBDUK["RAW_DATA"], "dna", p.BBDUK["VERSION"], p.BBDUK["RESULTS"], p.BBDUK["REF_ADAPTER"], p.BBDUK["OPTIONS"], p.BBDUK["MEMORY"], str(p.BBDUK["CPU"])])
        job_status(jobname = 'bbduk', resultspath = p.BBDUK["RESULTS"] + "/" + dna_normal_sample + "/stats/", SAMPLE = dna_normal_sample, outputfilename = dna_normal_sample + ".stats", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.RNA["TUMOR_EXT"]:
        rna_tumor_sample = sample + p.RNA["TUMOR_EXT"]
        spawn_job(jobname = 'bbduk', SAMPLE = rna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.BBDUK["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.BBDUK["NODES"], ppn = p.BBDUK["CPU"], memory = p.BBDUK["MEMORY"], script = "/bbduk.sh", args_list = [sample, p.RNA["TUMOR_EXT"], p.BBDUK["RAW_DATA"], "rna", p.BBDUK["VERSION"], p.BBDUK["RESULTS"], p.BBDUK["REF_ADAPTER"], p.BBDUK["OPTIONS"], p.BBDUK["MEMORY"], str(p.BBDUK["CPU"])])
        job_status(jobname = 'bbduk', resultspath = p.BBDUK["RESULTS"] + "/" + rna_tumor_sample + "/stats/", SAMPLE = rna_tumor_sample, outputfilename = rna_tumor_sample + ".stats", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    bbduk(sample, bbduk_flag)
    sys.exit(0)
