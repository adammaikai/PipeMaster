#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)

def fastqc(sample, fastqc_flag):
    '''QC check of raw .fastq files using FASTQC.
    
        input: 
            .fastq file
        output: 
            folder and zipped folder containing html, txt and image files
        citation: 
            Babraham Bioinformatics
        link: 
            http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
        parameters from parameters file: 
            RAW_DATA_DIR:
            
            QC_PATH:
            
            FASTQC_VERSION:
            
            COMPRESSION:
            ''' 
    print "sample name is: ", sample

    if p.DNA["TUMOR_EXT"]:
        dna_tumor_sample = sample + p.DNA["TUMOR_EXT"]       
        spawn_job(jobname = 'fastqc', SAMPLE = dna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.FASTQC["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.FASTQC["NODES"], ppn = p.FASTQC["CPU"], memory = p.FASTQC["MEMORY"], script = "/fastqc_drmaa.sh", args_list = [dna_tumor_sample, p.FASTQC["PATH"],p.FASTQC['RESULTS'], p.FASTQC["VERSION"], p.FASTQC["COMPRESSION"]])
        job_status(jobname = 'fastqc', resultspath = p.FASTQC["RESULTS"] + "/" + dna_tumor_sample, SAMPLE = dna_tumor_sample, outputfilename = dna_tumor_sample + "_2" + "_fastqc.html", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.DNA["NORMAL_EXT"]:
        dna_normal_sample = sample + p.DNA["NORMAL_EXT"]     
        spawn_job(jobname = 'fastqc', SAMPLE = dna_normal_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.FASTQC["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.FASTQC["NODES"], ppn = p.FASTQC["CPU"], memory = p.FASTQC["MEMORY"], script = "/fastqc_drmaa.sh", args_list = [dna_normal_sample, p.FASTQC["PATH"],p.FASTQC['RESULTS'], p.FASTQC["VERSION"], p.FASTQC["COMPRESSION"]])
        job_status(jobname = 'fastqc', resultspath = p.FASTQC["RESULTS"] + "/" + dna_normal_sample, SAMPLE = dna_normal_sample, outputfilename = dna_normal_sample + "_2" + "_fastqc.html", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])  
    if p.RNA["TUMOR_EXT"]:
        rna_tumor_sample = sample + p.RNA["TUMOR_EXT"]      
        spawn_job(jobname = 'fastqc', SAMPLE = rna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.FASTQC["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.FASTQC["NODES"], ppn = p.FASTQC["CPU"], memory = p.FASTQC["MEMORY"], script = "/fastqc_drmaa.sh", args_list = [rna_tumor_sample, p.FASTQC["PATH"],p.FASTQC['RESULTS'], p.FASTQC["VERSION"], p.FASTQC["COMPRESSION"]])
        job_status(jobname = 'fastqc', resultspath = p.FASTQC["RESULTS"] + "/" + rna_tumor_sample, SAMPLE = rna_tumor_sample, outputfilename = rna_tumor_sample + "_2" + "_fastqc.html", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    fastqc(sample, fastqc_flag)
    sys.exit(0)
