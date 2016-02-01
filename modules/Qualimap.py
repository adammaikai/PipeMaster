#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def Qualimap(sample, Qualimap_flag):
    '''Runs Qualimap on a processed .bam file.
        
        input:
        _gatk_recal.bam
        output:
    
        citation:
        parameters from parameters file:

        ALIGNMENT_DIR:
        
        RESULTS:
        
        QUALIMAP_VERSION:

        QUALIMAP_OPTIONS:
        '''
    if p.DNA["TUMOR_EXT"]:
        dna_tumor_sample = sample + p.DNA["TUMOR_EXT"]
        spawn_job(jobname = 'Qualimap', SAMPLE =  dna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.QUALIMAP["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.QUALIMAP["NODES"], ppn = p.QUALIMAP["CPU"], memory = p.QUALIMAP["MEMORY"], script = "/Qualimap.sh", args_list = [ dna_tumor_sample, p.QUALIMAP["ALIGNMENT_DIR"], p.QUALIMAP["RESULTS"], p.QUALIMAP["VERSION"], "DNA", p.QUALIMAP["DNA_OPTIONS"], p.CAPTURE_KIT_BED])
        job_status(jobname = 'Qualimap', resultspath = p.QUALIMAP["RESULTS"] + "/" + dna_tumor_sample, SAMPLE = dna_tumor_sample,  outputfilename = "genome_results.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.DNA["NORMAL_EXT"]:
        dna_normal_sample = sample + p.DNA["NORMAL_EXT"]
        spawn_job(jobname = 'Qualimap', SAMPLE = dna_normal_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.QUALIMAP["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.QUALIMAP["NODES"], ppn = p.QUALIMAP["CPU"], memory = p.QUALIMAP["MEMORY"], script = "/Qualimap.sh", args_list = [dna_normal_sample, p.QUALIMAP["ALIGNMENT_DIR"], p.QUALIMAP["RESULTS"], p.QUALIMAP["VERSION"], "DNA", p.QUALIMAP["DNA_OPTIONS"], p.CAPTURE_KIT_BED])
        job_status(jobname = 'Qualimap', resultspath = p.QUALIMAP["RESULTS"] + "/" + dna_normal_sample, SAMPLE = dna_normal_sample,  outputfilename = "genome_results.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    if p.RNA["TUMOR_EXT"]:
        rna_tumor_sample = sample + p.RNA["TUMOR_EXT"]
        spawn_job(jobname = 'Qualimap', SAMPLE = rna_tumor_sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.QUALIMAP["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.QUALIMAP["NODES"], ppn = p.QUALIMAP["CPU"], memory = p.QUALIMAP["MEMORY"], script = "/Qualimap.sh", args_list = [rna_tumor_sample, p.QUALIMAP["ALIGNMENT_DIR"], p.QUALIMAP["RESULTS"], p.QUALIMAP["VERSION"], "RNA", p.QUALIMAP["RNA_OPTIONS"]])
        job_status(jobname = 'Qualimap', resultspath = p.QUALIMAP["RESULTS"] + "/" +rna_tumor_sample, SAMPLE = rna_tumor_sample,  outputfilename = "genome_results.txt", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    Qualimap(sample, Qualimap_flag)
    sys.exit(0)
