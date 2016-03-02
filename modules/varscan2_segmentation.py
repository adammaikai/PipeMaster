#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def varscan2_segmentation(sample, varscan2_segmentation_flag):
    '''Runs DNAcopy from Bioconductor on Varscan2 copynumber output to apply circular binary segmentation.
        
        input:
        .varscan_cnv.copynumber
        output:
        .varscan_cnv.events.tsv
        citation:
        parameters from parameters file:
        CNV_DIR:
        
        R_VERSION:
        '''
    
    spawn_job(jobname = 'varscan2_segmentation', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], 
SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.CNV["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.CNV["NODES"], ppn = p.CNV["CPU"], 
memory = p.CNV["MEMORY"], script = "/varscan2_segmentation.sh", args_list = [sample, p.CNV["R_VERSION"], p.CNV["CNV_DIR"]])
    job_status(jobname = 'varscan2_segmentation', resultspath = p.CNV["CNV_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + 
".varscan_cnv.events.tsv", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    varscan2_segmentation(sample, varscan2_segmentation_flag)
    sys.exit(0)
