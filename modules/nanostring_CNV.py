#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def nanostringCNV(sample, nanostringCNV_flag):
    '''Performs normalization on Nanostring output data based on invariant regions and estimates copy numbers for ~ 90 genes
        
        input:
        -nanostringCNV-T.RCC
        output:
        .nanostringCNV.out
        citation:
        parameters from parameters file:
        CNV_DIR:
        
        R_VERSION:
        '''
    
    spawn_job(jobname = 'nanostringCNV', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], 
SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.CNV["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.CNV["NODES"], ppn = p.CNV["CPU"], 
memory = p.CNV["MEMORY"], script = "/nanostringCNV.sh", args_list = [sample, p.CNV["R_VERSION"], p.CNV["CNV_DIR"]])
    job_status(jobname = 'nanostringCNV', resultspath = p.CNV["CNV_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample + 
".nanostringCNV.out", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    nanostringCNV(sample, nanostringCNV_flag)
    sys.exit(0)
