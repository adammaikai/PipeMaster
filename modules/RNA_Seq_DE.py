#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def RNA_Seq_DE(sample, RNA_Seq_DE_flag):
    '''Perform n of 1 differential expression from RNA Seq counts
        
        input:
        
        output:
        
        citation:
        
        link:
        
        parameters from parameters file:
        
        sample:
        
        '''
    
    if p.RNA["TUMOR_EXT"]:
        sample = sample + p.RNA["TUMOR_EXT"]
    spawn_job(jobname = 'RNA_Seq_DE', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"], SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = "240:00:00", queue = p.OMICSPIPE["QUEUE"], nodes = 1, ppn = 32, memory = "58gb", script = "/RNA_Seq_DE.sh", args_list = [sample, p.GATK["R_VERSION"])
    job_status(jobname = 'RNA_Seq_DE', resultspath = p.ANNOTATION["VARIANT_DIR"], SAMPLE = sample,  outputfilename = sample + "", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
   RNA_Seq_DE(sample, RNA_Seq_DE_flag)
   sys.exit(0)

