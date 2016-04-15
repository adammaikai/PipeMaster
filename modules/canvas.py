#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def canvas(sample, canvas_flag):
    '''Runs Illumina's copy number caller, Canvas,(med)exome data and estimates copy numbers."

        input:
        -(T,B)1-DNA_gatk_recal.bam
        output:
        .CNV.vcf.gz
        citation:
        parameters from parameters file:
        CNV_DIR:
        '''

    spawn_job(jobname = 'canvas', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"],
SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.CNV["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.CNV["NODES"], ppn = p.CNV["CPU"],
memory = p.CNV["MEMORY"], script = "/canvas.sh", args_list = [sample, p.CNV["CNV_DIR"], p.VARSCAN["ALIGNMENT_DIR"], p.ANNOTATION["VARIANT_DIR"]])
    job_status(jobname = 'canvas', resultspath = p.CNV["CNV_DIR"] + "/canvas/" + sample, SAMPLE = sample,  outputfilename = sample +
".VCF.gz", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    canvas(sample, canvas_flag)
    sys.exit(0)
