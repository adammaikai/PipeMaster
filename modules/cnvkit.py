#!/usr/bin/env python

from omics_pipe.parameters.default_parameters import default_parameters
from omics_pipe.utils import *
p = Bunch(default_parameters)


def cnvkit(sample, cnvkit_flag):
    '''Runs the CNVkit pipeline from a docker container on (med)exome data and estimates copy numbers, produces scatter plots of logR ratio data across the exome"

        input:
        -T1-DNA_gatk_recal.bam
        output:
        .out.bed, -T1-DNA_gatk_recal-scatter.pdf
	intermediate output:
	-T1-DNA_gatk_recal.targetcoverage.cnn, -T1-DNA_gatk_recal.antitargetcoverage.cnn, -T1-DNA_gatk_recal.cnr, -T1-DNA_gatk_recal.cns, -T1-DNA_gatk_recal.call.cns
        citation:
        parameters from parameters file:
        CNV_DIR:
        '''

    spawn_job(jobname = 'cnvkit', SAMPLE = sample, LOG_PATH = p.OMICSPIPE["LOG_PATH"], RESULTS_EMAIL = p.OMICSPIPE["EMAIL"],
SCHEDULER = p.OMICSPIPE["SCHEDULER"], walltime = p.CNV["WALLTIME"], queue = p.OMICSPIPE["QUEUE"], nodes = p.CNV["NODES"], ppn = p.CNV["CPU"],
memory = p.CNV["MEMORY"], script = "/cnvkit.sh", args_list = [sample, p.CNV["CNV_DIR"]])
    job_status(jobname = 'cnvkit', resultspath = p.CNV["CNV_DIR"] + "/" + sample, SAMPLE = sample,  outputfilename = sample +
".out.bed", FLAG_PATH = p.OMICSPIPE["FLAG_PATH"])
    return

if __name__ == '__main__':
    cnvkit(sample, cnvkit_flag)
    sys.exit(0)
