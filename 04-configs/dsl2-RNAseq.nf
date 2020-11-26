#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2


println """\
      LIST OF PARAMETERS
================================
            GENERAL
Reads            : $params.reads
Results-folder   : $params.outdir/
Threads          : $params.threads
================================
          TRIMMOMATIC
Sliding window   : $params.slidingwindow
Average quality  : $params.avgqual
================================
             STAR
Reference genome : $params.dirgenome
Genome directory : $params.genome 
GTF-file         : $params.gtf
Length-reads     : $params.lengthreads
SAindexNbases    : $params.genomeSAindexNbases
================================
"""


// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

genome = file(params.genome)
gtf = file(params.gtf)

index = Channel
    .fromPath(params.indexpath, checkIfExists:true)

include { QC as fastqc_raw; QC as fastqc_trim } from "${launchDir}/modules/fastqc" //addParams(OUTPUT: fastqcOutputFolder)
include { IDX; MAP } from "${launchDir}/modules/star"
include { TRIM } from "${launchDir}/modules/trimmomatic"

// Running a workflow with the defined processes here.  
workflow {
	read_pairs_ch.view()
  index.view()
	fastqc_raw(read_pairs_ch) 
	//TRIM(read_pairs_ch)
  //trim_paired_fq = TRIM(read_pairs_ch)
	//fastqc_trim(trim_paired_fq)
	//index = IDX(genome, gtf)
  //trim_paired_fq.view()
	IDX(genome, gtf)
  //MAP(read_pairs_ch, index, gtf)
}
