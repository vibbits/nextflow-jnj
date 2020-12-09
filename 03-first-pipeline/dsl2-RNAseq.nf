#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Define project parameters needed for running the pipeline
params.reads = "$launchDir/data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"
params.threads = 2
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"
params.dirgenome = "$launchDir/data"
params.genome = "$launchDir/data/Drosophila_melanogaster.BDGP6.dna.fa"
params.gtf = "$launchDir/data/Drosophila_melanogaster.BDGP6.85.sample.gtf"
params.genomeSAindexNbases = 10
params.lengthreads = 98
params.indexpath = "$launchDir/data/index/"


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

include { QC as fastqc_raw; QC as fastqc_trim } from "${launchDir}/modules/fastqc" //addParams(OUTPUT: fastqcOutputFolder)
include { IDX; MAP } from "${launchDir}/modules/star"
include { TRIM } from "${launchDir}/modules/trimmomatic"
include { MULTIQC } from "${launchDir}/modules/multiqc"

// Running a workflow with the defined processes here.  
workflow {
	// QC on raw reads
  fastqc_raw_out = fastqc_raw(read_pairs_ch) 
	
  // Trimming & QC
  trim_fq = TRIM(read_pairs_ch)
	fastqc_trim_out = fastqc_trim(trim_fq)
	
  // Mapping
  index_dir = IDX(genome, gtf)
  MAP(trim_fq, index_dir, gtf)
  
  // Multi QC on all results
  MULTIQC(fastqc_raw_out.mix(fastqc_trim_out).collect())
}
