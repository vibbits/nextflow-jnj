#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Similar to DSL1, the input data is defined in the beginning.
params.reads = "$launchDir/data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"
params.threads = 2
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"


println """\
      LIST OF PARAMETERS
================================
            GENERAL
Reads            : $params.reads
Output-folder    : $params.outdir/

          TRIMMOMATIC
Threads          : $params.threads
Sliding window   : $params.slidingwindow
Avg quality      : $params.avgqual
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

// Definition of a process, notice the absence of the 'from channel'.
// A process being defined, does not mean it's invoked (see workflow)
process fastqc {
  publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
    
  input:
  tuple val(sample), path(reads)

  script:
  """
  mkdir -p $params.outdir/quality-control-$sample
  fastqc ${reads}
  """
}

// Process trimmomatic
process trimmomatic {
    publishDir "$params.outdir/trimmed-reads", mode: 'copy'

    // Same input as fastqc on raw reads, comes from the same channel. 
    input:
    tuple val(sample), path(reads) 

    output:
    tuple val("${sample}"), path("${sample}*_P.fq"), emit: paired_fq
    tuple val("${sample}"), path("${sample}*_U.fq"), emit: unpaired_fq

    script:
    """
    mkdir -p $params.outdir/trimmed-reads/
    trimmomatic PE -threads $params.threads ${reads[0]} ${reads[1]} ${sample}1_P.fq ${sample}1_U.fq ${sample}2_P.fq ${sample}2_U.fq $params.slidingwindow $params.avgqual 
    """
}

// Running a workflow with the defined processes here.  
workflow {
	read_pairs_ch.view()
	fastqc(read_pairs_ch) 
  paired_fq = trimmomatic(read_pairs_ch)
  trimmomatic.out.paired_fq.view()
  //fastqc(trimmomatic.out.paired_fq) // This will raise an error. Do you remember why?
}
