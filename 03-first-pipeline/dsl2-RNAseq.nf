#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Similar to DSL1, the input data is defined in the beginning.
params.reads   = "$launchDir/data/*{1,2}.fq.gz"
params.outdir  = "$launchDir/results"
params.threads = 2

params.slidingwindow   = "SLIDINGWINDOW:4:15"
params.avgqual         = "AVGQUAL:30"

params.dirgenome   = "$launchDir/data"
params.genome      = "$launchDir/data/Drosophila_melanogaster.BDGP6.dna.fa"
params.gtf         = "$launchDir/data/Drosophila_melanogaster.BDGP6.85.sample.gtf"
params.lengthreads = 98


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
================================
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

//dirgenome = file(params.dirgenome)
genome = file(params.genome)
gtf = file(params.gtf)

// Process trimmomatic
process trimmomatic {
    publishDir "$params.outdir/trimmed-reads", mode: 'copy'

    // Same input as fastqc on raw reads, comes from the same channel. 
    input:
    tuple val(sample), file(reads) 

    output:
    tuple val(sample), file("${sample}{1,2}_P.fq"), emit: paired_fq
    //tuple val(sample), file("${sample}_1U.fq"), file("${sample}_2U.fq"), emit: unpaired_fq

    script:
    """
    mkdir -p $params.outdir/trimmed-reads/
    trimmomatic PE -threads $params.threads ${reads[0]} ${reads[1]} ${sample}1_P.fq ${sample}1_U.fq ${sample}2_P.fq ${sample}2_U.fq $params.slidingwindow $params.avgqual 
    """
}

include { QC as fastqc_raw; QC as fastqc_trim } from "${launchDir}/modules/fastqc" 
include { IDX; MAP } from "${launchDir}/modules/star"

// Running a workflow with the defined processes here.  
workflow {
	read_pairs_ch.view()
	fastqc_raw(read_pairs_ch) 
  paired_fq = trimmomatic(read_pairs_ch)
  fastqc_trim(paired_fq)
  IDX(genome, gtf)
  MAP(paired_fq, gtf)
}
