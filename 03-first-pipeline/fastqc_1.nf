#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*.fq.gz"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromPath( params.reads )
    .view()
    
process fastqc_raw_reads {

    input:
    path read from reads_ch 
    
    script:
    """
    fastqc ${read}
    """
}

