#!/usr/bin/env nextflow

params.reads = "$projectDir/../data/*.fq.gz"
// alternatively
// params.reads = "$launchDir/data/*.fq.gz"


/**
 * Quality control fastq
 */

Channel
    .fromPath( params.reads )
    .set{ reads_ch }
    

process fastqc_raw_reads {

    input:
    file read from reads_ch
   
    script:
    """
    fastqc $read
    """
}