#!/usr/bin/env nextflow

params.reads = "$baseDir/../data/*.fq.gz"

/**
 * Quality control fastq
 */

read_ch = Channel
    .fromPath( params.reads )
    .view()

process fastqc_raw_reads {

    input:
    file read from read_ch
   
    script:
    """
    fastqc $read
    """
}