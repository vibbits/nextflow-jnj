#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*0.01_1.fq.gz"

println "reads: $params.reads"

/**
 * Quality control fastq
 */

read_ch = Channel.fromPath( params.reads )

process fastqc_raw_reads {

    input:
    file read from read_ch
   
    script:
    """
    fastqc $read
    """
}