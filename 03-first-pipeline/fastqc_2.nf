#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*{1,2}.fq.gz"

/**
 * Quality control fastq
 */

read_pairs_ch = Channel.fromFilePairs(params.reads)
                      .view()

process fastqc_raw_reads {

    input:
    tuple val(sample), file(reads) from read_pairs_ch

    script:
    """
    fastqc ${reads}
    """
}
