#!/usr/bin/env nextflow

params.reads = "$baseDir/../data/*{1,2}.fq.gz"

read_pairs_ch = Channel.fromFilePairs(params.reads)
                      .view()

process fastqc_raw_reads {

    input:
    tuple val(sample), file(reads) from read_pairs_ch

    script:
    """
    fastqc $reads
    """
}
