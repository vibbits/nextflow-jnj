#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*{1,2}.fq.gz"

/**
 * Quality control fastq
 */

read_pairs_ch = Channel.fromFilePairs(params.reads)
                      .view()

process fastqc_raw_reads {

    input:
    //  first element is the grouping key of the matching pair and the second element is the list of files 
    tuple val(sample), path(reads) from read_pairs_ch

    script:
    """
    fastqc ${reads}
    """
}
