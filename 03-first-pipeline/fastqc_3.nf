#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*{1,2}.fq.gz"

/**
 * Quality control fastq
 */


println """\
      LIST OF PARAMETERS
================================
Reads         : $params.reads
"""

read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true) // Alternatively: .ifEmpty { 'No read pairs found' } 
        .view()

process fastqc_raw_reads {

    input:
    tuple val(sample), file(reads) from read_pairs_ch

    script:
    """
    fastqc ${reads}
    """
}
