#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"

/**
 * Quality control fastq
 */

println """\
       LIST OF PARAMETERS
================================
Reads            : $params.reads
Output-folder    : $params.outdir/
"""


read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)
        .view()                     

process fastqc_raw_reads {
    // Copies the output files into the published directory. Overwrite in case data has changed (other raw reads etc.)
    publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), file(reads) from read_pairs_ch

    //output:
    //file("fastqc_${sample}_logs")

    script:
    """
    mkdir -p $params.outdir/quality-control-$sample/
    fastqc --outdir $params.outdir/quality-control-$sample/ ${reads}
    """
}

