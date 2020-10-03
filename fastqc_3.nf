#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*0.01_{1,2}.fq.gz"
params.outdir = "$baseDir/results"

println """\
         NEXTFLOW COURSE
================================
Reads            : $params.reads
Output-folder    : $params.outdir
"""

read_pairs_ch = Channel.fromFilePairs(params.reads)
                      

process fastqc_raw_reads {
    publishDir "$params.outdir/quality-control", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), file(reads) from read_pairs_ch

    script:
    """
    mkdir -p $params.outdir/quality-control
    fastqc --outdir $params.outdir/quality-control ${reads}
    """
}